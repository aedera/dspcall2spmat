#!/usr/bin/env python3
import os
import argparse
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

import gzip
import scipy
import scipy.sparse
from scipy.sparse import coo_array

import multiprocessing as mp

import numpy as np
import h5py

import math

import sys
import pdb

def read_chrlen_file(fname, chrname):
    """
    Read file specifying the chromosome lengths
    """
    with open(fname) as f:
        for line in f:
            if line[0] == '#':
                continue # skip header
            cname, length = line.rstrip().split('\t')
            length = int(length)
            if cname == chrname:
                return length

    # halt execution if chromosome is not found
    print('Chromosome name not found.')
    sys.exit(1)

    return -1

class ForkedPdb(pdb.Pdb):
    """A Pdb subclass that may be used
    from a forked multiprocessing child

    """
    def interaction(self, *args, **kwargs):
        _stdin = sys.stdin
        try:
            sys.stdin = open('/dev/stdin')
            pdb.Pdb.interaction(self, *args, **kwargs)
        finally:
            sys.stdin = _stdin

def seek_start_position(f):
    """.

    Move file cursor to the beginning of a SAM line containing the first mate
    of a read pair.

    Args
    ----
    f: Python file object

    Return
    ------
    new_pos: file cursor pointing to the beginning of a read pair.
    """
    # seek beginning of a new line
    char = f.read(1)
    if char != '\n':
        f.readline()

    return f.tell()

def parse_chunk(targetchr, fname, pid, chunk_sz, target_ctx='all'):
    """
    Convert a bed file containing ONT reads into entries of a COO sparse matrix.
    """
    f = open(fname, 'r')
    f.seek((pid * chunk_sz), os.SEEK_SET)
    fpos = f.tell()

    new_fpos = fpos if pid == 0 else seek_start_position(f)
    f.seek(new_fpos, os.SEEK_SET)

    #read2row = {}
    read2data = {}
    read_bytes = abs(fpos - new_fpos)
    while True:
        # line should present the following format
        # Chr3    25396495        +       25396495        cae576ce-2a7b-4699-87ba-84fc9720cf9e    t       1.0     0.0     0       AACAC
        line = f.readline()
        if line == '': # end of file
            break

        line = line.rstrip().split('\t')
        if len(line) < 10:
            print('Error: Unknown format')
            break

        chrname, pos, strand, _, read_id, _, _, dscore, call, kmer = line
        tri = kmer[2:]
        if tri[:2] == 'CG':
            ctx = 'CG'
        elif tri[0] == 'C' and tri[-1] == 'G':
            ctx = 'CHG'
        elif tri[0] == 'C' and tri[-1] != 'G':
            ctx = 'CHH'
        else:
            ctx = 'UNK'
        # ForkedPdb().set_trace()

        strand = +1 if strand == '+' else -1 # fw/+, rv/-
        # dscore is the deep-signal-plant score assigns for a METHYLATED cytosine
        call = int(call)
        call = 1 if call == 0 else 2 # 1 unmethylated, 2 methylated
        call = call * strand

        # process if line corresponds to target chromosome and ctx
        if target_ctx == ctx and chrname == targetchr: 
            if read_id not in read2data:
                read2data[read_id] = [[],[]]

            pos  = int(pos)
            read2data[read_id][0].append(pos)
            read2data[read_id][1].append(call)

        # break loop if chunk was read completely
        read_bytes += f.tell() - new_fpos
        new_fpos    = f.tell()
        if read_bytes > chunk_sz:
            break

    return read2data

def main(args):
    if args.nproc < 0:
        args.nproc = mp.cpu_count()

    # chunk size processed by a worker when header lines are excluded
    n_bytes  = os.path.getsize(args.file)
    chunk_sz = math.ceil(n_bytes / args.nproc)

    print(f"num processes: {args.nproc}, chunk_sz: {chunk_sz}")

    pool = mp.Pool(processes=args.nproc)
                   #initializer=ccc,
                   #initargs=asdfasdf)

    workers = []
    for idx in range(args.nproc):
        p = pool.apply_async(
            parse_chunk, (
            args.target, 
             args.file, 
             idx, 
             chunk_sz, 
             args.ctx,
            )
        )
        workers.append(p)

    read2row = {}
    rows, cols, data = [], [], []
    for i, p in enumerate(workers):
        print(i)
        maps = p.get()
        
        for read_id, (_cols, _data) in maps.items():
            if read_id not in read2row:
                read2row[read_id] = len(read2row)
            rows.append(np.array([read2row[read_id]]*len(_cols)))
            cols.append(np.array(_cols))
            data.append(np.array(_data))
            # if read_id not in read2str:
            #     read2str[read_id] = _strand[0]
    pool.close()
    pool.join()

    rows   = np.hstack(rows)
    cols   = np.hstack(cols)
    data   = np.hstack(data).astype(np.int8)

    nrows = len(read2row)
    ncols = read_chrlen_file(args.chrl, args.target) # chromosome length

    # create a M-times-N sparse matrix
    # M: number of reads
    # N: number of bases in the chromosome
    print('Creating index')
    spmat = coo_array((data, (rows, cols)),
                      shape=(nrows, ncols),
                      dtype=scipy.int8)
    print('Converting format')
    spmat = spmat.tocsr() # to CSR format

    print('Saving index')
    mode = 'r+' if os.path.exists(args.out) else 'w' # update h5 file if it exists
    h5   = h5py.File(args.out, mode)

    if args.target in h5:
        del h5[args.target] # delete group if it exists

    grp = h5.create_group(args.target)
    grp.attrs['nrows'], grp.attrs['ncols'] = spmat.shape
    grp.create_dataset('indptr',  data=spmat.indptr)
    grp.create_dataset('data',    data=spmat.data)
    grp.create_dataset('indices', data=spmat.indices)

    # save read names associated to each row
    read2row = list(read2row.keys())
    dt       = h5py.special_dtype(vlen=str) 
    read2row = np.array(read2row, dtype=dt)
    grp.create_dataset('read_ids', data=read2row)

    h5.close()

def argparser():
    # 'Indexing single fast5 files.',
    parser = argparse.ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, add_help=False)
    inp = parser.add_argument_group('INPUT')
    inp.add_argument('--file', '-i',
                     action='store',
                     type=str,
                     required=True,
                     help='file containing DSP calls')
    inp.add_argument('--target', '-t',
                     action='store',
                     type=str,
                     required=True,
                     help='Target chromosome')
    inp.add_argument('--chrl', '-l',
                     action='store',
                     type=str,
                     required=True,
                     help='tsv file containing chromosome lengths')
    inp.add_argument('--ctx', '-c',
                     action='store',
                     type=str,
                     required=False,
                     default='CG',
                     choices=('CG', 'CHG', 'CHH'),
                     help='filter entries according to context')
    inp.add_argument('--nproc', '-p',
                     action='store',
                     type=int,
                     required=False,
                     default=-1,
                     help='number of processes')
    out = parser.add_argument_group('OUTPUT')
    out.add_argument('--out', '-o',
                     action='store',
                     type=str,
                     required=True,
                     help='output file used for storing the index')

    return parser
