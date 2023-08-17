# dspcall2spmat

Small repository containing computational utilies to foster reproducibility.

## Install

Download repository, and go inside the downloaded repo and perform

```
pip install .
```

## File convertion

Convert the calls into a h5 file

```shell
dspcall2spmat pred2spmat --file dsp_call.tsv \
                         --nproc 100 \
			 --out mychr1.h5 \
			 --target Chr1 \
			 --chrl Col-CEN_v1.2.lengths \
			 --ctx CG
```

Load the created h5 file

```python
import h5py
import scipy

chrn = 'Chr1'

h5 = h5py.File('mychr1.h5')
read_ids = h5[chrn]['read_ids'][:]
nrows, ncols = h5[chrn].attrs['nrows'], h5[chrn].attrs['ncols']
spmat = scipy.sparse.csr_array((h5[chrn]['data'][:], h5[chrn]['indices'][:], h5[chrn]['indptr'][:]),
                              shape=(nrows, ncols),
                              dtype=scipy.int8)
spmat.shape
h5.close()
```

Now the loaded h5 file can be used to construct a sparse matrix from which the methylation levels can be derived in a strand-specific way

```python
profs = []
for strand in [1, -1]: # [forward, reverse]
    uC = (spmat == strand*1).sum(axis=0).flatten() # num of unmethylated CGs
    mC = (spmat == strand*2).sum(axis=0).flatten() # num of methylated CGs
    prof = mC / (uC + mC + 1e-7) # calculate methylated levels
    mask = np.logical_and(uC == 0, mC == 0) # assign nan values to uninformative positions
    prof[mask] = np.nan
    profs.append(prof)
```