from distutils.core import setup
from setuptools import find_packages

package_name = 'dspcall2spmat'

setup(name=package_name,
      description="Converting DSP calls into a sparse matrix",
      version='1.0',
      author="Alejandro A. Edera",
      author_email="edera.alejandro@gmail.com",
      packages=find_packages(),
      entry_points = {
          'console_scripts': [
              '{0} = {0}:main'.format(package_name)
          ]
      },
      setup_requires=['numpy', 'h5py'],
      install_requires=['numpy', 'h5py']
)
