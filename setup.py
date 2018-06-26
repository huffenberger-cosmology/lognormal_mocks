from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib
import os

numpy_inc = os.path.join(get_python_lib(plat_specific=1), 'numpy/core/include')

module2 =  Extension('lognormal_mocks.mod',
                     sources = ['lognormal_mocks/src/mod.c','lognormal_mocks/src/healpix_legendre_c.c','lognormal_mocks/src/lognormal_mocks_stats.c','lognormal_mocks/src/symmetric_legendre.c'],
                     include_dirs = [numpy_inc,'lognormal_mocks/include/'],
                     libraries=['gsl','gslcblas','fftw3','healpix','healpix_cxx','cfitsio','lapack','gomp'],
                     extra_compile_args=['-fPIC','-Wall','-g'])

setup (name = 'lognormal_mocks',
       version = '0.1',
       url='https://github.com/huffenberger-cosmology/lognormal_mock',
       description = 'lognormal mocks function',
       packages =['lognormal_mocks'],
       ext_modules = [module2]
       )
