from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib
import os

numpy_inc = os.path.join(get_python_lib(plat_specific=1), 'numpy/core/include')

module2 =  Extension('lognormal_mocks.lognormal_mocks_mod',
                     sources = ['code/src/mod.c','code/src/healpix_legendre_c.c','code/src/lognormal_mocks_stats.c','code/src/symmetric_legendre.c'],
                     include_dirs = [numpy_inc,'code/include/'],
                     libraries=['gsl','gslcblas','healpix','cfitsio','gomp','lapack'],
                     extra_compile_args=['-fPIC','-Wall','-g'])

setup (name = 'lognormal_mocks',
       version = '0.1',
       url='https://github.com/huffenberger-cosmology/lognormal_mocks',
       description = 'lognormal mocks function',
       packages =['code'],
       ext_modules = [module2]
       )

