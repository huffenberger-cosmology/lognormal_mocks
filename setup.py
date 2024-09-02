from setuptools import setup, Extension
import numpy
import os

numpy_inc = numpy.get_include()

module2 =  Extension('lognormal_mocks.lognormal_mocks_mod',
                     sources = ['lognormal_mocks/src/mod.c','lognormal_mocks/src/healpix_legendre_c.c','lognormal_mocks/src/lognormal_mocks_stats.c','lognormal_mocks/src/symmetric_legendre.c'],
                     include_dirs = [numpy_inc,'lognormal_mocks/include'],
                     libraries=['gsl','gslcblas','healpix','cfitsio','gomp','lapack','sharp'],
                     extra_compile_args=['-fPIC','-Wall','-g','-fopenmp']) #compile arg for openmp

setup (name = 'lognormal_mocks',
       version = '0.11',
       url='https://github.com/huffenberger-cosmology/lognormal_mocks',
       description = 'lognormal mocks function',
       packages =['lognormal_mocks'],
       ext_modules = [module2],
       headers = ['lognormal_mocks/include/healpix_legendre_c.h',
                  'lognormal_mocks/include/lognormal_mocks.h',
                  'lognormal_mocks/include/symmetric_legendre.h']
       )

