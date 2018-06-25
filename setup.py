from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib
import os

numpy_inc = os.path.join(get_python_lib(plat_specific=1), 'numpy/core/include')

module1 = Extension('ylm_tools',
                    sources = ['code/ylm_tools_mod.c'],
                    include_dirs = [numpy_inc,'include'],
                    libraries=['corrfun','wcs','healpix','cfitsio','gomp'],
                    library_dirs = ["lib"],
                    extra_compile_args=['-fPIC','-Wall','-g'],
                    depends=['lib/libcorrfun.a']
                    )

setup (name = 'ylm_tools',
       version = 'VERSION',
       description = 'Help for ylm functions',
       ext_modules = [module1],
       )



module2 =  Extension('cmbtools',
                    sources = ['code/cmbtools_mod.c'],
                     include_dirs = [numpy_inc,'include','/usr/include/cfitsio/'],
                     libraries=['corrfun','gsl','gslcblas','fftw3','healpix','cfitsio','lapack','gomp'],
                    library_dirs = ["lib"],
                    extra_compile_args=['-fPIC','-Wall','-g'],
                    depends=['lib/libcorrfun.a']
                    )
setup (name = 'cmbtools',
       version = 'VERSION',
       description = 'cmbtools functions',
       ext_modules = [module2],
       )

module3 = Extension('UPP',
                    sources = ['code/UPP_mod.c'],
                    include_dirs = [numpy_inc, 'include'],
                    libraries = ['gsl','gslcblas','m','corrfun'],
                    library_dirs = ["lib"],
                    extra_compile_args = ['-fPIC', '-Wall', '-g']
                    )
setup (name = 'UPP',
       version = '01/14/2016',
       description = 'scripts for universal pressure profile: P(r), P(l) for integration along line of sight, y_UPP from P(l)',
       ext_modules = [module3],
       )
