from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
    ext_modules = cythonize([Extension("python_wrapper", ["python_wrapper.pyx", 'weak_strong_4d_c.c'], 
								libraries=[])])
)
