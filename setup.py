from distutils.core import setup
from Cython.Build import cythonize

setup(
        name = "OrbFitLib",
        ext_modules = cythonize('OrbFitLib.pyx'),
        )
