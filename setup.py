from distutils.core import setup
from Cython.Build import cythonize
from os.path import isfile
from sys import exit

# Yell at the user to RTFM
def nospice():
    print('\n"He who controls the spice controls the universe."\n')
    print("You don't have SPICE downloaded and/or extracted!")
    print("Follow the directions in README.md and try again.\n")
    exit()

# Check that the SPICE header is where I need it
if not isfile('cspice/include/SpiceUsr.h'): nospice()

# Check that the SPICE library is where I need it
if not isfile('cspice/lib/cspice.a'): nospice()

# The Cython build call
setup(
        name = "OrbFitLib",
        ext_modules = cythonize(('OrbFitLib.pyx','SimSpice.pyx')
            # Don't need these with Cython 0.17+
            ,language='c++',extra_compile_args='-O3'
            ),
        )
