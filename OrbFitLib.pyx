# distutils: language = c++
# distutils: extra_compile_args = -std=c++11 -O3

import numpy as np

# C++ Standard Library vector class
from libcpp.vector cimport vector

# Import the parts of the C++ Conic class we need
cdef extern from "conic.hpp":
    cdef cppclass Conic:
        double rp,e,i,O,w,M0,t0,mu
        Conic() except +
        void setup_elements(vector[double])
        void setup_state(double,double,vector[double])
        vector[double] elements(double)
        vector[double] state(double)

# Python Orbit class
cdef class Orbit:

    # Contain a copy of the C++ Conic class
    cdef Conic *thisptr
    def __cinit__(self):
        self.thisptr = new Conic()
    def __del__(self):
        del self.thisptr

    # Directly get the orbital elements
    property rp:
        def __get__(self): return self.thisptr.rp
    property e:
        def __get__(self): return self.thisptr.e
    property i:
        def __get__(self): return self.thisptr.i
    property O:
        def __get__(self): return self.thisptr.O
    property w:
        def __get__(self): return self.thisptr.w
    property M0:
        def __get__(self): return self.thisptr.M0
    property t0:
        def __get__(self): return self.thisptr.t0
    property mu:
        def __get__(self): return self.thisptr.mu

    # Setup the Conic class
    def setup_elements(self, elts):
        self.thisptr.setup_elements(elts)
    def setup_state(self, t0,mu,sv):
        self.thisptr.setup_state(t0,mu,sv)
    def setup_rv(self, t0,mu,r,v):
        sv = list(r)+list(v)
        self.thisptr.setup_state(t0,mu,sv)

    def setup_helioline(self, line):
        AU = 149597870.700
        day = 24.*3600.
        mu = 132712440023.310 *(day*day)
        d2r = np.pi/180.
        a,e,i,O,w,M0,t0 = [float(i) for i in line.split()]
        self.setup_elements([a*(1-e)*AU,e,i*d2r,O*d2r,w*d2r,M0*d2r,t0,mu])

    # Extract elements or state from the class
    def elements(self, t):
        return np.array(self.thisptr.elements(t))
    def state(self, t):
        return np.array(self.thisptr.state(t))
    def rv(self, t):
        sv = self.state(t)
        return sv[:3],sv[3:]

####################################################

# Import the parts of the C++ Lambert class we need
cdef extern from "lambert.hpp":
    cdef cppclass Lambert:
        Lambert() except +
        vector[double] transfer(double,vector[double],vector[double],double)

# Python Lambert function
"""cdef class lambert_transfer:

    cdef Lambert *thisptr
    def __cinit__(self):
        self.thisptr = new Lambert()
    def __del__(self):
        del self.thisptr
    
    def velocity(self, mu,r1,r2,dt):
        return np.array( self.thisptr.transfer(mu,r1,r2,dt) )"""

def lambert_transfer(mu,r1,r2,dt):
    cdef Lambert *lamb = new Lambert()
    res = lamb.transfer(mu,r1,r2,dt)
    del lamb
    return np.array(res)
