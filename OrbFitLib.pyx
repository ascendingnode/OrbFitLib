# distutils: language = c++
# distutils: extra_compile_args = -std=c++11 -O3

# Numpy and math
import numpy as np
#cimport numpy as np
import math

# Spice wrapper library
import spice

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

    # Read a heliocentric orbit file
    def setup_helioline(self, line):
        cdef double AU = 149597870.700
        cdef double day = 24.*3600.
        cdef double mu = 132712440023.310 *(day*day)
        cdef double d2r = np.pi/180.
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
def lambert_transfer(mu,r1,r2,dt):
    cdef Lambert *lamb = new Lambert()
    res = lamb.transfer(mu,r1,r2,dt)
    del lamb
    return np.array(res)

####################################################

def ec2eq(ec):
    cdef double eta = (23.+(26/60.)+(21.406/3600.))*math.pi/180.
    cdef double se = math.sin(eta)
    cdef double ce = math.cos(eta)
    eq = np.array(ec)
    eq[1] =  ce*ec[1] - se*ec[2]
    eq[2] =  se*ec[1] + ce*ec[2]
    return eq

def eq2ec(eq):
    cdef double eta = (23.+(26/60.)+(21.406/3600.))*math.pi/180.
    cdef double se = math.sin(eta)
    cdef double ce = math.cos(eta)
    ec = np.array(eq)
    ec[1] =  ce*eq[1] + se*eq[2]
    ec[2] = -se*eq[1] + ce*eq[2]
    return ec

def ec2eq_orbit(orbit):
    r,v = orbit.rv(orbit.t0)
    o2 = Orbit()
    o2.setup_rv(orbit.t0,orbit.mu,ec2eq(r),ec2eq(v))
    return o2

def eq2ec_orbit(orbit):
    r,v = orbit.rv(orbit.t0)
    o2 = Orbit()
    o2.setup_rv(orbit.t0,orbit.mu,eq2ec(r),eq2ec(v))
    return o2

####################################################

def MPCdate2JD(ds):
    cdef int year,month,y,m,B
    cdef double dday
    ss = ds.split()
    year,month,dday = (int(ss[0]),int(ss[1]),float(ss[2]))
    if month>2: 
        y=year; m=month
    else:
        y=year-1; m=month+12
    B = int(y//400) - int(y//100)
    return 1720996.5 + int(365.25*y) + int(30.6001*(m+1)) + B + dday

def sex2rad(sex):
    cdef int d,m,ad
    cdef double s,dd
    ss = sex.split()
    d,m,s = (int(ss[0]),int(ss[1]),float(ss[2]))
    ad = abs(d)
    dd = ad + (m/60.) + (s/3600.)
    if d<0: dd *= -1.
    return dd*np.pi/180.;

def JD2etUTC(JD):
    return spice.str2et('JD {0:.15f} UTC'.format(JD))

class MPC_File:

    class MPC_Line:
        def __init__(self, l1=None,l2=None):
            if l1==None or l2==None: return
            name,date = (l1[ 6-1:12],l1[16-1:32])
            ras ,decs = (l1[33-1:44],l1[45-1:56])
            self.JD = MPCdate2JD(date)
            self.ra = sex2rad(ras)*15.
            self.dec = sex2rad(decs)
            x = float(l2[35-1:45])
            y = float(l2[47-1:57])
            z = float(l2[59-1:69])
            self.hst_geo_eq = np.array([x,y,z])

    lines = []

    def read_file(self, fn):
        with open(fn) as f:
            while True:
                line1 = f.readline()
                line2 = f.readline()
                if not line1 or not line2: break
                self.lines.append(self.MPC_Line(line1,line2))

    def add_spice(self, kernel,dist=None):
        spice.furnsh(kernel)
        for l in self.lines:
            l.et = JD2etUTC(l.JD)
            r2,lt = spice.spkpos("10",l.et,"J2000","LT","399")
            l.hst_eq = np.array(r2)-l.hst_geo_eq
            if dist==None:
                l.direction = np.array(spice.radrec(1,l.ra,l.dec))
            else:
                AU = 149597870.7
                d = np.array(spice.radrec(dist*AU,l.ra,l.dec)) + l.hst_eq;
                l.direction = d/np.linalg.norm(d);
        spice.unload(kernel)

    def chi2(self, orbit):
        cdef double r2a = 3600.*180./np.pi
        cdef double c = 299792.458*(24.*3600.)
        cdef double err = 0, ra, dec, dr, dd
        for l in self.lines:
            r1,v1 = orbit.rv(l.JD)
            r = r1 + l.hst_eq
            r1,v1 = orbit.rv(l.JD - np.linalg.norm(r)/c)
            r = r1 + l.hst_eq
            ra = math.atan2(r[1],r[0])
            dec = math.atan(math.sin(ra)*r[2]/r[1])
            if ra<0: ra += 2.*math.pi
            dr = (l.ra-ra)*math.sin(l.dec)*r2a
            dd = (l.dec-dec)*r2a
            err += dr*dr + dd*dd
        return err

    def rms(self, orbit):
        return math.sqrt( self.chi2(orbit) / float(len(self.lines)) )
