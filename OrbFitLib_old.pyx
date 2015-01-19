# distutils: language = c++
# distutils: extra_compile_args = -std=c++11 -O3

# Numpy and math
import numpy as np
cimport numpy as np
import math
import scipy.optimize as opt
from PyNBody import Conic
import SimSpice as spice

# C++ Standard Library vector class
from libcpp.vector cimport vector
from libcpp.string cimport string

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

def make_helioline(orb):
    AU = 149597870.700
    r2d = 180./math.pi
    a = orb.rp/(AU*(1-orb.e))
    return '{0:.14E} {1:.14E} {2:.14E} {3:.14E} {4:.14E} {5:.14E} {6:.3f}'.format(
            a,orb.e,orb.i*r2d,orb.O*r2d,orb.w*r2d,orb.M0*r2d,orb.t0)

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
    o2 = Conic()
    o2.setup_rv(orbit.t0,orbit.mu,ec2eq(r),ec2eq(v))
    return o2

def eq2ec_orbit(orbit):
    r,v = orbit.rv(orbit.t0)
    o2 = Conic()
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
    
    def __init__(self, fn=None):
        self.lines = []
        if fn!=None: self.read_file(fn)

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
            l.et = spice.str2et('JD {0:.15f} UTC'.format(l.JD))
            r2,lt = spice.spkpos("10",l.et,"J2000","LT","399")
            l.hst_eq = np.array(r2)-l.hst_geo_eq
        spice.kclear()

    def sumsq(self, orbit):
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

    def residuals(self, orbit):
        cdef double r2a = 3600.*180./np.pi
        cdef double c = 299792.458*(24.*3600.)
        cdef double ra, dec, dr, dd
        cdef double AU = 149597870.7
        ras = []
        des = []
        dists = []
        eras = []
        edes = []
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
            ras.append(ra)
            des.append(dec)
            dists.append(np.linalg.norm(r)/AU)
            eras.append(dr)
            edes.append(dd)
        return ras,des,dists,eras,edes

    def rms(self, orbit):
        return math.sqrt( self.sumsq(orbit) / float(2.*len(self.lines)) )

    def chi2(self, orbit,scatter):
        return self.sumsq(orbit) / float(scatter)

    def probability(self, orbit,scatter):
        return math.exp(-0.5*self.chi2(orbit,scatter))

####################################################

class deltav:

    def __init__(self, ssc,jd,orbit):
        self.jd0 = jd; self.orbit = orbit
        self.rsc,self.vsc = (ssc[:3], ssc[3:])

    def calc_vt(self, tf):
        day = 24.*3600.
        rf,vf = self.orbit.rv(tf)
        return lambert_transfer(self.orbit.mu,self.rsc,rf,tf-self.jd0)/day

    def calc_dv(self, tf):
        return np.linalg.norm(self.vsc-self.calc_vt(tf))

    def opt_tf(self, guess):
        return opt.fmin(self.calc_dv,guess,disp=False)[0]

