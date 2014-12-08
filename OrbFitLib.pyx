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
            self.spiced = False
    
    def __init__(self, fn=None):
        self.AU = 149597870.7
        self.mu = 132889724474.839203
        self.lines = []
        self.radec = []
        self.uncert = []
        if fn!=None: self.read_file(fn)
        self.maxe = 0.25
        self.mina = 25*self.AU
        self.maxa = 100*self.AU

    # Read an MPC file
    def read_file(self, fn):
        with open(fn) as f:
            while True:
                line1 = f.readline()
                line2 = f.readline()
                if not line1 or not line2: break
                self.lines.append(self.MPC_Line(line1,line2))

    # Add spice information
    def add_spice(self, kernel):
        s2r = np.pi/(180.*3600.)
        spice.furnsh(kernel)
        for l in self.lines:
            if l.spiced: continue
            l.et = spice.str2et('JD {0:.15f} UTC'.format(l.JD))
            r2,lt = spice.spkpos("10",l.et,"J2000","LT","399")
            l.hst_eq = np.array(r2)-l.hst_geo_eq
            self.radec = np.append(self.radec,[l.ra,l.dec])
            self.uncert = np.append(self.uncert,[0.001*s2r/15.,0.01*s2r])
            l.spiced = True
        spice.kclear()
        self.et0 = self.lines[0].et

    # Predict ra and dec for a given orbit
    def predict(self, orbit):
        c = 299792.458
        radec_calc = []
        for l in self.lines:
            r1,v1 = orbit.rv(l.et)
            r = r1 + l.hst_eq
            r1,v1 = orbit.rv(l.et - np.linalg.norm(r)/c)
            r = r1 + l.hst_eq
            ra = math.atan2(r[1],r[0])
            dec = math.atan(math.sin(ra)*r[2]/r[1])
            if ra<0: ra += 2.*math.pi
            radec_calc = np.append(radec_calc,[ra,dec])
        return radec_calc

    # Calculate RMS error of an orbit (correcting for cos(dec))
    def rms(self, state):
        orbit = Conic()
        orbit.setup_state(self.et0,self.mu,state)
        r2mas = 3600e3*180./np.pi
        radec_calc = self.predict(orbit)
        sumsq = 0
        for i in range(len(self.radec)//2):
            ra1,de1 = (self.radec[2*i],self.radec[2*i+1])
            ra2,de2 = (radec_calc[2*i],radec_calc[2*i+1])
            dr = ra1*math.sin(de1) - ra2*math.sin(de2)
            dd = de1 - de2
            sumsq += dr*dr + dd*dd
        return math.sqrt(sumsq/len(self.radec))*r2mas

    # Calculate log likelihood of a given orbit
    def lnlike(self, orbit):
        model = self.predict(orbit)
        inv_sigma2 = 1.0/(self.uncert**2 + model**2)
        return -0.5*(np.sum((self.radec-model)**2*inv_sigma2))

    # Apply priors
    def lnprior(self, orbit):
        if 0 < orbit.e < self.maxe and self.mina < orbit.rp/(1.-orbit.e) < self.maxa:
            return 0.0
        return -np.inf

    # Calculate log probability of a given state vector
    def lnprob(self, state):
        orb = Conic()
        orb.setup_state(self.et0,self.mu,state)
        lp = self.lnprior(orb)
        if not np.isfinite(lp): return -np.inf
        return lp + self.lnlike(orb)

####################################################

class deltav:

    def __init__(self, ssc,et,orbit):
        self.et0 = et; self.orbit = orbit
        self.rsc,self.vsc = (ssc[:3], ssc[3:])

    # Calculate a transfer velocity, given an encounter time
    def calc_vt(self, tf):
        rf,vf = self.orbit.rv(tf)
        return lambert_transfer(self.orbit.mu,self.rsc,rf,tf-self.et0)

    # Return a transfer orbit, given an encounter time
    def calc_orb(self, tf):
        vt = self.calc_vt(tf)
        ret = Conic()
        ret.setup_rv(self.orbit.t0,self.orbit.mu,self.rsc,vt)
        return ret

    # Calculate the delta v for an encounter time
    def calc_dv(self, tf):
        return np.linalg.norm(self.vsc-self.calc_vt(tf))

    # Find the encounter time with minimum delta v 
    def opt_tf(self, guess):
        return opt.fmin(self.calc_dv,guess,disp=False)[0]

