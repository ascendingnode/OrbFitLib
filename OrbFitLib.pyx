# distutils: language = c++
# distutils: extra_compile_args = -O3

####################################################

# Standard math
import math

# Numpy and Scipy
import numpy as np
from scipy.optimize import fmin

# C++ vector for Conic
from libcpp.vector cimport vector

# Emcee for cloudiness
import emcee

# Temporary SPICE wrapper library
import SimSpice as spice

####################################################

# Import the parts of the C++ Conic class we need
cdef extern from "conic.hpp":
    cdef cppclass _Conic "Conic":
        double rp,e,i,O,w,M0,t0,mu
        _Conic() except +
        void setup_elements(vector[double])
        void setup_state(double,double,vector[double])
        vector[double] elements(double)
        vector[double] state(double)

# Python Conic class
cdef class Conic:

    # Contain a copy of the C++ Conic class
    cdef _Conic *thisptr
    def __cinit__(self):
        self.thisptr = new _Conic()
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

# Make formated string of elements
def make_helioline(orb):
    AU = 149597870.700
    r2d = 180./math.pi
    a = orb.rp/(AU*(1-orb.e))
    return '{0:.14E} {1:.14E} {2:.14E} {3:.14E} {4:.14E} {5:.14E} {6:.3f}'.format(
            a,orb.e,orb.i*r2d,orb.O*r2d,orb.w*r2d,orb.M0*r2d,orb.t0)

# Convert an ecliptic vector to equitorial
def ec2eq(ec):
    cdef double eta = (23.+(26/60.)+(21.406/3600.))*math.pi/180.
    cdef double se = math.sin(eta)
    cdef double ce = math.cos(eta)
    eq = np.array(ec)
    eq[1] =  ce*ec[1] - se*ec[2]
    eq[2] =  se*ec[1] + ce*ec[2]
    return eq

# Convert an equitorial vector to ecliptic
def eq2ec(eq):
    cdef double eta = (23.+(26/60.)+(21.406/3600.))*math.pi/180.
    cdef double se = math.sin(eta)
    cdef double ce = math.cos(eta)
    ec = np.array(eq)
    ec[1] =  ce*eq[1] + se*eq[2]
    ec[2] = -se*eq[1] + ce*eq[2]
    return ec

# Convert an ecliptic orbit to equitorial
def ec2eq_orbit(orbit):
    r,v = orbit.rv(orbit.t0)
    o2 = Conic()
    o2.setup_rv(orbit.t0,orbit.mu,ec2eq(r),ec2eq(v))
    return o2

# Convert an equitorial orbit to ecliptic
def eq2ec_orbit(orbit):
    r,v = orbit.rv(orbit.t0)
    o2 = Conic()
    o2.setup_rv(orbit.t0,orbit.mu,eq2ec(r),eq2ec(v))
    return o2

####################################################

# Convert an MPC file's date string into a decimal JD UTC
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

# Convert sexigesimal ddd:mm:ss.ssss string to decimal dd.ddddd float
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

    # Object to process and store an astrometric point from an MPC file
    # *** Currently only works with two-line spacecraft format ***
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
            self.obs_geo_eq = np.array([x,y,z])
            self.spiced = False
    
    # Constructor
    def __init__(self, fn=None):
        self.AU = 149597870.7
        # Includes the masses of the Sun, Jupiter, and Saturn
        self.mu = 132877093391.27286 
        self.lines = []
        self.radec = []
        self.inv_sigma2 = []
        if fn!=None: self.read_file(fn)
        # Rough priors for TNOs
        self.maxe = 0.9
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
        self.kernel = kernel
        s2r = np.pi/(180.*3600.)
        ra_round = 0.001*s2r*15.
        de_round = 0.01*s2r
        spice.furnsh(kernel)
        for l in self.lines:
            if l.spiced: continue
            l.et = spice.str2et('JD {0:.15f} UTC'.format(l.JD))
            # Cache lighttime corrected position of the geocenter relative to solar system barycenter
            r2 = spice.spkpos("0",l.et,"J2000","LT","399")[0]
            # Offset position of observer from solar system barycenter
            l.obs_eq = r2-l.obs_geo_eq
            # The radec array is [RA0 cos(Dec0), Dec0, RA1 cos(Dec1), Dec1, ...] in radians
            self.radec = np.append(self.radec,[l.ra*math.cos(l.dec),l.dec])
            # Default uncertainty = rounding error in MPC format
            self.inv_sigma2 = np.append(self.inv_sigma2,[ra_round**-2,de_round**-2])
            l.spiced = True
        spice.kclear()
        self.et0 = self.lines[0].et

    # Convert a state vector to an orbit
    def to_orbit(self, state):
        orb = Conic()
        orb.setup_state(self.et0,self.mu,state)
        return orb

    # Predict arbitrary geocentric RA & Dec for a given orbit
    def predict_dates(self, orbit,dates):
        if not hasattr(dates, "__len__"): dates = [dates]
        #if isinstance(dates, str): dates = [dates]
        c = 299792.458
        spice.furnsh(self.kernel)
        ret = []
        for date in dates:
            et = spice.str2et(date)
            earth = spice.spkpos("0",et,"J2000","LT","399")[0]
            r1 = orbit.rv(et)[0]
            r2 = r1 + earth
            r3 = orbit.rv(et - np.linalg.norm(r2)/c)[0]
            r = r3 + earth
            ra = math.atan2(r[1],r[0])
            dec = math.atan(math.sin(ra)*r[2]/r[1])
            if ra<0: ra += 2.*math.pi
            ret.append([ra,dec])
        spice.kclear()
        return ret

    # Predict observational RA & Dec for a given orbit
    def predict(self, orbit):
        # Seems small when you put in km/s
        c = 299792.458
        radec_calc = []
        for l in self.lines:
            # Make initial calculation of position
            r1 = orbit.rv(l.et)[0]
            r2 = r1 + l.obs_eq
            # Calculate lighttime-corrected position
            r3 = orbit.rv(l.et - np.linalg.norm(r2)/c)[0]
            r = r3 + l.obs_eq
            # Convert to RA & Dec (in radians)
            ra = math.atan2(r[1],r[0])
            dec = math.atan(math.sin(ra)*r[2]/r[1])
            if ra<0: ra += 2.*math.pi
            radec_calc = np.append(radec_calc,[ra*math.cos(dec),dec])
        return radec_calc

    # Calculate RMS error of an orbit
    def rms(self, state):
        model = self.predict(self.to_orbit(state))
        r2mas = 3600e3*180./np.pi
        #return np.sum((self.radec-model)**2/len(self.radec))*r2mas
        return np.sqrt(np.mean(np.square(self.radec-model)))*r2mas

    # Calculate Chi-Squared of a given orbit
    def X2(self, orbit):
        #return np.sum((self.radec-self.predict(orbit))**2*self.inv_sigma2)
        return np.sum(np.square(self.radec-self.predict(orbit))*self.inv_sigma2)

    # Calculate log likelihood of a given orbit
    def lnlike(self, orbit):
        # P ~ exp( -X2 / 2 )
        return -0.5*self.X2(orbit)

    # Apply priors
    def lnprior(self, orbit):
        if 0 < orbit.e < self.maxe and self.mina < orbit.rp/(1.-orbit.e) < self.maxa:
            return 0.0
        return -np.inf

    def save_kernel(self, g,fn,objid,overwrite=False):
        year = 365.25*24.*3600.
        spice.write_spk5(fn,objid,0,self.et0,self.mu,g,10*year,'J2000',overwrite)

    # Calculate log probability of a given state vector
    def lnprob(self, state):
        orb = self.to_orbit(state)
        lp = self.lnprior(orb)
        if not np.isfinite(lp): return -np.inf
        return lp + self.lnlike(orb)

    # Make an initial guess within the bounds
    def guess(self):
        a = np.random.uniform(self.mina,self.maxa)
        e = np.random.uniform(0,self.maxe)
        i = np.random.uniform(0,0.5*np.pi)
        O = np.random.uniform(0,2.*np.pi)
        w = np.random.uniform(0,2.*np.pi)
        M = np.random.uniform(0,2.*np.pi)
        orb = Conic()
        orb.setup_elements([a*(1-e),e,i,O,w,M,self.et0,self.mu])
        return orb.state(self.et0)

    # Generate a guess with less than the given RMS error
    def good_guess(self, maxrms):
        def weight(f): return -self.lnprob(f)
        while True:
            g = self.guess()
            for i in range(10):
                g = fmin(weight,g,disp=False)
            if self.rms(g)<maxrms: return g

    # Create a group of walkers for emcee
    def make_walkers(self, g,maxrms,nwalkers):
        pert = 1e-7
        rm,vm = ( np.linalg.norm(g[:3]), np.linalg.norm(g[3:]) )
        walkers = []
        for nw in range(nwalkers):
            while True:
                p = np.zeros(6)
                for i in range(3):
                    p[i  ] = np.random.normal(g[i  ],rm*pert)
                    p[i+3] = np.random.normal(g[i+3],vm*pert)
                if self.rms(p)<maxrms:
                    walkers.append(p)
                    break
        return np.array(walkers)

    # Use emcee to make an unbiased cloud of solutions
    def make_cloud(self, maxrms,nwalkers,niter1,niter2,verbose=False,threads=1):

        if verbose and threads>1: print("# Using {:} threads".format(threads))
        
        if verbose: print("# Generating initial good guess")
        g = self.good_guess(maxrms)

        if verbose: print("# Generating {:} walkers".format(nwalkers))
        walkers = self.make_walkers(g,maxrms,nwalkers)

        if verbose: print("# Running emcee for {:} iterations to burn in".format(niter1))
        sampler = emcee.EnsembleSampler(nwalkers, 6, self.lnprob, threads=threads)
        pos, prob, state = sampler.run_mcmc(walkers, niter1)

        if verbose: print("# Running emcee for {:} iterations".format(niter2))
        sampler.run_mcmc(pos, niter2)
        
        return sampler

####################################################

# Class to make delta v calculations
class deltav:

    # Constructor
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
        ret.setup_rv(self.et0,self.orbit.mu,self.rsc,vt)
        return ret

    # Calculate the delta v for an encounter time
    def calc_dv(self, tf):
        return np.linalg.norm(self.vsc-self.calc_vt(tf))

    # Find the encounter time with minimum delta v 
    def opt_tf(self, guess):
        return fmin(self.calc_dv,guess,disp=False)[0]

