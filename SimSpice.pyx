# distutils: language = c++
# distutils: extra_link_args = cspice/lib/cspice.a

# Import numpy 
cimport numpy as np
import numpy as np

# Import C++ string library to minimize char* sillyness
from libcpp.string cimport string

# Import NAIF SPICE functions
cdef extern from "cspice/include/SpiceUsr.h":
    void furnsh_c(char *)
    void unload_c(char *)
    void kclear_c()
    void spkgps_c(int,double,char *,int,double *,double *)
    void spkgeo_c(int,double,char *,int,double *,double *)
    void spkpos_c(char *,double,char *,char *,char *,double *,double *)
    void spkezr_c(char *,double,char *,char *,char *,double *,double *)
    void str2et_c(char *,double *)
    void timout_c(double,char *,int,char *)
    void tpictr_c(char *,int,int,char *,int *,char *)
    void recrad_c(double *,double *,double *,double *)
    void radrec_c(double,double,double,double *)
    void conics_c(double *,double,double *)
    void oscelt_c(double *,double,double,double *)

def furnsh(kernel):
    cdef string k2 = kernel.encode('UTF-8')
    furnsh_c(k2.c_str())

def unload(kernel):
    cdef string k2 = kernel.encode('UTF-8')
    unload_c(k2.c_str())

def kclear():
    kclear_c()

def spkgps(int targ,double et,ref,int obs):
    cdef string ref2 = ref.encode('UTF-8')
    cdef double lt        
    cdef np.ndarray[double, ndim=1, mode="c"] state = np.zeros(3)
    spkgps_c(targ,et,ref2.c_str(),obs,&state[0],&lt)
    return state,lt

def spkgeo(int targ,double et,ref,int obs):
    cdef string ref2 = ref.encode('UTF-8')
    cdef double lt
    cdef np.ndarray[double, ndim=1, mode="c"] state = np.zeros(6)
    spkgeo_c(targ,et,ref2.c_str(),obs,&state[0],&lt)
    return state,lt

def spkpos(targ,double et,ref,abcorr,obs):
    cdef string targ2 = targ.encode('utf-8')
    cdef string ref2 = ref.encode('utf-8')
    cdef string abcorr2 = abcorr.encode('utf-8')
    cdef string obs2 = obs.encode('utf-8')
    cdef double lt
    cdef np.ndarray[double, ndim=1, mode="c"] state = np.zeros(3)
    spkpos_c(targ2.c_str(),et,ref2.c_str(),abcorr2.c_str(),obs2.c_str(),&state[0],&lt)
    return state,lt

def spkezr(targ,double et,ref,abcorr,obs):
    cdef string targ2 = targ.encode('utf-8')
    cdef string ref2 = ref.encode('utf-8')
    cdef string abcorr2 = abcorr.encode('utf-8')
    cdef string obs2 = obs.encode('utf-8')
    cdef double lt
    cdef np.ndarray[double, ndim=1, mode="c"] state = np.zeros(6)
    spkezr_c(targ2.c_str(),et,ref2.c_str(),abcorr2.c_str(),obs2.c_str(),&state[0],&lt)
    return state,lt

def str2et(s):
    cdef string s2 = s.encode('UTF-8')
    cdef double et
    str2et_c(s2.c_str(),&et)
    return et

def timout(double et,picture):
    cdef string picture2 = picture.encode('UTF-8')
    cdef char buff[256]
    timout_c(et,picture2.c_str(),256,buff)
    return buff.decode('UTF-8')

def tpictr(example):
    cdef string example2 = example.encode('UTF-8')
    cdef char buff1[256]
    cdef char buff2[256]
    cdef int ok = 0
    tpictr_c(example2.c_str(),256,256,buff1,&ok,buff2)
    return buff1.decode('UTF-8'),ok,buff2.decode('UTF-8')

def recrad(rectan):
    cdef double ran,ra,dec
    cdef np.ndarray[double, ndim=1, mode="c"] state = rectan
    recrad_c(&state[0],&ran,&ra,&dec)
    return ran,ra,dec

def radrec(double ran,double ra,double dec):
    cdef np.ndarray[double, ndim=1, mode="c"] state = np.zeros(3)
    radrec_c(ran,ra,dec,&state[0])
    return state

def conics(elts,double et):
    cdef np.ndarray[double, ndim=1, mode="c"] elts2 = elts
    cdef np.ndarray[double, ndim=1, mode="c"] state = np.zeros(6)
    conics_c(&elts2[0],et,&state[0])
    return state

def oscelt(state,double et,double mu):
    cdef np.ndarray[double, ndim=1, mode="c"] state2 = state
    cdef np.ndarray[double, ndim=1, mode="c"] elts = np.zeros(8)
    oscelt_c(&state2[0],et,mu,&elts[0])
    return elts
