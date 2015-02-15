# distutils: language = c++
# distutils: extra_link_args = cspice/lib/cspice.a

import os

# Import numpy 
cimport numpy as np
import numpy as np

# Import C++ string and vector libraries
from libcpp.string cimport string
from libcpp.vector cimport vector

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
    void pxform_c(char *,char *,double,double[3][3])

# Import our helper functions for SPICE's archaic preprocessor-defined objects
cdef extern from "spice_helper.hpp":
    vector[int] spkobjects(string kernel)
    vector[ vector[double] ] spkcoverage(string, int)
    void write_spk5_c(string filename,int objid,int cntrid,double epoch,
            double cntrgm,vector[double] state,double rang, string cframe)
    void write_spk3_c(string filename,int body,int center,string frame,
            double first,double last,double intlen,unsigned n,unsigned polydg,
            double* cdata, double btime)

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
    cdef np.ndarray[double, ndim=1, mode="c"] elts2 = np.array(elts)
    cdef np.ndarray[double, ndim=1, mode="c"] state = np.zeros(6)
    conics_c(&elts2[0],et,&state[0])
    return state

def oscelt(state,double et,double mu):
    cdef np.ndarray[double, ndim=1, mode="c"] state2 = np.array(state)
    cdef np.ndarray[double, ndim=1, mode="c"] elts = np.zeros(8)
    oscelt_c(&state2[0],et,mu,&elts[0])
    return elts

def spkobj(kernel):
    cdef string k2 = kernel.encode('UTF-8')
    return np.array(spkobjects(k2),dtype=int)

def spkcov(kernel,int obj):
    cdef string k2 = kernel.encode('UTF-8')
    return np.array(spkcoverage(k2,obj))

def write_spk5(filename,int objid,int cntrid,double epoch,double cntrgm,state,
        rang=None,cframe=None,overwrite=False):
    if os.path.isfile(filename):
        if overwrite: os.remove(filename)
        else:
            print('Error! Not clobbering existing file '+filename)
            return
    cdef string fn2 = filename.encode('UTF-8')
    cdef np.ndarray[double, ndim=1, mode="c"] state2 = state
    if rang==None: rang=365.25*24.*3600.
    if cframe==None: cframe='J2000'
    cdef string cf2 = cframe.encode('UTF-8')
    write_spk5_c(fn2,objid,cntrid,epoch,cntrgm,state2,rang,cf2)

def write_spk3(filename,int body,int center,frame,double first,double last,
        double intlen,unsigned n,unsigned polydg,cdata,double btime,overwrite=False):
    if os.path.isfile(filename):
        if overwrite: os.remove(filename)
        else:
            print('Error! Not clobbering existing file '+filename)
            return
    cdef string fn2 = filename.encode('UTF-8')
    cdef string cf2 = frame.encode('UTF-8')
    cdef np.ndarray[double, ndim=3, mode="c"] cdata2 = np.ascontiguousarray(np.zeros((n,6,polydg+1)))
    for i in range(n):
        for j in range(6):
            for k in range(polydg+1):
                cdata2[i][j][k] = cdata[i][j][k]
    write_spk3_c(fn2,body,center,cf2,first,last,intlen,n,polydg,&cdata2[0,0,0],btime)

def pxform(fr_frame,to_frame,double et):
    cdef string f2 = fr_frame.encode('UTF-8')
    cdef string t2 = to_frame.encode('UTF-8')
    #cdef np.ndarray[double, ndim=2, mode="c"] mat = np.zeros((3,3))
    cdef double mat[3][3]
    pxform_c(f2.c_str(),t2.c_str(),et,mat)
    mat2 = np.zeros((3,3))
    for i in range(3):
        for j in range(3): mat2[i][j] = mat[i][j]
    return mat2
