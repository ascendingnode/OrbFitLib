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
    void reclat_c(double *,double *,double *,double *)
    void conics_c(double *,double,double *)
    void oscelt_c(double *,double,double,double *)
    void pxform_c(char *,char *,double,double[3][3])
    void sxform_c(char *,char *,double,double[6][6])
    void getfov_c(int,int,int,int,char *,char *,double[3],int *,double[][3])
    void bodn2c_c(char *,int *,int *)
    void subpnt_c(char *,char *,double,char *,char *,char *,double[3],double *,double[3])
    void subslr_c(char *,char *,double,char *,char *,char *,double[3],double *,double[3])
    void recgeo_c(double[3],double,double,double *,double *,double *)
    void georec_c(double,double,double, double,double, double[3])

# Import our helper functions for SPICE's archaic preprocessor-defined objects
cdef extern from "spice_helper.hpp":
    vector[int] spkobjects(string kernel)
    vector[ vector[double] ] spkcoverage(string, int)
    void write_spk5_c(string filename,int objid,int cntrid,string cframe,double etbeg,
            double etend,double cntrgm,int nstate,double *cstate,double *cepoch)
    void write_spk3_c(string filename,int body,int center,string frame,
            double first,double last,double intlen,unsigned n,unsigned polydg,
            double* cdata, double btime)
    void write_spk2_c(string filename,int body,int center,string frame,
            double first,double last,double intlen,unsigned n,unsigned polydg,
            double* cdata, double btime, string segid, char append)
    void write_spk9_c(string filename,int body,int center,string frame,
            double first,double last,int n,int degree,double *states,double *epochs)

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

def reclat(rectan):
    cdef double ran,ra,dec
    cdef np.ndarray[double, ndim=1, mode="c"] state = rectan
    reclat_c(&state[0],&ran,&ra,&dec)
    return ran,ra,dec

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

def write_spk5(filename,int objid,int cntrid,epochs,double cntrgm,states,cframe=None,overwrite=False):
    if os.path.isfile(filename):
        if overwrite: os.remove(filename)
        else:
            print('Error! Not clobbering existing file '+filename)
            return
    cdef string fn2 = filename.encode('UTF-8')
    n = len(epochs)
    cdef np.ndarray[double, ndim=1, mode="c"] cepoch = epochs
    cdef np.ndarray[double, ndim=2, mode="c"] cstate = np.ascontiguousarray(np.zeros((n,6)))
    for i in range(n):
        for j in range(6):
            cstate[i][j] = states[i][j]
    if cframe==None: cframe='J2000'
    cdef string cf2 = cframe.encode('UTF-8')
    write_spk5_c(fn2,objid,cntrid,cf2,epochs[0],epochs[-1],cntrgm,n,&cstate[0,0],&cepoch[0])

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

def write_spk2(filename,int body,int center,frame,double first,double last,
        double intlen,unsigned n,unsigned polydg,cdata,double btime,overwrite=False,segid=None,append=False):
    if segid is None: segid = "1"
    cdef string segid2 = segid.encode('UTF-8')
    if os.path.isfile(filename):
        if overwrite: 
            if not append:
                os.remove(filename)
        else:
            print('Error! Not clobbering existing file '+filename)
            return
    cdef string fn2 = filename.encode('UTF-8')
    cdef string cf2 = frame.encode('UTF-8')
    cdef np.ndarray[double, ndim=3, mode="c"] cdata2 = np.ascontiguousarray(np.zeros((n,3,polydg+1)))
    for i in range(n):
        for j in range(3):
            for k in range(polydg+1):
                cdata2[i][j][k] = cdata[i][j][k]
    write_spk2_c(fn2,body,center,cf2,first,last,intlen,n,polydg,&cdata2[0,0,0],btime,segid2,append)

def write_spk9(filename,int body,int center,epochs,states,degree=15,cframe=None,overwrite=False):
    if os.path.isfile(filename):
        if overwrite: os.remove(filename)
        else:
            print('Error! Not clobbering existing file '+filename)
            return
    cdef string fn2 = filename.encode('UTF-8')
    n = len(epochs)
    cdef np.ndarray[double, ndim=1, mode="c"] cepoch = epochs
    cdef np.ndarray[double, ndim=2, mode="c"] cstate = np.ascontiguousarray(np.zeros((n,6)))
    for i in range(n):
        for j in range(6):
            cstate[i][j] = states[i][j]
    if cframe==None: cframe='J2000'
    cdef string cf2 = cframe.encode('UTF-8')
    write_spk9_c(fn2,body,center,cf2,epochs[0],epochs[-1],n,degree,&cstate[0,0],&cepoch[0])

def pxform(fr_frame,to_frame,double et):
    cdef string f2 = fr_frame.encode('UTF-8')
    cdef string t2 = to_frame.encode('UTF-8')
    cdef double mat[3][3]
    pxform_c(f2.c_str(),t2.c_str(),et,mat)
    mat2 = np.zeros((3,3))
    for i in range(3):
        for j in range(3): mat2[i][j] = mat[i][j]
    return mat2

def sxform(fr_frame,to_frame,double et):
    cdef string f2 = fr_frame.encode('UTF-8')
    cdef string t2 = to_frame.encode('UTF-8')
    cdef double mat[6][6]
    sxform_c(f2.c_str(),t2.c_str(),et,mat)
    mat2 = np.zeros((6,6))
    for i in range(6):
        for j in range(6): mat2[i][j] = mat[i][j]
    return mat2

def getfov(int instid):
    cdef int room = 256
    cdef int shapelen = 256
    cdef int framelen = 256
    cdef char shape[256]
    cdef char frame[256]
    cdef np.ndarray[double, ndim=1, mode="c"] bsight = np.zeros(3)
    cdef int n
    cdef double bounds[256][3]
    getfov_c(instid,room,shapelen,framelen,shape,frame,&bsight[0],&n,bounds)
    bounds2 = np.zeros((n,3))
    for i in range(n):
        for j in range(3):
            bounds2[i][j] = bounds[i][j]
    return shape.decode('UTF-8'),frame.decode('UTF-8'),bsight,bounds2

def bodn2c(name):
    cdef int objid=0, found
    cdef string n2 = name.encode('UTF-8')
    bodn2c_c(n2.c_str(),&objid,&found)
    if found>0: f2 = True
    else: f2 = False
    return objid,f2

def subpnt(method,target,double et,fixref,abcorr,obsrvr):
    cdef string method2 = method.encode('UTF-8')
    cdef string target2 = target.encode('UTF-8')
    cdef string fixref2 = fixref.encode('UTF-8')
    cdef string abcorr2 = abcorr.encode('UTF-8')
    cdef string obsrvr2 = obsrvr.encode('UTF-8')
    cdef np.ndarray[double, ndim=1, mode="c"] spoint = np.zeros(3)
    cdef double trgepc
    cdef np.ndarray[double, ndim=1, mode="c"] srfvec = np.zeros(3)
    subpnt_c(method2.c_str(),target2.c_str(),et,fixref2.c_str(),abcorr2.c_str(),
            obsrvr2.c_str(),&spoint[0],&trgepc,&srfvec[0])
    return spoint,trgepc,srfvec

def subslr(method,target,double et,fixref,abcorr,obsrvr):
    cdef string method2 = method.encode('UTF-8')
    cdef string target2 = target.encode('UTF-8')
    cdef string fixref2 = fixref.encode('UTF-8')
    cdef string abcorr2 = abcorr.encode('UTF-8')
    cdef string obsrvr2 = obsrvr.encode('UTF-8')
    cdef np.ndarray[double, ndim=1, mode="c"] spoint = np.zeros(3)
    cdef double trgepc
    cdef np.ndarray[double, ndim=1, mode="c"] srfvec = np.zeros(3)
    subslr_c(method2.c_str(),target2.c_str(),et,fixref2.c_str(),abcorr2.c_str(),
            obsrvr2.c_str(),&spoint[0],&trgepc,&srfvec[0])
    return spoint,trgepc,srfvec

def recgeo(rectan,double re,double f):
    cdef double lon, lat, alt
    cdef np.ndarray[double, ndim=1, mode="c"] state = np.array(rectan,dtype=float)
    recgeo_c(&state[0],re,f,&lon,&lat,&alt)
    return lon, lat, alt

def georec(double lon,double lat,double alt,double re,double f):
    cdef np.ndarray[double, ndim=1, mode="c"] rectan = np.zeros(3)
    georec_c(lon,lat,alt,re,f,&rectan[0])
    return rectan
