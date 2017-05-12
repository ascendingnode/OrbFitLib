#include "cspice/include/SpiceUsr.h"
#include <vector>
#include <string>

// Wrapper functions for SPICE's strange anachronistic SpiceCell objects

// Extract the objects in an SPK kernel
std::vector<int> spkobjects(std::string kernel) {
    SPICEINT_CELL(ids,100);
    scard_c(0,&ids);
    std::vector<int> ret;
    spkobj_c(kernel.c_str(),&ids);
    for(int i=0;i<card_c(&ids);i++) 
        ret.push_back(SPICE_CELL_ELEM_I(&ids,i));
    return ret;
}

// Extract the temporal coverages for a given object in an SPK kernel
std::vector< std::vector<double> > spkcoverage(std::string kernel, int obj) {
    SPICEDOUBLE_CELL(cover,2000);
    scard_c(0,&cover);
    spkcov_c(kernel.c_str(),obj,&cover);
    int niv = wncard_c(&cover);
    std::vector< std::vector<double> > ret(niv,std::vector<double>(2,0.0));
    for(int j=0; j<niv; j++) {
        double b,e; wnfetd_c(&cover,j,&b,&e);
        ret[j][0] = b; ret[j][1] = e;
    }
    return ret;
}

// Make a simple Type 05 (two-body) SPK file
void write_spk5_c(const std::string &filename,int objid,int cntrid,const std::string cframe,
        double etbeg,double etend,double cntrgm,int nstate,double *cstate,double *cepoch  
        /*int objid,int cntrid,double epoch,
        double cntrgm,const std::vector<double> &state,double rang,
        std::string cframe,int nstate*/) {
    std::string ifname="simple spk file", segmid="1";
    int ncomch=5000, handle;
    //double cstate[1][6]; for(unsigned i=0;i<6;i++) cstate[0][i] = state[i];
    //double cepoch[1] = {epoch}, etbeg=epoch-rang, etend=epoch+rang;
    spkopn_c(filename.c_str(),ifname.c_str(),ncomch,&handle);
    spkw05_c(handle,objid,cntrid,cframe.c_str(),etbeg,etend,segmid.c_str(),
            cntrgm,nstate,cstate,cepoch);
    spkcls_c(handle);
}

// Make a single-segment Type 03 (Double Chebyshev) SPK file
void write_spk3_c(const std::string &filename,int body,int center,const std::string &frame,
        double first,double last,double intlen,unsigned n,unsigned polydg,
        ConstSpiceDouble *cdata, double btime) {
    std::string ifname="type 3 spk file", segid="1";
    int ncomch=5000,handle;
    spkopn_c(filename.c_str(),ifname.c_str(),ncomch,&handle);
    spkw03_c(handle,body,center,frame.c_str(),first,last,segid.c_str(),intlen,n,polydg,cdata,btime);
    spkcls_c(handle);
}

// Make a single-segment Type 02 (Single Chebyshev) SPK file
void write_spk2_c(const std::string &filename,int body,int center,const std::string &frame,
        double first,double last,double intlen,unsigned n,unsigned polydg,
        ConstSpiceDouble *cdata, double btime) {
    std::string ifname="type 2 spk file", segid="1";
    int ncomch=5000,handle;
    spkopn_c(filename.c_str(),ifname.c_str(),ncomch,&handle);
    spkw02_c(handle,body,center,frame.c_str(),first,last,segid.c_str(),intlen,n,polydg,cdata,btime);
    spkcls_c(handle);
}

// Make a single-segment Type 09 (Lagrange Interpolation) SPK file
void write_spk9_c(const std::string &filename,int body,int center,const std::string frame,
        double first,double last,int n,int degree,double *states,double *epochs) {
    std::string ifname="type 9 spk file", segid="1";
    int ncomch=5000, handle;
    spkopn_c(filename.c_str(),ifname.c_str(),ncomch,&handle);
    spkw09_c(handle,body,center,frame.c_str(),first,last,segid.c_str(),degree,n,states,epochs);
    spkcls_c(handle);
}
