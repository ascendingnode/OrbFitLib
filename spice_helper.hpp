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
    for(int i=0;i<card_c(&ids);i++)
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
void write_spk5_c(const std::string &filename,int objid,int cntrid,double epoch,
        double cntrgm,const std::vector<double> &state,double rang,
        std::string cframe) {
    std::string ifname="simple spk file", segmid="1";
    int ncomch=5000, nstate=1, handle;
    double cstate[1][6]; for(unsigned i=0;i<6;i++) cstate[0][i] = state[i];
    double cepoch[1] = {epoch}, etbeg=epoch-rang, etend=epoch+rang;
    spkopn_c(filename.c_str(),ifname.c_str(),ncomch,&handle);
    spkw05_c(handle,objid,cntrid,cframe.c_str(),etbeg,etend,segmid.c_str(),
            cntrgm,nstate,cstate,cepoch);
    spkcls_c(handle);
}
