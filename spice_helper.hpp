#include "cspice/include/SpiceUsr.h"
#include <vector>
#include <string>

// Wrapper functions for SPICE's strange anachronistic SpiceCell objects

// Extract the objects in an SPK kernel
std::vector<int> spkobjects(std::string kernel) {
    SPICEINT_CELL(ids,10000);
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
