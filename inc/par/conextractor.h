#ifndef _CONTOUREXTRACTOR_H
#define _CONTOUREXTRACTOR_H

#include "surface3d.h"
#include "reg3data.h"

class ContourExtractor {
public:    
	
    ContourExtractor() {}
	
    virtual ~ContourExtractor() {}
	
    virtual Surface3D* contourReg3Data(Reg3Data* p_reg3d, float isoval) =0;

private:    

};
#endif //_CONTOUREXTRACTOR_H

