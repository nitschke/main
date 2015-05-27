#ifndef DEFECTTRACKER_H
#define DEFECTTRACKER_H

#include "AMDiS.h"
#include "BP_pFieldOnSurf.h"

class DefectTracker
{
public:
	DefectTracker( DOFVector<double> *phi_ );
	
	~DefectTracker();
	
	void track( AdaptInfo *adaptInfo , BP_pFieldOnSurf *problem, int k = -1 );
    
private:
	int mode , enabled;
	double startingTime ;
	double actualTime ;
	DOFVector<double> *phi ;
	std::string folder;
	
	double radiusDefectTypeIntegration;
	
protected:
	void trackDefects( AdaptInfo *adaptInfo );	
	void trackDefectsAndTypes( AdaptInfo *adaptInfo, BP_pFieldOnSurf *problem);	
};


#endif
