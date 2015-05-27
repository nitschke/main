#include "DefectTracker.h"

#include "kdtree_nanoflann_dof.h"
#include "ExtremeValues.h"
#include "Expressions.h"


DefectTracker::DefectTracker( DOFVector<double> *phi_ ) : phi(phi_)
{
	startingTime = -1.0 ;
	Initfile::get( "user parameter->DefectTracker->start time" , startingTime );
	TEST_EXIT( startingTime != -1.0 )( "startingTime not set (user parameter->DefectTracker->start time)!" );
	
	mode = -1 ;
	Initfile::get( "user parameter->DefectTracker->mode" ,mode );
	TEST_EXIT( mode != -1 )( "user parameter->DefectTracker->mode not set!\n use > 1: provided scalar field \n" );
	
	if ( !phi )
		ERROR_EXIT( "DOFVector is NULL!\n" );
	
	folder = "";
	Initfile::get( "user parameter->output->folder" , folder );
	
	enabled = 0;
	Initfile::get( "user parameter->DefectTracker->enabled" , enabled );
	
	actualTime = 0.0 ;
	
	radiusDefectTypeIntegration = 0.0;
	Initfile::get( "user parameter->DefectTracker->radiusDefectTypeIntegration" , radiusDefectTypeIntegration );
}

DefectTracker::~DefectTracker() {}

void DefectTracker::track( AdaptInfo *adaptInfo , BP_pFieldOnSurf *problem, int k)
{
	FUNCNAME("DefectTracker::track()");
	actualTime = adaptInfo->getTime() ;
	if ( enabled = 1 && adaptInfo->getTime() >= startingTime )
	{
		MSG( " > TRACK DEFECTS...\n" );
		int i = ( k == -1 ? mode : k ) ;
		switch ( i ) 
		{
			case 1 :
				trackDefects(adaptInfo);
				break;
			case 2 :
				trackDefectsAndTypes(adaptInfo,problem);
				break;
				
			default :
				ERROR_EXIT("Wrong mode for DefectTracker!");
				break;
		}
		MSG( "    > done\n" );
	}
}
   
void DefectTracker::trackDefects( AdaptInfo *adaptInfo )
{
	FUNCNAME( "DefectTracker::trackDefects()" );
	
	std::vector<WorldVector<double> >defectLocations ;
	std::vector<double> defectValues ;
			
	// DETECT MINIMA
	int nDefect = maxima( phi , defectLocations , &defectValues ) ;
	
	// WRITE CURRENT VORTICES TO A FILE IN XYZ-FORMAT (USED IN MATLAB)
	stringstream timeStringTmp ;
	timeStringTmp << adaptInfo->getTime() ;
	string timeString = timeStringTmp.str();
	std::string strDefectLocations = folder + "/DefectTracker/defectLocations_" + timeString + ".dat";
	#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	if ( MPI::COMM_WORLD.Get_rank() == 0 )
	{
	#endif
		// WRITE MINIMA INFORMATION IN XYZ FORMAT
		ofstream foutDefectLocations( strDefectLocations.c_str() , ios::app );
		foutDefectLocations << defectLocations.size() << endl;
		foutDefectLocations << "time: " << timeString << endl;
		std::vector<WorldVector<double> >::iterator miniIter;
		int i = 0 ;
		for ( miniIter = defectLocations.begin() ; miniIter != defectLocations.end() ; ++i , ++miniIter )
			foutDefectLocations << i << " " << (*miniIter)[0] << " " << (*miniIter)[1] << " " << (*miniIter)[2] << " " <<defectValues[i] << endl;
	#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	}
	#endif
}

void DefectTracker::trackDefectsAndTypes( AdaptInfo *adaptInfo, BP_pFieldOnSurf *problem)
{
	FUNCNAME( "DefectTracker::trackDefectsAndTypes()" );
	
	std::vector<WorldVector<double> >defectLocations ;
	std::vector<double> defectValues ;
			
	// DETECT MINIMA
	int nDefect = maxima( phi , defectLocations , &defectValues ) ;
	
	//determine Type of Defects via integration over divergence of p*doubleWell
	std::vector<double> defectTypes; defectTypes.resize(nDefect); 
	for(int n=0; n<nDefect; n++)
	  defectTypes[n]=0.0;
	
	if(radiusDefectTypeIntegration > 0.0)
	{
	  
	  DOFVector<double> *phaseToDefect = new DOFVector<double>(phi->getFeSpace(), "distToDefect" );
	  DOFVector<double> *divP = new DOFVector<double>(phi->getFeSpace(), "divP" );
	  
	  for(int n=0; n<nDefect; n++)
	  {
	    MSG("get defect Type for %d defect \n",n);
	    //get phase function for defect vicinity
	    (*phaseToDefect) << 0.5*tanh( (two_norm(defectLocations[n] - X()) - radiusDefectTypeIntegration)/(radiusDefectTypeIntegration/5.0) );
	    
	    //get restricted divergence p
	    (*divP) << valueOf(problem->getDoubleWell()) * valueOf(phaseToDefect) * 
			( derivativeOf( problem->getProblem()->getSolution(0), 0) 
			+ derivativeOf( problem->getProblem()->getSolution(1), 1) 
			+ derivativeOf( problem->getProblem()->getSolution(2), 2) );
			
	    //get integral
	    defectTypes[n] = divP->Int();
	   
	    //clean up
	    divP->set(0.0);
	    phaseToDefect->set(0.0);
	  }
	  
	  delete phaseToDefect,divP;
	  
	}
	else
	{
	  MSG("provided radiusDefectTypeIntegration=0.0 -> types will not be dertmined \n");
	}
			
	// WRITE CURRENT VORTICES TO A FILE IN XYZ-FORMAT (USED IN MATLAB)
	stringstream timeStringTmp ;
	timeStringTmp << adaptInfo->getTime() ;
	string timeString = timeStringTmp.str();
	std::string strDefectLocations = folder + "/DefectTracker/defectLocations_" + timeString + ".dat";
	#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	if ( MPI::COMM_WORLD.Get_rank() == 0 )
	{
	#endif
		// WRITE MINIMA INFORMATION IN XYZ FORMAT
		ofstream foutDefectLocations( strDefectLocations.c_str() , ios::app );
		foutDefectLocations << defectLocations.size() << endl;
		foutDefectLocations << "time: " << timeString << endl;
		std::vector<WorldVector<double> >::iterator miniIter;
		int i = 0 ;
		for ( miniIter = defectLocations.begin() ; miniIter != defectLocations.end() ; ++i , ++miniIter )
			foutDefectLocations << i << " " << (*miniIter)[0] << " " << (*miniIter)[1] << " " << (*miniIter)[2] 
					      << " " << defectValues[i] << " " << defectTypes[i] << endl;
	#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	}
	#endif
}

