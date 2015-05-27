#ifndef EXTREMEVALUES_H
#define EXTREMEVALUES_H

#include "VectorOperations.h"
#include "Helpers.h"
#include "GenericOperatorTerm.h"

using namespace AMDiS;

/// functor for min assign operator
template<typename T>
struct min_assign : FunctorBase
{
  typedef T result_type;
  
  static result_type& apply(T& v, T const& v0) { return (v = std::min(v,v0) ); }
  result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
};

  
inline void unique_value(std::vector<WorldVector<double> > &vec, std::vector<double> &val, double tol, std::vector<unsigned> *ind = NULL)
{
	compareTol<WorldVector<double> > comp(tol);
	size_t newVec = 0;
	for (size_t i = 0; i < vec.size(); ++i) 
	{
		bool inNew = false;
		for (size_t j = 0; j < newVec; ++j) 
		{
			if (comp(vec[i], vec[j])) 
			{
				if (val[i] > val[j]+FLT_TOL) 
				{
					vec[j] = vec[i];
					val[j] = val[i];
				}
				inNew = true;
				break;
			}
		}
		if (!inNew) 
		{
			vec[newVec] = vec[i]; 
			val[newVec] = val[i]; 
			newVec++; 
			if (ind)
				ind->push_back(i);
		}
	}
	vec.erase(vec.begin()+newVec,vec.end());
}

/**
*	Traverses local Mesh to identify maximas, in parallel extremas on interior boundarys(between ranks) are removed if they only exist on single rank
*	results are transfered to rank 0 => only rank 0 has COMPLETE set of extremas
*/
inline int maxima(const DOFVector<double> *rho, std::vector<WorldVector<double> >&x, std::vector<double> *values = NULL)
{ 
	FUNCNAME("maxima()");

	int num= 0;
	double tf= 0.7, tolF= 2.0;
	Parameters::get("user parameter->VortexTracker->Maxima->threshold factor",tf);
	Parameters::get("user parameter->VortexTracker->Maxima->maxima tolerance",tolF);

	double minH, maxH;
	int minLevel, maxLevel;
	Helpers::calcMeshSizes(rho->getFeSpace()->getMesh(), minH, maxH, minLevel, maxLevel) ;

	double h= 2.0*minH;
	double tol= tolF*h;

	DOFVector<int> selectDOF(rho->getFeSpace(), "selection");
	selectDOF.set(1); // possibly unnecessary when using standard-value of DOFVector

	double minNorm= min(*rho), maxNorm= max(*rho);
	#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	Parallel::mpi::globalMin(minNorm);
	Parallel::mpi::globalMax(maxNorm);
	#endif
	double threshold= (1.0-tf) * minNorm + tf * maxNorm;

	const DOFAdmin* admin= rho->getFeSpace()->getAdmin();
	const BasisFunction *basFcts= rho->getFeSpace()->getBasisFcts();
	int nBasFcts= basFcts->getNumber();
	std::vector<DegreeOfFreedom> localIndices(nBasFcts);

	///
	/// Traverse Local Mesh and mark( 0/1 ) Dofs as Extremas
	///
	TraverseStack stack;
	ElInfo *elInfo= stack.traverseFirst( rho->getFeSpace()->getMesh() , -1 , Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS );
	while (elInfo) 
	{
		basFcts->getLocalIndices( elInfo->getElement() , admin , localIndices );
		double maxValue= -1.e15;
		int maxIdx= -1;
		// find vertex with maximal value that fulfills some conditions for each element
		FixVec<WorldVector<double>,VERTEX> coords = elInfo->getMacroElement()->getCoord();
		for (int i= 0; i < nBasFcts; ++i) {
			double tempValue= (*rho)[localIndices[i]];
			if (tempValue > std::max(maxValue, threshold)) {
				maxIdx= i;
				maxValue= tempValue;
			}
		}
		// mark other vertices than minimal vertex, to be definitiv no minima
		for (int i= 0; i < nBasFcts; ++i) {
			if (i != maxIdx) {
				selectDOF[localIndices[i]]= 0;
			}
		}
		elInfo = stack.traverseNext(elInfo);
	}

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	///
	/// In parallel Case check Extremas on internal Boundarys(boundary Extrema)
	Parallel::MeshDistributor::globalMeshDistributor->synchVector(selectDOF,min_assign<int>());
	
#endif	
	
	/// 
	/// get Positions, Values and DOFindex of extremal Dofs on local Mesh
	///
	x.clear();
	if (values)
	values->clear();

	std::vector<DegreeOfFreedom> dofs;

	DOFVector<WorldVector<double> > coords(rho->getFeSpace(), "coords");
	rho->getFeSpace()->getMesh()->getDofIndexCoords(coords);

	DOFIterator<WorldVector<double> > coordsIter(&coords, USED_DOFS);
	DOFIterator<int> selectIter(&selectDOF, USED_DOFS);
	DOFConstIterator<double> energyIter(rho, USED_DOFS);

	for ( coordsIter.reset() , selectIter.reset() , energyIter.reset() ; !coordsIter.end() ; ++coordsIter , ++selectIter , ++energyIter )
	{
		if (*selectIter == 1) 
		{
			bool alreadyFound= false;
			// filter minima that are too close to an other minima
			for ( std::vector<WorldVector<double> >::iterator it= x.begin() ; it != x.end() ; ++it ) 
			{
				if (norm(*it - (*coordsIter))<tol) 
				{ 
					alreadyFound= true; 
					break; 
				}
			}
			if (!alreadyFound) 
			{
				x.push_back(*coordsIter);
				dofs.push_back(coordsIter.getDOFIndex());

				if (values)
					values->push_back(*energyIter);
			}
		}
	}
	
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS	
	///send boundary extrema list to rank 0
	AMDiS::Parallel::StdMpi< std::vector<WorldVector<double> > > coordCom(MPI::COMM_WORLD);
	AMDiS::Parallel::StdMpi< std::vector<double > > valueCom(MPI::COMM_WORLD);

	if ( MPI::COMM_WORLD.Get_rank() == 0 )
	{ 
		//recieve data on rank 0
		for(int r=1; r<MPI::COMM_WORLD.Get_size(); r++)
		{
			coordCom.recv(r);
			if (values)
			  valueCom.recv(r);
		}		
	}
	else
	{
		//send data to rank 0
		coordCom.send(0,x);
		if (values)
		  valueCom.send(0,(*values) );
	}	
	coordCom.startCommunication();
	valueCom.startCommunication();
	
	//add extremas to list at rank 0
	if ( MPI::COMM_WORLD.Get_rank() == 0 )
	{
	  //add recieved data to lists
	  std::vector<WorldVector<double> > recieveCoords;
	  std::vector< double > recieveValue;
	  for(int r=1; r<MPI::COMM_WORLD.Get_size(); r++)
	  {
		  recieveCoords = coordCom.getRecvData(r);
		  if (values)
		    recieveValue = valueCom.getRecvData(r);

		  for (int m=0; m<recieveCoords.size(); m++){
		    x.push_back(recieveCoords[m]);
		    if (values)
		      values->push_back(recieveValue[m]);
		  }
	  }
	}
#endif
	//remove multiple entries from list
	unique_value(x,*values,1e-4);

	num = x.size();
	
	MSG("found %d values on rank\n",num);
	for(int d=0; d<num; d++)
	  MSG("\t defect %d : %f %f %f \n",d,x[d][0],x[d][1],x[d][2]);
	

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	/// send final list of defects from Rank 0 to all ranks
	AMDiS::Parallel::StdMpi< std::vector<WorldVector<double> > > coordComFinal(MPI::COMM_WORLD);
	AMDiS::Parallel::StdMpi< std::vector<double > > valueComFinal(MPI::COMM_WORLD);

	if ( MPI::COMM_WORLD.Get_rank() == 0 )
	{ 
	  //send values to all ranks
	  for(int r=1; r<MPI::COMM_WORLD.Get_size(); r++)
	    {
	      coordComFinal.send(r,x);
	      if (values)
		valueComFinal.send(r,(*values) );
	    }		
	}
	else
	{
	  coordComFinal.recv(0);
	  if (values)
	    valueComFinal.recv(0);
	}	
	coordComFinal.startCommunication();
	valueComFinal.startCommunication();
	
	if ( MPI::COMM_WORLD.Get_rank() != 0 )
	{
	  ///get final values to return valiables
	  std::vector<WorldVector<double> > recieveCoordsFinal;
	  std::vector< double > recieveValueFinal;
	  
	  recieveCoordsFinal = coordComFinal.getRecvData(0);
	  if (values)
	    recieveValueFinal = valueComFinal.getRecvData(0);
	  //clear work results
	  x.clear();
	  if (values)
	    values->clear();
	  
	  //copy recieved values to output variables
	  for(int d=0; d<recieveCoordsFinal.size(); d++)
	  {
	    x.push_back(recieveCoordsFinal[d]);
	    if (values)
	      values->push_back(recieveValueFinal[d]);
	  }
	  
	  num = x.size();
	}
	
#endif	
	MSG("final List %d values overall\n",num);
	for(int d=0; d<num; d++)
	  MSG("\t defect %d : %f %f %f \n",d,x[d][0],x[d][1],x[d][2]);
	
	return num;
}
  
inline int minima( DOFVector<double> *rho , std::vector<WorldVector<double> >&x , std::vector<double> *values = NULL )
{
	FUNCNAME("minima()");
	DOFVector<double> minusRho( rho->getFeSpace() , "minusRho" );
	minusRho << -1.0 * valueOf(rho) ;
	return maxima( &minusRho , x , values ) ;
}
  
#endif
