#ifndef NORMALSAPPROXIMATOR_H
#define NORMALSAPPROXIMATOR_H

#include "AMDiS.h"
#include "kdtree_nanoflann_dof.h"
#include "QubicGeometryQuantities.h"
#include "QuadricGeometryQuantities.h"
#include "MyHelpers.h"
#include "MeshHelper.h"


using namespace std;
using namespace mtl;
using namespace AMDiS;
using namespace itl;

class NormalsApproximator
{
public:
	NormalsApproximator( DOFVector<WorldVector<double> > *normals_ , 
                        WorldVector<DOFVector<double>* >* normalsX_ , 
                        DOFVector<double> *conditions_ = NULL, 
                        Projection *proj_ = NULL ) : normals(normals_) ,normalsX(normalsX_), proj(proj_), conditions(conditions_)
	{
		mode = 3 ;
		Initfile::get( "user parameter->NormalsApproximator->mode" , mode );
		weighting = false ;
		Initfile::get( "user parameter->NormalsApproximator->weighting" , weighting );
		weightingFactor = -99.0 ;
		Initfile::get( "user parameter->NormalsApproximator->weighting->factor" , weightingFactor );
		TEST_EXIT( weightingFactor != -99.0 )( "user parameter->NormalsApproximator->weighting->factor not set!" );

		// BUILD STRUCTURE FOR COEFFICIENTS OF APPROXIMATED NORMALS
		coefficientDOFVector = new DOFVector<mtl::dense_vector<double> >( normals->getFeSpace() , "coefficientDOFVector" );
		coords = new DOFVector<WorldVector<double> >( normals->getFeSpace() , "coords" );

	}
	
	/// DESTRUCTOR
	~NormalsApproximator() 
	{
		delete coefficientDOFVector;
		delete coords;
	}

	void fillNormals()
	{
		fillNormalsByMode( mode ) ;
	}


private:
  DOFVector<double> *conditions;

	DOFVector<mtl::dense_vector<double> > *coefficientDOFVector;
	DOFVector<WorldVector<double> > *normals , *coords;
	WorldVector<DOFVector<double>* > *normalsX ;
	Projection *proj ;
	int mode ;
	bool weighting ;
	double weightingFactor;
	
	void fillNormalsByMode( int i )
	{
		FUNCNAME("NormalsApproximator::fillNormals()\n");
		//TEST_EXIT( gradN0 && gradN1 && gradN2 && laplaceN && meanCurvature && gaussianCurvature )( "Forgot to call NormalsApproximator::setData() or change calculation mode in init file!" );
		switch (i)
		{
			case 0 :
				if ( proj != NULL )
					fillNormalsTraverse();
				else
					fillNormalsQuadric();
				break;
			case 1 :
				fillNormalsQuadric();
				break;
			case 2 :
				fillNormalsWeighted();
				break;
			case 3 :
				fillNormalsWeighted();
				break;
			case 4 :
				fillNormalsQubic();
				break;
			case 5 :
				fillNormalsQubicExt();
				break;
			default :
				fillNormalsWeighted();
				break;
		}
	}

	/*
	 * UPDATES THE COORDS
	 */
	void updateCoords()
	{
		normals->getFeSpace()->getMesh()->getDofIndexCoords( *coords );
	}

	/*
	 * CALCULATE NORMALS ACCORDING TO AN APPROXIMATION BY USING A QUBIC OF THE FORM
	 *     \sum_{i,j,k=0}^3 a_{ijk} * x_i * x_j * x_k + c = 0
	 * HERE, THE POINTS ARROUND EACH GRID POINT IS A RESULT OF THE PROJECTION AND MUST
	 * NOT BE IN THE SAME PARTITION.
	 */
	void fillNormalsTraverse()
	{
		FUNCNAME( "NormalsApproximator::fillNormalsTraverse()" );
		MSG("NormalsApproximator::fillNormalsTraverse()\n");
		using namespace mtl;
		using namespace itl;
		//using namespace experimental;
		
		DOFVector<WorldVector<double> > oldNormals( normals->getFeSpace() , "oldNormals" ) ;
		oldNormals = *normals ;
		DOFVector<int> normalsCalc( normals->getFeSpace() , "normalsCalc" ) ;
		normalsCalc.set( 0 );
		
		Mesh *mesh = normals->getFeSpace()->getMesh();
		TraverseStack stack;
		ElInfo *elInfo = stack.traverseFirst(mesh , -1 , Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS );
		int nPreDOFs = normals->getFeSpace()->getAdmin()->getNumberOfPreDofs(0);
		while (elInfo)
		{
			double det1 = 1.0/elInfo->getDet();
			const DegreeOfFreedom **dof = elInfo->getElement()->getDof();
			WorldVector<WorldVector<double> > edges ;
			for ( int i = 0 ; i < 3 ; i++)
				edges[i] = elInfo->getCoord( i ) ;

			/// RESULT --> NORMALS
			DegreeOfFreedom dofIndex;
			dofIndex = dof[0][nPreDOFs];
			if ( normalsCalc[dofIndex] == 0 )
			{
				(*normals)[dofIndex] = approxNormal( edges[0] , edges[1] , edges[2] ) ;
				normalsCalc[dofIndex] = 1 ;
			}
			dofIndex = dof[1][nPreDOFs];
			if ( normalsCalc[dofIndex] == 0 )
			{
				(*normals)[dofIndex] = approxNormal( edges[1] , edges[2] , edges[0] ) ;
				normalsCalc[dofIndex] = 1 ;
			}
			dofIndex = dof[2][nPreDOFs];
			if ( normalsCalc[dofIndex] == 0 )
			{
				(*normals)[dofIndex] = approxNormal( edges[2] , edges[0] , edges[1] ) ;
				normalsCalc[dofIndex] = 1 ;
			}
			elInfo = stack.traverseNext(elInfo);
		}

		// NORM NORMALS
		DOFIterator<WorldVector<double> > oldNormalsIter( &oldNormals , USED_DOFS ) ; oldNormalsIter.reset();
		DOFIterator<WorldVector<double> > normalsIter( normals , USED_DOFS ) ; normalsIter.reset();
		for (  ; !normalsIter.end() , !oldNormalsIter.end() ; ++normalsIter , ++oldNormalsIter )
		{
			double nrm = std::sqrt( (*normalsIter) * (*normalsIter) ) ;
			//     double diff = norm(*oldNormalsIter - *normalsIter);
			double signFactor = 1.0;
			//     if (diff > nrm)
			//       signFactor = -1.0;
			*normalsIter *= signFactor / nrm ;
		}
		transform(normals,normalsX);
		MSG("NormalsApproximator::fillNormalsTraverse() done\n");
	}
	
	/*
	 * CALCULATE NORMALS ACCORDING TO AN APPROXIMATION BY USING A QUADRIC OF THE FORM
	 *     \sum_{i,j=0}^3 a_{ij} * x_i * x_j + c = 0
	 * HERE, THE POINTS AROUND EACH GRID POINT IS SEARCHED BY USING A KDTREE STRUCTURE.
	 * AT THE END OTHER GEOMETRIC QUANTITIES LIKE
	 *     - \nabla_\Gamma n0
	 *     - \nabla_\Gamma n1
	 *     - \nabla_\Gamma n2
	 *     - \Delta_\Gamma n0
	 *     - Mean curvature
	 *     - Gaussian curvatue
	 * ARE FILLED ACCORDING TO THIS APPROXIMATION
	 */
	void fillNormalsQuadric()
	{
		FUNCNAME( "NormalsApproximator::fillNormalsQuadric()" );
		
		const FiniteElemSpace *feSpace = normals->getFeSpace();
		MSG("NormalsApproximator::fillNormalsQuadric()\n");
		using namespace mtl;
		using namespace itl;
		//using namespace experimental;
		
		fillNormalsWeighted();
		DOFVector<WorldVector<double> > oldNormals( normals->getFeSpace() , "oldNormals" ) ;
		oldNormals = *normals ;
		
		int nBasisFct = 10 ;
		int additional_points = 90;
		int degree = 2;
		Initfile::get( "user parameter->NormalsApproximator->additional_points" , additional_points );
		
		// UPDATE COORDS DOF
		updateCoords();
		
		DOFIterator<WorldVector<double> > coordsIter( coords, USED_DOFS);
		DOFIterator<WorldVector<double> > normalsIter(normals, USED_DOFS);
		DOFIterator<mtl::dense_vector<double> > coeffIter( coefficientDOFVector , USED_DOFS );
		coordsIter.reset();
		normalsIter.reset();
		coeffIter.reset();
		
		KD_Tree_Dof *tree = provideKDTree( feSpace );
		tree->reinit();
		
		for (  ; !coordsIter.end() , !normalsIter.end() , !coeffIter.end() ; ++coordsIter , ++normalsIter , ++coeffIter )
		{
			
			// SEARCH FOR EVERY COORD THE numPoints NEAREST POINTS IN COORDS SET (DEPENDS ON MESH REFINEMENT LEVEL)
			WorldVector<double> p = (*coordsIter) ;
			size_t numPoints = nBasisFct + additional_points ; // 10 basisFunctions and 40 additional Points
			
			std::vector<DegreeOfFreedom> indices(numPoints);
			std::vector<double> distances(numPoints);
			tree->query(p.begin(), numPoints, &indices[0], &distances[0]);
			
			// CALCULATE MAX. DISTANCES FOR WEIGHTING AND BUILD WEIGHTING VECTOR
			double dMax = 0.0 ;
			for ( size_t j = 0; j < indices.size(); ++j )
			{
				double d = std::sqrt( ( tree->m_data[indices[j]] - p ) * ( tree->m_data[indices[j]] - p ) );
				if ( d > dMax )
					dMax = d ;
			}
			
			std::vector<WorldVector<double> > points;
			std::vector<double> weights;
			for ( size_t j = 0; j < indices.size(); ++j ) 
			{
				points.push_back( tree->m_data[indices[j]] );
				double d = std::sqrt( ( tree->m_data[indices[j]] - p ) * ( tree->m_data[indices[j]] - p ) ) ;
				weights.push_back( pow( weightingFactor , (d/dMax)*(d/dMax) ) );
			}
			
			// PROVIDE MATRICES FOR LEAST SQUARE APPROACH
			TEST_EXIT( points.size() >= numPoints )( "Nicht genügend Stützstellen gefunden!\n" );
			size_t i = 0;
			mtl::dense2D<double> A(numPoints,6), B(6,6);
			mtl::dense_vector<double> b(numPoints, 1.0), c(6), y(6, 1.0);
			std::vector<WorldVector<double> >::const_iterator pointsIter;
			std::vector<double>::const_iterator weightsIter;
			for ( pointsIter = points.begin() , weightsIter = weights.begin() ; pointsIter != points.end() && i < numPoints ,  weightsIter != weights.end() ; ++pointsIter , ++weightsIter , ++i ) 
			{
				WorldVector<double> x = *pointsIter ;
				double w = 1.0 ;
				if ( weighting )
					w = *weightsIter ;
				A[i][0] = w * sqr(x[0]) ; 
				A[i][1] = w * sqr(x[1]) ; 
				A[i][2] = w * sqr(x[2]) ;
				A[i][3] = w * 2.0 * x[0]*x[1] ; 
				A[i][4] = w * 2.0 * x[0]*x[2] ; 
				A[i][5] = w * 2.0 * x[1]*x[2] ;
				b[i] = w ;
			}
			
			// SOLVE LEAST SQUARE PROBLEM
			//	y[0] = a[0][0] ;
			//	y[1] = a[1][1] ;
			//	y[2] = a[2][2] ;
			//	y[3] = a[0][1] ;
			//	y[4] = a[0][2] ;
			//	y[5] = a[1][2] ;
			c = trans(A)*b;
			B = trans(A)*A;
			solveLGS( B , c , y );
			
			// UPDATE COEFFICIENT DOF VECTOR
			*coeffIter = y ;
			
			// SAVE OLD NORMAL
 			WorldVector<double> oldNormal = (*normalsIter) ;
			
			// BUILD NORMALS
			(*normalsIter)[0] = QuadricGeometryQuantities::getN0( y , p );
			(*normalsIter)[1] = QuadricGeometryQuantities::getN1( y , p );
			(*normalsIter)[2] = QuadricGeometryQuantities::getN2( y , p );
			
		}
		
		/// SYNCH NORMALS VECTOR AND NORMALIZE AFTERWARDS -> NEEDED, BECAUSE SYNCHRONIZATION CANNOT
		/// ASSUME THAT  THE NEW VALUE AT INTERANL BOUNDARIES IS OF VALUE 1
		#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
		MyHelpers::mySynchAddVector( *normals );
		for ( normalsIter.reset() ; !normalsIter.end() ; ++normalsIter )
			*normalsIter *= 1.0 / std::sqrt( (*normalsIter) * (*normalsIter) ) ;
		#endif
		transform(normals,normalsX);

		// OUTPUT
		MSG("NormalsApproximator::fillNormalsQuadric() done\n");
	}
	

  /*
	 * CALCULATE NORMALS ACCORDING TO AN APPROXIMATION BY USING A QUBIC OF THE FORM
	 *     \sum_{i,j,k=0}^3 a_{ijk} * x_i * x_j * x_k + c = 0
	 * HERE, THE POINTS AROUND EACH GRID POINT IS SEARCHED BY USING A KDTREE STRUCTURE
	 * AT THE END OTHER GEOMETRIC QUANTITIES LIKE
	 *     - \nabla_\Gamma n0
	 *     - \nabla_\Gamma n1
	 *     - \nabla_\Gamma n2
	 *     - \Delta_\Gamma n
	 *     - Mean curvature
	 *     - Gaussian curvatue
	 * ARE FILLED ACCORDING TO THIS APPROXIMATION
	 */
	void fillNormalsQubic()
	{
		FUNCNAME( "NormalsApproximator::fillNormalsQubic()" );
		const FiniteElemSpace *feSpace = normals->getFeSpace();
		MSG("NormalsApproximator::fillNormalsQubic()\n");
		using namespace mtl;
		using namespace itl;
		//using namespace experimental;
		
		fillNormalsWeighted();
		
		int nBasisFct = 10 ;
		int additional_points = 18;
		int degree = 2;
		Initfile::get( "user parameter->NormalsApproximator->additional_points" , additional_points );
		
		// CREATE COORDS DOF
		updateCoords();
		
		// DOF ITERATORS
		DOFIterator<WorldVector<double> > coordsIter(coords, USED_DOFS);
		DOFIterator<WorldVector<double> > normalsIter(normals, USED_DOFS);
		DOFIterator<mtl::dense_vector<double> > coeffIter( coefficientDOFVector , USED_DOFS );
		coordsIter.reset();
		normalsIter.reset();
		coeffIter.reset();
		
		// PROVIDE KD TREE
		KD_Tree_Dof *tree = provideKDTree( feSpace );
		tree->reinit();
		
		for (  ; !coordsIter.end() , !normalsIter.end() , !coeffIter.end() ; ++coordsIter , ++normalsIter , ++coeffIter )
		{
			// SEARCH FOR EVERY COORD THE numPoints NEAREST POINTS IN COORDS SET (DEPENDS ON MESH REFINEMENT LEVEL)
			WorldVector<double> p = (*coordsIter) ;
			size_t numPoints = nBasisFct + additional_points ; // 10 basisFunctions and 40 additional Points
			std::vector<DegreeOfFreedom> indices(numPoints);
			std::vector<double> distances(numPoints);
			tree->query(p.begin(), numPoints, &indices[0], &distances[0]);
			
			// CALCULATE MAX. DISTANCES FOR WEIGHTING AND BUILD WEIGHTING VECTOR
			double dMax = 0.0 ;
			for ( size_t j = 0; j < distances.size(); ++j )
				dMax = ( dMax < distances[j] ? distances[j] : dMax ) ;
			std::vector<WorldVector<double> > points;
			std::vector<double> weights;
			double r = 0.01;
			for ( size_t j = 0; j < indices.size(); ++j )
			{
				points.push_back( tree->m_data[indices[j]] );
				weights.push_back( std::pow( r , (distances[j]/dMax)*(distances[j]/dMax) ) );
			}
			
			// PROVIDE MATRICES FOR LEAST SQUARE APPROACH
			TEST_EXIT( points.size() >= numPoints )( "Nicht genügend Stützstellen gefunden!\n" );
			size_t i = 0;
			size_t dim = 10 ;
			mtl::dense2D<double> A(numPoints,dim), B(dim,dim);
			mtl::dense_vector<double> b(numPoints, 1.0), c(dim), y(dim, 1.0);
			std::vector<WorldVector<double> >::const_iterator pointsIter;
			std::vector<double>::const_iterator weightsIter;
			for ( pointsIter = points.begin() , weightsIter = weights.begin() ; pointsIter != points.end() && i < numPoints ,  weightsIter != weights.end() ; ++pointsIter , ++weightsIter , ++i )
			{
				WorldVector<double> x = *pointsIter ;
				double w = 1.0 ;
				if ( weighting )
					w = *weightsIter ;
				A[i][0] = w * x[0]*x[0]*x[0] ;
				A[i][1] = w * x[1]*x[1]*x[1] ;
				A[i][2] = w * x[2]*x[2]*x[2] ;
				A[i][3] = w * 3.0 * x[0]*x[0]*x[1] ;
				A[i][4] = w * 3.0 * x[0]*x[0]*x[2] ;
				A[i][5] = w * 3.0 * x[0]*x[1]*x[1] ;
				A[i][6] = w * 6.0 * x[0]*x[1]*x[2] ;
				A[i][7] = w * 3.0 * x[0]*x[2]*x[2] ;
				A[i][8] = w * 3.0 * x[1]*x[1]*x[2] ;
				A[i][9] = w * 3.0 * x[1]*x[2]*x[2] ;
				b[i] = w ;
			}
			
			// SOLVE LEAST SQUARE PROBLEM
			c = trans(A)*b;
			B = trans(A)*A;
			solveLGS( B , c , y );
			
			// FILL IN NORMAL VECTOR WITH THE APPROXIMATION OF TEH SURFACE BY A QUADRIK
			WorldVector<double> x = *coordsIter ;
			(*normalsIter)[0] = QubicGeometryQuantities::getN0( y , p );
			(*normalsIter)[1] = QubicGeometryQuantities::getN1( y , p );
			(*normalsIter)[2] = QubicGeometryQuantities::getN2( y , p );
			
			// UPDATE COEFFICIENT DOF VECTOR
			*coeffIter = y ;
			
		}
		
		#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
		/// SYNCH NORMALS VECTOR!!
		Parallel::MeshDistributor::globalMeshDistributor->synchVector( *normals ) ;
		#endif
		transform(normals,normalsX);
		
		MSG("NormalsApproximator::fillNormalsQubic() done\n");
	}

	void fillNormalsQubicExt()
	{
		FUNCNAME( "NormalsApproximator::fillNormalsQubic()" );
		const FiniteElemSpace *feSpace = normals->getFeSpace();
		MSG("NormalsApproximator::fillNormalsQubic()\n");
		using namespace mtl;
		using namespace itl;
		//using namespace experimental;
		
		fillNormalsWeighted();
		
		int nBasisFct = 10 ;
    nBasisFct += 6;
		int additional_points = 18;
		int degree = 2;
		Initfile::get( "user parameter->NormalsApproximator->additional_points" , additional_points );

    bool adaptive = false;
    Initfile::get( "user parameter->NormalsApproximator->adaptive" , adaptive );
    double adaptiveTol = 0.1;
    Initfile::get( "user parameter->NormalsApproximator->adaptive tol" , adaptiveTol );

    bool average = false;
    Initfile::get( "user parameter->NormalsApproximator->average" , average );
    double nAv = 0.1;
    Initfile::get( "user parameter->NormalsApproximator->average n" , nAv );


    double c2 = 1.0;

		
		// CREATE COORDS DOF
		updateCoords();
		
    // normals precalculation for orientation preserving and adaptivity
    DOFVector<WorldVector<double> > normalsSimple = AMDiS::getNormals(feSpace, true);


		// DOF ITERATORS
    DOFIterator<double> *condIter = NULL; 
    if (conditions) {
      condIter = new DOFIterator<double>(conditions, USED_DOFS); 
      condIter->reset();
    }
    DOFIterator<WorldVector<double> > normalsSimpleIter(&normalsSimple, USED_DOFS);
		DOFIterator<WorldVector<double> > coordsIter(coords, USED_DOFS);
		DOFIterator<WorldVector<double> > normalsIter(normals, USED_DOFS);
		DOFIterator<mtl::dense_vector<double> > coeffIter( coefficientDOFVector , USED_DOFS );
    normalsSimpleIter.reset();
		coordsIter.reset();
		normalsIter.reset();
		coeffIter.reset();
		
		// PROVIDE KD TREE
		KD_Tree_Dof *tree = provideKDTree( feSpace );
		tree->reinit();
		
    size_t numPointsBasic = nBasisFct + additional_points ;
    size_t numPoints = numPointsBasic;
    int iPoint = 0;
    int nTry = 0;
    std::vector<WorldVector<double> > NTmp(nAv);
		for (  ; !coordsIter.end() , !normalsIter.end() , !coeffIter.end(), !normalsSimpleIter.end() ; ++coordsIter , ++normalsIter , ++coeffIter, ++normalsSimpleIter, ++iPoint )
		{
			// SEARCH FOR EVERY COORD THE numPoints NEAREST POINTS IN COORDS SET (DEPENDS ON MESH REFINEMENT LEVEL)
			WorldVector<double> p = (*coordsIter) ;
			std::vector<DegreeOfFreedom> indices(numPoints);
			std::vector<double> distances(numPoints);
			tree->query(p.begin(), numPoints, &indices[0], &distances[0]);
			
			// CALCULATE MAX. DISTANCES FOR WEIGHTING AND BUILD WEIGHTING VECTOR
			double dMax = 0.0 ;
			for ( size_t j = 0; j < distances.size(); ++j )
				dMax = ( dMax < distances[j] ? distances[j] : dMax ) ;
			std::vector<WorldVector<double> > points;
			std::vector<double> weights;
			double r = 0.01;
			for ( size_t j = 0; j < indices.size(); ++j )
			{
				points.push_back( tree->m_data[indices[j]] );
				weights.push_back( std::pow( r , (distances[j]/dMax)*(distances[j]/dMax) ) );
			}
			
			// PROVIDE MATRICES FOR LEAST SQUARE APPROACH
			TEST_EXIT( points.size() >= numPoints )( "Nicht genügend Stützstellen gefunden!\n" );
			size_t i = 0;
			size_t dim = 10 ;
      dim += 6;
			mtl::dense2D<double> A(numPoints,dim), B(dim,dim);
			mtl::dense_vector<double> b(numPoints, 1.0), c(dim), y(dim, 1.0);
			std::vector<WorldVector<double> >::const_iterator pointsIter;
			std::vector<double>::const_iterator weightsIter;
			for ( pointsIter = points.begin() , weightsIter = weights.begin() ; pointsIter != points.end() && i < numPoints ,  weightsIter != weights.end() ; ++pointsIter , ++weightsIter , ++i )
			{
				WorldVector<double> x = *pointsIter ;
        //cout << "*** " <<  x << endl;
				double w = 1.0 ;
				if ( weighting )
					w = *weightsIter ;
				A[i][0] = w * x[0]*x[0]*x[0] ;
				A[i][1] = w * x[1]*x[1]*x[1] ;
				A[i][2] = w * x[2]*x[2]*x[2] ;
				A[i][3] = w * 3.0 * x[0]*x[0]*x[1] ;
				A[i][4] = w * 3.0 * x[0]*x[0]*x[2] ;
				A[i][5] = w * 3.0 * x[0]*x[1]*x[1] ;
				A[i][6] = w * 6.0 * x[0]*x[1]*x[2] ;
				A[i][7] = w * 3.0 * x[0]*x[2]*x[2] ;
				A[i][8] = w * 3.0 * x[1]*x[1]*x[2] ;
				A[i][9] = w * 3.0 * x[1]*x[2]*x[2] ;
        
        A[i][10] = c2 * w * 3.0 * x[0]*x[0] ;
        A[i][11] = c2 * w * 6.0 * x[0]*x[1] ;
        A[i][12] = c2 * w * 6.0 * x[0]*x[2] ;
        A[i][13] = c2 * w * 3.0 * x[1]*x[1] ;
        A[i][14] = c2 * w * 6.0 * x[1]*x[2] ;
        A[i][15] = c2 * w * 3.0 * x[2]*x[2] ;

				b[i] = w ;
			}
      //cout << endl;
			
			// SOLVE LEAST SQUARE PROBLEM
			c = trans(A)*b;
			B = trans(A)*A;
			solveLGS( B , c , y );

      //condition
      if (conditions) {
        mtl::dense_vector<double> eig(dim) ;
        eig = qr_algo(B,100);
        *(*condIter) = abs(eig[0]/eig[dim-1]);
        //cout << *(*condIter) << endl;
        ++(*condIter);
      }
			
			// FILL IN NORMAL VECTOR WITH THE APPROXIMATION OF TEH SURFACE BY A QUADRIK
			WorldVector<double> x = *coordsIter ;
			//(*normalsIter)[0] = QubicGeometryQuantities::getN0( y , p );
			//(*normalsIter)[1] = QubicGeometryQuantities::getN1( y , p );
			//(*normalsIter)[2] = QubicGeometryQuantities::getN2( y , p );
      double x2 = x[0]*x[0];
      double y2 = x[1]*x[1];
      double z2 = x[2]*x[2];
      double xy = 2.0*x[0]*x[1];
      double xz = 2.0*x[0]*x[2];
      double yz = 2.0*x[1]*x[2];
      (*normalsIter)[0] = y[0]*x2 + y[5]*y2 + y[7]*z2 + y[3]*xy + y[4]*xz + y[6]*yz;
      (*normalsIter)[1] = y[3]*x2 + y[1]*y2 + y[9]*z2 + y[5]*xy + y[6]*xz + y[8]*yz;
      (*normalsIter)[2] = y[4]*x2 + y[8]*y2 + y[2]*z2 + y[6]*xy + y[7]*xz + y[9]*yz;

      (*normalsIter)[0] += c2 * 2.0 * (y[10]*x[0] + y[11]*x[1] + y[12]*x[2]); 
      (*normalsIter)[1] += c2 * 2.0 * (y[11]*x[0] + y[13]*x[1] + y[14]*x[2]); 
      (*normalsIter)[2] += c2 * 2.0 * (y[12]*x[0] + y[14]*x[1] + y[15]*x[2]); 


      //ensure orientation
      if (AMDiS::dot(*normalsSimpleIter, *normalsIter) < 0) (*normalsIter) *= -1.0;

      //normalize
      (*normalsIter) *= 1./sqrt(AMDiS::dot(*normalsIter, *normalsIter));

      //tol test
      double res = 0.0;
      if (adaptive) {
      //res = AMDiS::norm((*normalsIter)-(*normalsSimpleIter));
      //res =  mtl::two_norm(B*y - c);
      //res =  mtl::two_norm(A*y - b) / numPoints;
      mtl::dense_vector<double> eig(dim);
      eig = qr_algo(B,2);
      res = abs(eig[0]/eig[dim-1]);
      }
      if (adaptive && res > adaptiveTol && !average) {
        numPoints += 100;
        if (nTry!= 0 && !(nTry%100))
          MSG("NormalsApproximator::fillNormalsQubic at %d: TOL(%e > %e) not reached , increase additional_points locally to %d\n",
                  iPoint, res, adaptiveTol, numPoints-nBasisFct);
        --coordsIter; --normalsIter; --coeffIter; --normalsSimpleIter; --iPoint;
        ++nTry;
        if (conditions) --(*condIter);
      } else if (average && nTry < nAv) {
        NTmp[nTry] = *normalsIter;
        numPoints += 1;
        --coordsIter; --normalsIter; --coeffIter; --normalsSimpleIter; --iPoint;
        ++nTry;
      } else {
        // UPDATE COEFFICIENT DOF VECTOR
			  *coeffIter = y ;
        numPoints = numPointsBasic;
        nTry = 0;
        if (average) {
          (*normalsIter) = vectorMedian(NTmp);
          //(*normalsIter) = vectorHuberKMean(NTmp);
          (*normalsIter) *= 1./sqrt(AMDiS::dot(*normalsIter, *normalsIter));
        }
        if (!(iPoint%1000)) cout << iPoint << endl;
      }


			
		}
		
		#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
		/// SYNCH NORMALS VECTOR!!
		Parallel::MeshDistributor::globalMeshDistributor->synchVector( *normals ) ;
		#endif
    if (normalsX) {
      MSG("NormalsApproximator::fillNormalsQubic::transform(normals,normalsX)\n");
		  transform(normals,normalsX);
      MSG("NormalsApproximator::fillNormalsQubic::transform(normals,normalsX) done\n");
    }
		
		MSG("NormalsApproximator::fillNormalsQubic() done\n");
	}

  WorldVector<double> vectorMedian(std::vector<WorldVector<double> > vecs){
    int n = vecs.size();
    for (int k = n ; k > 1; k--) {
      for (int i = 0; i < k-1; i++) {
        for (int p= 0; p < 3; p++) {
          if ((vecs[i])[p] > (vecs[i+1])[p]) {
            double tmp = (vecs[i])[p];
            (vecs[i])[p] = (vecs[i+1])[p];
            (vecs[i+1])[p] = tmp;
          }
        }
      }
    }
    return vecs[n/2];
  }

  inline WorldVector<double> vectorHuberKMean(std::vector<WorldVector<double> > vecs){
    //cout << vecs << endl;
    double k = 1.28;
    int maxIter = 10;
    double tol = 1.e-6;

    int n = vecs.size();
    WorldVector<double> median = vectorMedian(vecs);
    std::vector<WorldVector<double> > vecsTmp(n);
    for (int i = 0; i < n; ++i) 
      for (int j = 0; j < 3; ++j) 
        (vecsTmp[i])[j] = abs((vecs[i])[j] - median[j]); 
    WorldVector<double> MAD = vectorMedian(vecsTmp);
    
    int nIter = 0;
    double res;
    WorldVector<double> z;
    WorldVector<double> w;
    WorldVector<double> mean = median;
    WorldVector<double> meanOld;
    WorldVector<double> num;
    WorldVector<double> den;
    while ( (res > tol && nIter < maxIter) || !nIter) {
    //cout << mean << endl;
      for (int i = 0; i < n; ++i) {
        z = (vecs[i] - mean);
        for (int j = 0; j < 3; ++j) z[j] /= MAD[j];
        w = vectorHuberKWeight(z, k);
        num.set(0.0);
        den.set(0.0);
        for (int j = 0; j < 3; ++j) { 
          num[j] += (vecs[i])[j] * w[j];
          den[j] += w[j];
        }
      }
      meanOld = mean;
      for (int j = 0; j < 3; ++j) mean[j] = num[j] / den[j];
      res = AMDiS::norm(mean - meanOld);
      nIter++;
    }
    //cout << mean << endl << endl;
    if (nIter == maxIter) cout << res << endl;
    return mean;
  }

  inline WorldVector<double> vectorHuberKWeight(WorldVector<double> z, double k) {
    WorldVector<double> rval;
    for (int i = 0; i < 3; i++) {
      double absZi = abs(z[i]);
      rval[i] = (absZi > k) ? (k/absZi) : 1.0;
    }
    return rval;
  }
	
	/*
	 * SIMPLE CALCULATION OF THE NORMALS WITH LOCAL WEIGHTING
	 */
	void fillNormalsWeighted()
	{
		const FiniteElemSpace *feSpace = normals->getFeSpace();
		Mesh *mesh = feSpace->getMesh();
		const BasisFunction *basisFcts = feSpace->getBasisFcts();
		int nBasisFcts = basisFcts->getNumber();
		std::vector<DegreeOfFreedom> localIndices(nBasisFcts);
		WorldVector<double> normal; normal.set(0.0);
		normals->set(normal);
		TraverseStack stack;
		ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
		while (elInfo)
		{
			elInfo->getElementNormal(normal);
			Element *el = elInfo->getElement();
			double det = elInfo->getDet();
			if (det <= DBL_TOL)
				det = elInfo->calcDet();
			basisFcts->getLocalIndices(el, feSpace->getAdmin(), localIndices);
			for (int i = 0; i < nBasisFcts; i++)
				(*normals)[localIndices[i]] += normal*(1.0/det);
			elInfo = stack.traverseNext(elInfo);
		}
		DOFIterator<WorldVector<double> > normalsIter(normals, USED_DOFS); normalsIter.reset();
		for (; !normalsIter.end(); normalsIter++)
			*normalsIter *= 1.0/std::sqrt( (*normalsIter) * (*normalsIter) ) ;

		transform(normals,normalsX);
	}

	/*
	 * CALCULATES THE NORMAL OF A TRIANGLE WITH VERTICES @param edge0, @param edge1 and @param edge2
	 * BY USING THE PROJECTION
	 */
	WorldVector<double> approxNormal( WorldVector<double> edge0 , WorldVector<double> edge1 , WorldVector<double> edge2 )
	{
		FUNCNAME("approxNormal()");
		int nMu = 5 ;
		int nLambda = 5 ;
		
		/// BUILD POINTS
		std::vector<WorldVector<double> > points ;
		std::vector<double> weights ;
		std::vector<double> muVec;
		std::vector<double> lambdaVec;
		for ( int i = 0 ; i <= nMu ; ++i )
		muVec.push_back( -0.5 + i * 1.0/nMu );
		for ( int i = 0 ; i <= nLambda ; ++i )
		lambdaVec.push_back( -0.5 + i * 1.0/nLambda );
		std::vector<double>::iterator muVecIter;
		std::vector<double>::iterator lambdaVecIter;
		for ( lambdaVecIter = lambdaVec.begin() ; lambdaVecIter != lambdaVec.end() ; ++lambdaVecIter )
		for ( muVecIter = muVec.begin() ; muVecIter != muVec.end() ; ++muVecIter )
		{
			WorldVector<double> tmp = edge0 + (*lambdaVecIter) * (edge1-edge0) + (*muVecIter) * (edge2-edge0) ;
			proj->project( tmp ) ;
			points.push_back( tmp ) ;
			double weight = std::sqrt( (tmp-edge0) * (tmp-edge0) )<=DBL_TOL ? 1.e6 : 1.0/std::sqrt( (tmp-edge0)*(tmp-edge0) ) ;
			weights.push_back( weight );
		}
		
		/// CALCULATION FOR edge0
		size_t i = 0;
		int numPoints = points.size() ; 
		mtl::dense2D<double> A(numPoints,6), B(6,6);
		mtl::dense_vector<double> b(numPoints, 1.0), c(6), y(6, 1.0);
		std::vector<WorldVector<double> >::const_iterator pointsIter;
		std::vector<double>::const_iterator weigthsIter;
		for ( pointsIter = points.begin() , weigthsIter = weights.begin() ; pointsIter != points.end() , weigthsIter != weights.end() ; ++pointsIter, ++weigthsIter , ++i ) 
		{
			WorldVector<double> x = *pointsIter ;
			double w = 1.0 ;
			if ( weighting )
				w = *weigthsIter ;
			A[i][0] = w * sqr(x[0]); 
			A[i][1] = w * sqr(x[1]); 
			A[i][2] = w * sqr(x[2]);
			A[i][3] = w * x[0]*x[1]; 
			A[i][4] = w * x[0]*x[2]; 
			A[i][5] = w * x[1]*x[2];
			b[i] = w ;
		}
		
		/// BUILD MATRICES FOR LEAST SQUARE APPROACH AND SOLVE SYSTEM
		c = trans(A)*b;
		B = trans(A)*A;
		solveLGS( B , c , y );
		
		/// FILL IN NORMAL VECTOR WITH THE APPROXIMATION OF TEH SURFACE BY A QUADRIK
		WorldVector<double> result;
		result[0] = 2.0*y[0]*(edge0)[0] + y[3]*(edge0)[1] + y[4]*(edge0)[2] ;
		result[1] = 2.0*y[1]*(edge0)[1] + y[3]*(edge0)[0] + y[5]*(edge0)[2] ;
		result[2] = 2.0*y[2]*(edge0)[2] + y[4]*(edge0)[0] + y[5]*(edge0)[1] ;
		return result;
	}
	
	/*
	 * SOLVES THE LINEAR SYSTEM
	 *     A x = b
	 * BY USING A DIRECT SOLVER BASED ON LU FACTORIZATION. IF IT CANNOT BE SOLVED,
	 * THE BICGSTAB SOLVER IS USED INSTEAD.
	 */
	void solveLGS( mtl::dense2D<double> &A , mtl::dense_vector<double> &b , mtl::dense_vector<double> &x )
	{
		FUNCNAME("NormalsApproximator::solve()");
		try
		{
			x = mtl::matrix::lu_solve( A , b ) ;
		}
		catch (exception &e)
		{
			WARNING("Could not solve normal equation directly!\nUsing BiCGStab instead!\n");
			pc::identity<mtl::dense2D<double> > P(A);
			basic_iteration<double> iter(b, 1000, 1.e-10);
			bicgstab( A , x , b , P , iter );
		}
	}
	
};

#endif
