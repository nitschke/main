#ifndef NORMALSAPPROXIMATOR_H
#define NORMALSAPPROXIMATOR_H

#include "AMDiS.h"
#include "kdtree_nanoflann_dof.h"
#include "QubicGeometryQuantities.h"
#include "QuadricGeometryQuantities.h"
#include "MyHelpers.h"


using namespace std;
using namespace mtl;
using namespace AMDiS;
using namespace itl;

class NormalsApproximator
{
public:
	NormalsApproximator( DOFVector<WorldVector<double> > *normals_ , WorldVector<DOFVector<double>* >* normalsX_ , Projection *proj_ = NULL ) : normals(normals_) ,normalsX(normalsX_), proj(proj_)
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
		TEST_EXIT( gradN0 && gradN1 && gradN2 && laplaceN && meanCurvature && gaussianCurvature )( "Forgot to call NormalsApproximator::setData() or change calculation mode in init file!" );
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
		using namespace experimental;
		
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
		using namespace experimental;
		
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
		using namespace experimental;
		
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
