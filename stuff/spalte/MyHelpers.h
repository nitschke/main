#ifndef HELPERSSEBASTIAN_H
#define HELPERSSEBASTIAN_H

#include "AMDiS.h"
#include "ExtendedProblemStat.h"


#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/StdMpi.h"
using namespace AMDiS::Parallel;
#endif
struct MyHelpers 
{
	/// Works quite similar to the function \ref synchVector, but instead the 
	/// values of subdomain vectors are add along the boundaries.
	#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	template<typename T>
	static void mySynchAddVector(DOFVector<T> &vec )
	{
		const FiniteElemSpace *fe = vec.getFeSpace();
		int nLevels = Parallel::MeshDistributor::globalMeshDistributor->getMeshLevelData().getNumberOfLevels();
		for (int level = nLevels - 1; level >= 0; level--) 
		{
			StdMpi<std::vector<T> > stdMpi(Parallel::MeshDistributor::globalMeshDistributor->getMeshLevelData().getMpiComm(level));

			for (DofComm::Iterator it(Parallel::MeshDistributor::globalMeshDistributor->getDofComm(level).getRecvDofs(), fe);!it.end(); it.nextRank()) 
			{
				std::vector<T> dofs;
				dofs.reserve(it.getDofs().size());

				for (; !it.endDofIter(); it.nextDof())
					dofs.push_back(vec[it.getDofIndex()]);

				stdMpi.send(it.getRank(), dofs);
			}

			for (DofComm::Iterator it(Parallel::MeshDistributor::globalMeshDistributor->getDofComm(level).getSendDofs()); !it.end(); it.nextRank())
				stdMpi.recv(it.getRank());

			stdMpi.startCommunication();

			for (DofComm::Iterator it(Parallel::MeshDistributor::globalMeshDistributor->getDofComm(level).getSendDofs(), fe); !it.end(); it.nextRank())
				for (; !it.endDofIter(); it.nextDof())
					vec[it.getDofIndex()] = 0.5 * ( vec[it.getDofIndex()] + stdMpi.getRecvData(it.getRank())[it.getDofCounter()] );
		}

		Parallel::MeshDistributor::globalMeshDistributor->synchVector(vec);
	}
	#endif
	
	static double surfaceArea( Mesh *mesh ) 
	{
		FUNCNAME("Helpers::surfaceArea()");
		// TRAVERSE MESH
		double A = 0.0 ;
		TraverseStack stack;
		ElInfo *elInfo = stack.traverseFirst(mesh , -1 , Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS );
		while (elInfo)
		{
			// COORDS OF VERTICES
			std::vector<WorldVector<double> > vertices(3) ;
			for ( int i = 0 ; i < 3 ; i++)
				vertices[i] = elInfo->getCoord( i ) ;
			elInfo = stack.traverseNext(elInfo);
			WorldVector<double> u = vertices[1] - vertices[0] ;
			WorldVector<double> v = vertices[2] - vertices[0] ;
			// CROSSPRODUCT OF EDGES
			WorldVector<double> uCrossV;
			uCrossV[0] = u[1]*v[2] - u[2]*v[1];
			uCrossV[1] = u[2]*v[0] - u[0]*v[2];
			uCrossV[2] = u[0]*v[1] - u[1]*v[0];
			A += 0.5 * std::sqrt( uCrossV*uCrossV );
		}
		#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
		double ASum = 0.0 ;
		MPI_Allreduce( &A , &ASum , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
		A = ASum ;
		MPI_Barrier( MPI_COMM_WORLD );
		#endif
		return A ;
	}
	static void normalizePhaseField( DOFVector<double> *c_ , DOFVector<double> *normalizedC_ )
	{
		FUNCNAME("CHBaseProblem::normalizePhaseField()");
		// FILL NORMALIZED PHASEFIELD
		// THIS DOFVECTOR IS SIMPLY THE SAME AS THE PHASEFIELD c BUT CUTTED IN REGIONS WHERE c > 1.0 OR c < -1.0
		DOFIterator<double> cIter( c_ , USED_DOFS );
		DOFIterator<double> normCIter( normalizedC_ , USED_DOFS );
		cIter.reset();
		normCIter.reset();
		for (  ; !cIter.end() , !normCIter.end() ; ++cIter , ++normCIter )
		{
			if ( *cIter > 1.0 )
				*normCIter = 1.0 ;
			else if ( *cIter < -1.0 )
				*normCIter = -1.0 ;
			else
				*normCIter = *cIter ;
		}
	}
	/// projectToTangentSpace()
	static void projectToTangentSpace( DOFVector<WorldVector<double> > *v , DOFVector<WorldVector<double> > *normals_ )
	{
		FUNCNAME("MyCouplingBaseProblem::projectToTangentSpace()");
		DOFVector<double> helpDOF( v->getFeSpace() , "helpDOF" );
		DOFIterator<WorldVector<double> > vIter( v , USED_DOFS ); vIter.reset();
		DOFIterator<WorldVector<double> > normalsIter( normals_ , USED_DOFS ); normalsIter.reset();
		DOFIterator<double> helpDOFIter( &helpDOF , USED_DOFS ); helpDOFIter.reset();
		for (  ; !vIter.end() , !normalsIter.end() , !helpDOFIter.end() ; ++normalsIter , ++vIter ,++helpDOFIter )
			(*helpDOFIter) = (*normalsIter) * (*vIter) ;
		vIter.reset();
		normalsIter.reset();
		helpDOFIter.reset();
		for (  ; !vIter.end() , !normalsIter.end() , !helpDOFIter.end() ; ++normalsIter , ++vIter ,++helpDOFIter )
			(*vIter) = (*vIter) - ( (*helpDOFIter) * (*normalsIter) ) ;
	}
	// INITIAL CONDITION FOR \psi
	struct Psi0 : public AbstractFunction<double,WorldVector<double> >
	{
		Psi0() : AbstractFunction<double, WorldVector<double> >(6) {};
		double operator()(const WorldVector<double> &x) const { return x[2] ; }
	};

	// INITIAL CONDITION FOR \varphi
	struct Phi0 : public AbstractFunction<double,WorldVector<double> >
	{
		Phi0() : AbstractFunction<double, WorldVector<double> >(6) {};
		double operator()(const WorldVector<double> &x) const { return -2.0 * x[2] ; }
	};

	struct SubExtendedProblemStat : public ExtendedProblemStat
	{
		SubExtendedProblemStat( const std::string &name_ , Mesh *mesh , int n ) : ExtendedProblemStat( name_ )
		{
			FUNCNAME( "SubExtendedProblemStat::SubExtendedProblemStat()" );
			meshes.resize(1);
			componentMeshes.resize( n );
			meshes[0] = mesh ;
			for (int i = 0 ; i < n ; ++i )
				setComponentMesh( i , mesh );
		}
	};

	// CLASS WIHCH REPRESENTS F(\nabla v) := n \times \grad \psi
	struct Rotation : public TertiaryAbstractFunction<WorldVector<double>,WorldVector<double>,std::vector<double>,std::vector<WorldVector<double> > >
	{
		Rotation( double factor_ = 1.0 ) : TertiaryAbstractFunction<WorldVector<double>,WorldVector<double>,std::vector<double>,std::vector<WorldVector<double> > >(6) , factor(factor_) {};
		WorldVector<double> operator()( const WorldVector<double> &x , const std::vector<double> &n , const std::vector<WorldVector<double> > &gradPsi ) const
		{
			WorldVector<double> res;
			res.set( 0.0 );
			res[0] = factor * ( n[1]*gradPsi[0][2] - n[2]*gradPsi[0][1] ) ;
			res[1] = factor * ( n[2]*gradPsi[0][0] - n[0]*gradPsi[0][2] ) ;
			res[2] = factor * ( n[0]*gradPsi[0][1] - n[1]*gradPsi[0][0] ) ;
			return res ;
		}
	private:
		double factor;
	};

	struct NablaKNablaPsiFun : public AbstractFunction<double,std::vector<WorldVector<double>*> >
	{
		NablaKNablaPsiFun( double alpha_ = 1.0 ) : AbstractFunction<double,std::vector<WorldVector<double>*> >(6) , alpha(alpha_) {} ;
		double operator()( const std::vector<WorldVector<double>*> &tmpVec ) const 
		{
			double res = alpha ;
			if ( tmpVec.size() >= 2 )
				res = (*tmpVec[0]) * (*tmpVec[1]) ;
			else
				ERROR_EXIT( "WRONG NUMBER OF ITEMS IN std::vector !\n" );
			return res ;
		}
	private:
		double alpha ;
	};
	
	// HELPCLASSES
	struct SqrFun : public AbstractFunction<double,double>
	{
		SqrFun() : AbstractFunction<double,double>(6) {};
		double operator()(const double &x) const { return x * x ; }
	};
	struct MultFun : public AbstractFunction<double,WorldVector<double> >
	{
		MultFun( double a_=1.0 ) : AbstractFunction<double,WorldVector<double> >(6) , a(a_) {};
		double operator()(const WorldVector<double> &x) const { return a * ( x * x ) ; }
	private:
		double a;
	};
	
	template <typename T=double>
	struct ID : public AbstractFunction< T , T >
	{
		ID() : AbstractFunction< T , T >(2) {};
		T operator()(const T& x) const { return x; }
	};

	// Constant function
	struct ConstantFct : public AbstractFunction<double, WorldVector<double> >
	{
		ConstantFct(const double& constant_) : AbstractFunction<double, WorldVector<double> >(1) , constant(constant_) {}
		double operator()(const WorldVector<double>& x) const  { return constant; }
	private:
		double constant;
	};

	struct OneMinusFct : public AbstractFunction<double,double>
	{
		OneMinusFct() : AbstractFunction<double,double>(6) {}
		double operator()(const double& x) const  { return 1.0 - x ; }
	};

	struct NormGradXi : public AbstractFunction<double,WorldVector<double> >
	{
		NormGradXi() : AbstractFunction<double,WorldVector<double> >(6) {}
		double operator()(const WorldVector<double>& gradXi) const  { return std::sqrt( gradXi*gradXi ) ; }
	};
};

#endif
