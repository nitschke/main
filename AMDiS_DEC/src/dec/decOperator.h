#include "AMDiS.h"

namespace AMDiS {

  using namespace std;

  class DecOperator : public Operator {

    public:
      DecOperator(const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : Operator(rowFeSpace, colFeSpace), opMat(3,3), opVec(3) {
        opFactor = 1.0;
       }

      void setFactor(double factor) {
        opFactor = factor;
      }

      void getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor = 1.0) {
          FUNCNAME("DecOperator::getElementMatrix");
          ERROR_EXIT("Sorry, this operation is not implemented...that makes me sad :( ");
        }

      void getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor = 1.0) {
            FUNCNAME("DecOperator::getElementVector");
            ERROR_EXIT("Sorry, this operation is not implemented...that makes me sad :( ");
        }


    protected:
      inline void updateUserMat(ElementMatrix& userMat, ElementMatrix& opMat);
      inline void updateUserVec(ElementVector& userVec, ElementVector& opVec);

    protected:
      double opFactor;

      ElementMatrix opMat;
      ElementVector opVec;
  };
  
  class LBeltramiDEC : public DecOperator {
    
    public:
      LBeltramiDEC(const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : DecOperator(rowFeSpace, colFeSpace) {}

      void getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor = 1.0);
    
  };

  
  // Laplace(q*f)
  class LBeltramiInteriorFunctionDEC : public DecOperator {
    
    public:
      LBeltramiInteriorFunctionDEC(AbstractFunction<double, WorldVector<double> > *qFun, const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : DecOperator(rowFeSpace, colFeSpace), q(qFun) {}

      void getElementVector(const ElInfo *elInfo, 
		    ElementVector& userVec, 
				double factor = 1.0);
    
     protected:
      
      AbstractFunction<double, WorldVector<double> > *q;

  };


  
  // rot(k*rot(f)) ~ "div(k*grad(f))"
  class LBeltramiWeightedDEC : public DecOperator {
    
    public:
      LBeltramiWeightedDEC(AbstractFunction<double, WorldVector<double> > *kFun, const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : DecOperator(rowFeSpace, colFeSpace), kappa(kFun) {}

      void getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor = 1.0);
    
    protected:
      
      AbstractFunction<double, WorldVector<double> > *kappa;

  };

  class PrimalPrimalGradFunctionDEC : public DecOperator {
    
    public:
      PrimalPrimalGradFunctionDEC(int direction, AbstractFunction<double, WorldVector<double> > *fun, const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : DecOperator(rowFeSpace, colFeSpace), f(fun), l(direction) {}

      void getElementVector(const ElInfo *elInfo, 
		    ElementVector& userVec, 
				double factor = 1.0);
    
     protected:
      
      AbstractFunction<double, WorldVector<double> > *f;
      int l;

  };

  class PrimalPrimalGradDEC : public DecOperator {
    
    public:
      PrimalPrimalGradDEC(int direction, const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : DecOperator(rowFeSpace, colFeSpace), l(direction) {}

      void getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor = 1.0);
    
     protected:
      
      int l;

  };

  class DualPrimalNormalDEC : public DecOperator {
    
    public:
      DualPrimalNormalDEC(int direction, const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : DecOperator(rowFeSpace, colFeSpace), l(direction) {}

      void getElementVector(const ElInfo *elInfo, 
		    ElementVector& userVec, 
				double factor = 1.0);
    
     protected:
      
      int l;

  };


  class JacobianDEC : public DecOperator {
    
    public:
      
      JacobianDEC(AbstractFunction<double, WorldVector<double> > *phiFun, const FiniteElemSpace *rowFeSpace, 
	     const FiniteElemSpace *colFeSpace = NULL) : DecOperator(rowFeSpace, colFeSpace), phi(phiFun) {}

      //void getElementVector(const ElInfo *elInfo, 
			//	  ElementVector& userVec, 
			//	  double factor = 1.0);

      void getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor = 1.0);
    
    protected:
      
      AbstractFunction<double, WorldVector<double> > *phi;

  };

  class FunctionDEC : public DecOperator {
    
    public:
      
      FunctionDEC(AbstractFunction<double, WorldVector<double> > *fun, const FiniteElemSpace *rowFeSpace, 
	     const FiniteElemSpace *colFeSpace = NULL) : DecOperator(rowFeSpace, colFeSpace), f(fun) {}

      void getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor = 1.0);

      void getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor = 1.0);
    
    protected:
      
      AbstractFunction<double, WorldVector<double> > *f;

  };

  class SimpleDEC : public DecOperator {
    
    public:
      SimpleDEC(const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : DecOperator(rowFeSpace, colFeSpace) {}

      void getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor = 1.0);

      void getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor = 1.0);
    
  };


  class GaussCurvatureDEC : public DecOperator {
    
    public:
      GaussCurvatureDEC(const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : DecOperator(rowFeSpace, colFeSpace) {}

      void getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor = 1.0);

  };

}
