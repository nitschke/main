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


}
