#include "AMDiS.h"

namespace AMDiS {

  using namespace std;

  class DecOperator : public Operator {

    public:
      DecOperator(const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : Operator(rowFeSpace, colFeSpace) {
        opFactor = 1.0;
       }

      void setFactor(double factor) {
        opFactor = factor;
      }

    protected:
      double opFactor;
  
  };
  
  class LBeltramiDEC : public DecOperator {
    
    public:
      LBeltramiDEC(const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : DecOperator(rowFeSpace, colFeSpace) {}

      void getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor = 1.0);
    
  };

  class FunctionDEC : public DecOperator {
    
    public:
      
      FunctionDEC(const FiniteElemSpace *rowFeSpace, AbstractFunction<double, WorldVector<double> > *fun,
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

      void getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor = 1.0);
    
  };


}
