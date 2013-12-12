#include "AMDiS.h"

namespace AMDiS {

  using namespace std;

  class SandboxOperator : public Operator {
    
    public:

      SandboxOperator(const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : Operator(rowFeSpace, colFeSpace) {}

      void getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor) {
        //if (!assembler) initAssembler(NULL, NULL, NULL, NULL);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            userMat(i,j) = factor;
          }
        }

         cout << userMat << endl;
      }
  };

  
  class LBeltramiDEC : public Operator {
    
    public:
      LBeltramiDEC(const FiniteElemSpace *rowFeSpace,
	     const FiniteElemSpace *colFeSpace = NULL) : Operator(rowFeSpace, colFeSpace) {}

      void getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor = 1.0);
    
  };

  class FunctionDEC : public Operator {
    
    public:
      
      FunctionDEC(const FiniteElemSpace *rowFeSpace, AbstractFunction<double, WorldVector<double> > *fun,
	     const FiniteElemSpace *colFeSpace = NULL) : Operator(rowFeSpace, colFeSpace), f(fun) {}

      void getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor = 1.0);
    
    protected:
      
      AbstractFunction<double, WorldVector<double> > *f;

  };

}
