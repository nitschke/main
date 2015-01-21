#include "AMDiS.h"


using namespace AMDiS;

class QuadraticSpline {

  public:
    QuadraticSpline(const FiniteElemSpace *finiteElemSpace);

    DOFVector<WorldVector<double> > getNormals();

  private: 

    FiniteElemSpace *feSpace;
    DOFVector<mtl::dense_vector<double> > coeffs;
};
