#include "AMDiS.h"

using namespace AMDiS;
using namespace mtl;

class Quadratic {

  public:
    Quadratic(const FiniteElemSpace *finiteElemSpace);

  private: 

    FiniteElemSpace *feSpace;
    dense_vector<double> coeffs;

};
