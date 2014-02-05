#include "AMDiS.h"

namespace AMDiS {

class MeshCorrector {

  public:
    MeshCorrector(const FiniteElemSpace *finiteElemSpace) : feSpace(const_cast<FiniteElemSpace *>(finiteElemSpace)) {}

    void oneIteration(double h);

    const FiniteElemSpace *getFeSpace() {return feSpace;}

  private:

    FiniteElemSpace *feSpace;

};

}
