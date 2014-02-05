#include "AMDiS.h"

namespace AMDiS {

class MeshCorrector {

  public:
    MeshCorrector(const FiniteElemSpace *finiteElemSpace) : feSpace(const_cast<FiniteElemSpace *>(finiteElemSpace)) {
      for (int i = 0; i < 3; i++) {
        newCoords[i] = new DOFVector<double>(feSpace, "paraCoordComponent");
      }
    }

    void oneIteration(double h);

    void iterate(int n, double h);

    const FiniteElemSpace *getFeSpace() {return feSpace;}

  private:

    FiniteElemSpace *feSpace;
    WorldVector<DOFVector<double> * > newCoords;

};

}
