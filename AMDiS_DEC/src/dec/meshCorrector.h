#include "AMDiS.h"

namespace AMDiS {

class MeshCorrector {

  public:
    MeshCorrector(const FiniteElemSpace *finiteElemSpace);

    void oneIteration(double h);
    
    void oneHeunIteration(double h);

    void iterate(int n, double h);

    const FiniteElemSpace *getFeSpace() {return feSpace;}

  protected:

    WorldVector<DOFVector<double> * > coordsDeepCopy();
    

  private:

    FiniteElemSpace *feSpace;
    WorldVector<DOFVector<double> * > coords;
    DOFVector<WorldVector<double> > F;

};

}
