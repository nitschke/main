#include "AMDiS.h"
#include "MeshHelper.h"

namespace AMDiS {

class MeshCorrector {

  public:
    MeshCorrector(const FiniteElemSpace *finiteElemSpace, std::string name = "");

    void oneIteration(double h);
    
    void oneHeunIteration(double h);

    bool iterate(int n, double h, std::string name = "");

    const FiniteElemSpace *getFeSpace() {return feSpace;}

  protected:

    WorldVector<DOFVector<double> * > coordsDeepCopy();
    

  private:

    FiniteElemSpace *feSpace;
    WorldVector<DOFVector<double> * > coords;
    DOFVector<WorldVector<double> > F;

    int ii;

    MeshInfoCSVWriter infowriter;

};

}
