#include "AMDiS.h"

namespace AMDiS {

class PhiProject : public Projection {

public:

  PhiProject(int id, 
             ProjectionType type, 
             AbstractFunction<double, WorldVector<double> > *phi_,
             AbstractFunction<WorldVector<double>, WorldVector<double> > *gradPhi_,
             double eps_ = 1.0e-6, int nMax_ = 100) : Projection(id, type), phi(phi_), gradPhi(gradPhi_), eps(eps_), nMax(nMax_) {}

  
  void project(WorldVector<double> &x);

  void project(const WorldVector<double> &x, WorldVector<double> &v);

  WorldVector<double> getNormal(const WorldVector<double> &x);

  ~PhiProject() {
    FUNCNAME("PhiProject::~PhiProject()")
    //MSG("Destroying Projection\n");
    projectionMap[projectionID] = NULL;
  }

private:

  double eps;
  int nMax;
  
  AbstractFunction<double, WorldVector<double> > *phi;
  AbstractFunction<WorldVector<double>, WorldVector<double> > *gradPhi;

};

}
