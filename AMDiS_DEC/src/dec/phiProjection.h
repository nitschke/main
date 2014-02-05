#include "AMDiS.h"

namespace AMDiS {

class PhiProject : public Projection {

public:

  PhiProject(int id, 
             ProjectionType type, 
             AbstractFunction<double, WorldVector<double> > *phi_,
             AbstractFunction<WorldVector<double>, WorldVector<double> > *gradPhi_,
             double eps_ = 1.0e-6) : Projection(id, type), phi(phi_), gradPhi(gradPhi_), eps(eps_) {}

  
  void project(WorldVector<double> &x);

private:

  double eps;
  
  AbstractFunction<double, WorldVector<double> > *phi;
  AbstractFunction<WorldVector<double>, WorldVector<double> > *gradPhi;

};

}
