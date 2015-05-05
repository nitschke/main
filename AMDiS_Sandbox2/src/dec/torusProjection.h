#include "AMDiS.h"

namespace AMDiS {

class TorusProject : public Projection
{
public:
  /// Constructor.
  TorusProject(int id, 
	ProjectionType type,
	double R_= 2.0,
  double r_= 0.5) 
    : Projection(id, type),
  R(R_),
  r(r_)
  {}

  /// Destructor.
  virtual ~TorusProject() {}

  /// Implementation of Projection::project();
  void project(WorldVector<double> &x);

protected:
  double R;
  double r;
};

}
