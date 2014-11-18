#include "AMDiS.h"

namespace AMDiS {

class DummyProject : public Projection {

public:

  DummyProject(int id, 
             ProjectionType type 
             ) : Projection(id, type){}

  
  void project(WorldVector<double> &x);

  void project(const WorldVector<double> &x, WorldVector<double> &v);



private:

};

}
