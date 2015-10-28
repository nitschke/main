#include "Dec.h"
#include "io/VtkVectorWriter.h"
#include "MeshMover.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

class GaussCurvNonic :  public AbstractFunction<double, WorldVector<double> > {
  public:

    GaussCurvNonic(double c_, double r_) : AbstractFunction<double, WorldVector<double> >(), c(c_), r(r_) {}

  double operator()(const WorldVector<double>& coords) const {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    return -256*(-16*pow(y,2) - (4*x + c*pow(z,2)*(-(r*(4 + 3*z)*pow(-1 + z,2)) + (-4 + 3*z)*pow(1 + z,2)))*
      (4*x + c*(-(r*(8 + z*(-14 + 3*z*(-24 + z*(-5 + 4*z*(8 + 5*z)))))*pow(-1 + z,2)) + 
           (-8 + z*(-14 + 3*z*(24 + z*(-5 + 4*z*(-8 + 5*z)))))*pow(1 + z,2))))*pow(-1 + pow(z,2),2)*
   (16*pow(y,2) + pow(4*x + c*pow(z,2)*(-(r*(4 + 3*z)*pow(-1 + z,2)) + (-4 + 3*z)*pow(1 + z,2)),2))*
   pow(256*pow(y,4)*pow(z,2) + 32*pow(y,2)*(8 + 
        pow(z,2)*(8*(-2 + 2*pow(x,2) + pow(z,2)) - 
           pow(c,2)*pow(z,2)*(r*(4 + 3*z)*pow(-1 + z,2) - (-4 + 3*z)*pow(1 + z,2))*
            (r*(-8 + z*(-1 + 6*z*(3 + 2*z)))*pow(-1 + z,2) - (8 + z*(-1 + 6*z*(-3 + 2*z)))*pow(1 + z,2)) + 
           4*c*x*(r*(-8 + z*(-1 + z*(14 + 9*z)))*pow(-1 + z,2) - (8 + z*(-1 + z*(-14 + 9*z)))*pow(1 + z,2)))) + 
     (16 + pow(z,2)*(16*(-2 + pow(x,2) + pow(z,2)) + 
           8*c*x*(r*(-8 + z*(-1 + 6*z*(3 + 2*z)))*pow(-1 + z,2) - (8 + z*(-1 + 6*z*(-3 + 2*z)))*pow(1 + z,2)) + 
           pow(c,2)*pow(-(r*(-8 + z*(-1 + 6*z*(3 + 2*z)))*pow(-1 + z,2)) + 
              (8 + z*(-1 + 6*z*(-3 + 2*z)))*pow(1 + z,2),2)))*
      pow(4*x + c*pow(z,2)*(-(r*(4 + 3*z)*pow(-1 + z,2)) + (-4 + 3*z)*pow(1 + z,2)),2),-2);
  }

  private:
    
    double c;
    double r;
};

class GaussCurvSphere :  public AbstractFunction<double, WorldVector<double> > {
  public:

    GaussCurvSphere() : AbstractFunction<double, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& coords) const {
    return 1.0;
  }

 };


int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  double c = -1.0;
  Parameters::get("nonicParameter->c", c);
  TEST_EXIT(c >= 0.0)("stretch factor c must be positive");

  double r = -1.0;
  Parameters::get("nonicParameter->r", r);
  TEST_EXIT(r >= 0.0)("stretch ratio factor r must be positive");

  string bn;
  Parameters::get("sphere->output->filename", bn);

  DOFVector<double> K(sphere.getFeSpace(), "GaussCurv");
  //K.interpol(new GaussCurvNonic(c, r));
  K.interpol(new GaussCurvSphere());
  io::VtkVectorWriter::writeFile(K, bn + "GausCurv.vtu");

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());
  DecProblemStat decSphere(&sphere, edgeMesh);

  VertexOperator LB;
  LB.addTerm(new LaplaceBeltramiAtVertices());
  decSphere.addMatrixOperator(LB, 0, 0);

  VertexOperator RHS;
  RHS.addTerm(new VertexVecAtVertices(&K));
  decSphere.addVectorOperator(RHS, 0);

  decSphere.assembleSystem();
  decSphere.solve();
  decSphere.writeSolution();
}
