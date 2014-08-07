#include "AMDiS.h"
#include "decOperator.h"
#include "MatrixHelper.h"
#include "MeshHelper.h"
#include "WorldVectorHelper.h"
#include "DOFVHelper.h"
#include "torusProjection.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

class MC : public AbstractFunction<double, WorldVector<double> >
{
public:
  MC() : AbstractFunction<double, WorldVector<double> >(1) {}

  /// Implementation of AbstractFunction::operator().
  //double operator()(const WorldVector<double>& x) const 
  double operator()(const WorldVector<double>& coords) const 
  {
    double r = 0.5;
    double R = 2.0;
    //double zeta =  (abs(x[1]) < r) ? (sqrt(r*r - x[1]*x[1])) : 0.0;
    //if (x[0]*x[0] + x[2]*x[2] < R*R) zeta *= -1.0;
    //double eta = zeta / (R + zeta);
    //return 0.5 * (1.0 + eta) / r;

    double x = coords[0];
    double y = coords[2]; //coords change 2 <-> 1
    double z = coords[1];    
    double sx2py2 =  sqrt(x*x + y*y);
    return 0.5 * (2.0 - R / sx2py2) / r;
  }
};

class GC : public AbstractFunction<double, WorldVector<double> >
{
public:
  GC() : AbstractFunction<double, WorldVector<double> >(1) {}

  /// Implementation of AbstractFunction::operator().
  //double operator()(const WorldVector<double>& x) const 
  double operator()(const WorldVector<double>& coords) const 
  {
    double r = 0.5;
    double R = 2.0;
    //double zeta =  (true) ? (sqrt(r*r - x[1]*x[1])) : 0.0;
    //if (x[0]*x[0] + x[2]*x[2] < R*R) zeta *= -1.0;
    //double eta = zeta / (R + zeta);
    //return eta / (r * r);
    double x = coords[0];
    double y = coords[2]; //coords change 2 <-> 1
    double z = coords[1];    
    double sx2py2 =  sqrt(x*x + y*y);
    return (1.0 - R / sx2py2) / (r*r);
  }
};

class Normal : public AbstractFunction<double, WorldVector<double> >
{
public:
  Normal(int i_) : AbstractFunction<double, WorldVector<double> >(1), i(i_) {}

  double operator()(const WorldVector<double>& coords) const 
  {
    double r = 0.5;
    double R = 2.0;
    double x = coords[0];
    double y = coords[2]; //coords change 2 <-> 1
    double z = coords[1];    
    double rval;
    if (i == 1) {
      rval = z;
    }
    else {
      double sx2py2 =  sqrt(x*x + y*y);
      rval = (sx2py2 - R) / sx2py2;
      rval *= (i==0) ? x : y;
    }
    return rval / r;
  }

private:
  int i;
};

class NormalVec : public AbstractFunction< WorldVector<double>, WorldVector<double> >
{
public:
  NormalVec() : AbstractFunction< WorldVector<double>, WorldVector<double> >(0), n0(0), n1(1), n2(2) {}

  WorldVector<double> operator()(const WorldVector<double>& x) const
  {
    WorldVector<double> rval;
    rval[0] = n0(x);
    rval[1] = n1(x);
    rval[2] = n2(x);
    return rval;
  }

private:
  Normal n0;
  Normal n1;
  Normal n2;
};


// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("torus main");

  AMDiS::init(argc, argv);

  new TorusProject(1, VOLUME_PROJECTION, 2.0, 0.5);

  // ===== create and init the scalar problem ===== 
  ProblemStat torus("torus");
  torus.initialize(INIT_ALL);


  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("torus->adapt", torus.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("torus->adapt",
					       &torus,
					       adaptInfo);
  
  //DOFVector<WorldVector<double> > avNormalDOFV = getNormals(torus.getFeSpace());
  DOFVector<WorldVector<double> > avNormalDOFV = getNormalsConnectionAverage(torus.getFeSpace());

  int oh = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
       int pos = matIndex(i,j);
       SimpleDEC *II = new SimpleDEC(torus.getFeSpace(pos), torus.getFeSpace(pos));
       torus.addMatrixOperator(II, pos, pos);

       GradDofWorldVecDEC *dN = new GradDofWorldVecDEC(i, j, &avNormalDOFV, torus.getFeSpace(pos), torus.getFeSpace(pos));
       torus.addVectorOperator(dN, pos);
    }
  }
  

  // ===== start adaption loop =====
  adapt->adapt();

  WorldMatrix<DOFVector<double> * > IIDV;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      IIDV[i][j] = torus.getSolution(matIndex(i,j) + oh);
    }
  }

  DOFVector<WorldVector<double> > eigDofVector = getEigenVals(IIDV);
  VtkVectorWriter::writeFile(eigDofVector, string("output/eigenVals.vtu"));

  DOFVector<double> gcDOFV(torus.getFeSpace(),"GaussCurvExact");
  gcDOFV.interpol(new GC());
  VtkVectorWriter::writeFile(gcDOFV, string("output/gaussExact.vtu"));

  DOFVector<double> mcDOFV(torus.getFeSpace(),"MeanCurvExact");
  mcDOFV.interpol(new MC());
  VtkVectorWriter::writeFile(mcDOFV, string("output/meanExact.vtu"));

  DOFVector<double> gcWeingarten = prod01(eigDofVector);
  VtkVectorWriter::writeFile(gcWeingarten, string("output/gaussWeingarten.vtu"));
  printError(gcWeingarten, gcDOFV, "GaussWeingarten");

  DOFVector<double> mcWeingarten = halfSum01(eigDofVector);
  VtkVectorWriter::writeFile(mcWeingarten, string("output/meanWeingarten.vtu"));
  printError(mcWeingarten, mcDOFV, "MeanWeingarten");


  MeshInfoCSVWriter mwriter("/dev/null/nonaynever.csv");
  mwriter.appendData(torus.getFeSpace(),true);

  torus.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


