#include "AMDiS.h"
#include "decOperator.h"
#include "phiProjection.h"
#include "meshCorrector.h"
#include "MeshHelper.h"
#include "WorldVectorHelper.h"
#include "MatrixHelper.h"
#include "DOFVHelper.h"
#include "quadraticSpline.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

class X : public AbstractFunction<double, WorldVector<double> >
{
public:
  X(int i_) : AbstractFunction<double, WorldVector<double> >(1), i(i_) {}

  double operator()(const WorldVector<double>& x) const 
  {
    return x[i];
  }

protected:
  int i;
};

class Phi : public AbstractFunction<double, WorldVector<double> >
{
public:
  Phi() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& x) const 
  {
    double z2 = x[2] * x[2];
    double xMz2 = x[0] - z2;
    double yMz2 = x[1] - z2;
    return 0.5 * (xMz2 * xMz2 + yMz2 * yMz2 + z2 - 1.0);
  }
};

class GradPhi : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  GradPhi() : AbstractFunction<WorldVector<double>, WorldVector<double> >(1) {}

  WorldVector<double> operator()(const WorldVector<double>& x) const 
  {
    WorldVector<double> rval(x);
    double z2 = x[2] * x[2];
    rval[0] -= z2;
    rval[1] -= z2;
    rval[2] += 2.0 * x[2]* (2.0 * z2 - x[0] - x[1]);
    return rval;
  }
};

class Normal : public AbstractFunction<double, WorldVector<double> >
{
public:
  Normal(int i_) : AbstractFunction<double, WorldVector<double> >(1), dPhi(), i(i_) {}

  double operator()(const WorldVector<double>& x) const 
  {
    WorldVector<double> rval = dPhi(x);
    return ((1.0/sqrt(dot(rval,rval))) * rval)[i];
  }

private:
  GradPhi dPhi;
  int i;
  
};

class GC : public AbstractFunction<double, WorldVector<double> >
{
public:
  GC() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& coord) const 
  {
    double x = coord[0];
    double y = coord[1];
    double z = coord[2];
    //double c = 81.0 + 972.0*y*y - 20.0*z*z;
    //return 11664.0 / (c * c);
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double oben = 1.0 + 2.0*x - 2.0*x2 + 2.0*y - 2.0*y2 + 2.0*(-3.0 + x + y)*z2;
    double sqUnten = 1.0 - 4.0*(-2.0 + x + x2 + y - 2.0*x*y + y2)*z2;
    return  - oben / (sqUnten*sqUnten);
  }
};

class MC : public AbstractFunction<double, WorldVector<double> >
{
public:
  MC() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& coord) const 
  {
    double x = coord[0];
    double y = coord[1];
    double z = coord[2];
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double oben = x - x2 + y - y2 + (-7.0 + 3.0*x + 3.0*y)*z2;
    double unten23 = 1.0 - 4.0*(-2.0 + x + x2 + y - 2.0*x*y + y2)*z2;
    return - oben / sqrt(unten23*unten23*unten23);
  }
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
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  new PhiProject(1, VOLUME_PROJECTION, new Phi(), new GradPhi(), 1.0e-8);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);


  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
					       &sphere,
					       adaptInfo);
  
  //DOFVector<WorldVector<double> > vertexNormals = getNormalsAngleEdgeReciprocalAverage(sphere.getFeSpace());  

  QuadraticSpline quad(sphere.getFeSpace());
  DOFVector<WorldVector<double> > vertexNormals = quad.getNormals();

  int oh = -3;
  //int oh = 0;
  for (int i = 0; i < 3; i++) {
    // N
    //sphere.addMatrixOperator(new SimpleDEC(sphere.getFeSpace(i+oh), sphere.getFeSpace(i+oh)), i+oh, i+oh);
    //sphere.addVectorOperator(new DualPrimalNormalDEC(i, sphere.getFeSpace(i+oh)), i+oh);
    for (int j = 0; j < 3; j++) {
       int pos = matIndex(i,j) + oh + 3;
       //int pos = matIndex(i,j) + oh;
       // -II_ij
       SimpleDEC *II = new SimpleDEC(sphere.getFeSpace(pos), sphere.getFeSpace(pos));
       II->setFactor(-1.0);
       sphere.addMatrixOperator(II, pos, pos);
       // [d(N_j)]_i
       //PrimalPrimalGradDEC *dN = new PrimalPrimalGradDEC(i, sphere.getFeSpace(pos), sphere.getFeSpace(j+oh));
       //sphere.addMatrixOperator(dN, pos, j+oh);
       //PrimalPrimalGradFunctionDEC *dN = new PrimalPrimalGradFunctionDEC(i, new Normal(j), sphere.getFeSpace(pos), sphere.getFeSpace(pos));
       GradDofWorldVecDEC *dN = new GradDofWorldVecDEC(i, j, &vertexNormals, sphere.getFeSpace(pos), sphere.getFeSpace(pos));       
       dN->setFactor(-1.0);
       sphere.addVectorOperator(dN, pos);
    }
  }
  

  // ===== start adaption loop =====
  adapt->adapt();
  
  WorldMatrix<DOFVector<double> * > IIDV;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      IIDV[i][j] = sphere.getSolution(matIndex(i,j) + oh + 3);
      //IIDV[i][j] = sphere.getSolution(matIndex(i,j) + oh);
    }
  }

  DOFVector<WorldVector<double> > eigDofVector = getEigenVals(IIDV);
  VtkVectorWriter::writeFile(eigDofVector, string("output/eigenVals.vtu"));


  DOFVector<double> gcDOFV(sphere.getFeSpace(),"GaussCurvExact");
  gcDOFV.interpol(new GC());
  VtkVectorWriter::writeFile(gcDOFV, string("output/gaussExact.vtu"));

  DOFVector<double> mcDOFV(sphere.getFeSpace(),"MeanCurvExact");
  mcDOFV.interpol(new MC());
  VtkVectorWriter::writeFile(mcDOFV, string("output/meanExact.vtu"));

  DOFVector<double> gcWeingarten = prod01(eigDofVector);
  printError(gcWeingarten, gcDOFV, "GaussWeingarten");
  VtkVectorWriter::writeFile(gcWeingarten, "output/GaussWeingarten.vtu");

  DOFVector<double> mcWeingarten = halfSum01(eigDofVector);
  printError(mcWeingarten, mcDOFV, "MeanWeingarten");
  VtkVectorWriter::writeFile(mcWeingarten, "output/MeanWeingarten.vtu");

  DOFVector<WorldVector<double> > normalDOFV(sphere.getFeSpace(),"Normal");
  normalDOFV.interpol(new NormalVec());
  VtkVectorWriter::writeFile(normalDOFV, string("output/exactNormals.vtu"));
  VtkVectorWriter::writeFile(vertexNormals, string("output/avNormals.vtu"));
  printError(normalDOFV, vertexNormals , "Normal");

  //MeshInfoCSVWriter mwriter("/dev/null/nonaynever.csv");
  //mwriter.appendData(sphere.getFeSpace(),true);

  //quad.writeVTK(string("output/approxSurvace.vtu"));
  
  
  //sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}



