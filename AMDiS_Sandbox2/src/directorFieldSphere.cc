#include <iostream>
#include <fstream>
#include "Dec.h"
#include "ExtremeValueTracker.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

class Michael : public AbstractFunction<WorldVector<double>, WorldVector<double> > {

public:
  Michael(double lambda) : AbstractFunction<WorldVector<double>, WorldVector<double> >(), l(lambda) {}

  WorldVector<double> operator()(const WorldVector<double>& coords) const {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    double cp4 = cos(M_PI / 4.0);
    WorldVector<double> e;
    if (abs(y) >= cp4) {
      e[0] = -x;
      e[1] = 0.0;
      e[2] = -z;
    } else if (x >= cp4) {
      e[0] = 0.0;
      e[1] = y;
      e[2] = z;
    } else if (x <= -cp4) {
      e[0] = 0.0;
      e[1] = sin(M_PI * (y - l));
      e[2] = -sin(M_PI *z);
    } else {
      double c = y / cp4;
      e[0] = abs(c) - 1.0;
      e[1] = c;
      e[2] = 0.0; 
    }
    return e;
  }

private:
  double l;
};


//length of resulting vec depends on the local edge metric
class Noise_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Noise_d(int seed) : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {
    srand(seed);
  }

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return myrand();
  }

private:
  inline double myrand() const{
    return 2.0 * ((double)rand()) / RAND_MAX - 1.0;
  }
};


// <(*d)(z),[p,q]> 
class RotZ_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  RotZ_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    double p1 = p[0];
    double p2 = p[1];
    double p3 = p[2];
    double q1 = q[0];
    double q2 = q[1];
    double q3 = q[2];
    return 2.*(p2*q1 - 1.*p1*q2)*atan((-1. + p1*q1 + p2*q2 + p3*q3)*pow(-2.*p1*p3*q1*q3 - 2.*p2*q2*(p1*q1 + p3*q3) + 
       pow(p3,2.)*(pow(q1,2.) + pow(q2,2.)) + pow(p2,2.)*(pow(q1,2.) + pow(q3,2.)) + pow(p1,2.)*(pow(q2,2.) + pow(q3,2.)),-0.5))*
   pow(-2.*p1*p3*q1*q3 - 2.*p2*q2*(p1*q1 + p3*q3) + pow(p3,2.)*(pow(q1,2.) + pow(q2,2.)) + pow(p2,2.)*(pow(q1,2.) + pow(q3,2.)) + 
     pow(p1,2.)*(pow(q2,2.) + pow(q3,2.)),-0.5);
  }
};

// <d(xyz),[p,q]>
class DXYZ_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DXYZ_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[0]*q[1]*q[2] - p[0]*p[1]*p[2];
  }
};


class ValSquaredMinusOne : public AbstractFunction<double,double> {
public:
  ValSquaredMinusOne() : AbstractFunction<double,double>() {}

  double operator()(const double &c) const {
    return c * c - 1.0;
  }
};


//projection
class Proj : public AbstractFunction<WorldVector<double>, WorldVector<double> > {

  public:
  Proj() : AbstractFunction<WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  WorldVector<double> operator()(const WorldVector<double>& x) const 
  {
    return x / sqrt(x * x);
  }
};

class Delta : public AbstractFunction<double,double> {
public:
  Delta(double eps_) : AbstractFunction<double,double>(), eps(eps_) {}

  double operator()(const double &c) const {
    return exp(- 0.5 * c * c / eps);
  }
private:
  double eps;
};



class MyInstat : public DecProblemInstat {
public:
  MyInstat(DecProblemStat *probStat, const DofEdgeVectorPD &initSol)
      : DecProblemInstat(probStat), 
        normAlpha(initSol.getNormOnEdges()),
        tracker(probStat),
        delta(0.01) {
    string csvfn;
    Parameters::get(probStat->getName() + "->output->filename", csvfn);
    csvfn += "Energies.csv";
    csvout.open(csvfn.c_str(), ios::out);
    csvout << "Time,Div,Rot,Norm,Full" << endl;

    
    double K0 = -1.0;
    Parameters::get("userParameter->K0", K0);
    TEST_EXIT(K0 >= 0.0)("K0 must be positive");
    MinusK0 = -K0;

    double K1 = -1.0;
    Parameters::get("userParameter->K1", K1);
    TEST_EXIT(K1 >= 0.0)("K1 must be positive");
    MinusK1 = -K1;

    double K3 = -1.0;
    Parameters::get("userParameter->K3", K3);
    TEST_EXIT(K3 >= 0.0)("K3 must be positive");
    MinusK3 = -K3;

    Kn = -1.0;
    Parameters::get("userParameter->Kn", Kn);
    TEST_EXIT(Kn >= 0.0)("Kn must be positive");

    
 }


  void closeTimestep() {
    double time = t;
    DecProblemInstat::closeTimestep();
    DofEdgeVectorPD evecPD(statProb->getSolution(0), statProb->getSolution(1));
    normAlpha = evecPD.getNormOnEdges();

    DofEdgeVector scaledNorm(normAlpha);
    scaledNorm.evalFunction(&delta);
    tracker.trackdownMaxima(scaledNorm, time, 0.1);

    DofEdgeVector normDeviat = normAlpha * normAlpha;
    normDeviat += -1.0;

    double ne =  0.25 * Kn * normDeviat.L2NormSquared(); 
    double dive = evecPD.getDirichletEnergy(-0.5 * MinusK1, 0.0);
    double rote = evecPD.getDirichletEnergy(0.0           , -0.5 * MinusK3);
    double energy = dive + rote + ne;
    csvout << time << "," << dive << "," << rote << "," << ne << ","<< energy << endl;

    //if (t > 0.005) {
    //  double eder = (oldEnergy - energy) / oldEnergy;
    //  if (eder < 1.E-4 && eder > 0.0 && tau < 0.1 ) tau *= 2.0;
    //}
    oldEnergy = energy;
  }

  DofEdgeVector* getNormPtr() {
    return &normAlpha;
  }

  double* getMinusK0Ptr() {
    return &MinusK0;
  }

  double* getMinusK1Ptr() {
    return &MinusK1;
  }

  double* getMinusK3Ptr() {
    return &MinusK3;
  }

  double* getKnPtr() {
    return &Kn;
  }

  ~MyInstat() {csvout.close();}

private:
  DofEdgeVector normAlpha;

  double MinusK0;
  double MinusK1;
  double MinusK3;
  double Kn;

  double oldEnergy;

  ofstream csvout;

  Delta delta;
  ExtremeValueTracker tracker;
};



int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  WorldVector<double> ballCenter;
  ballCenter.set(0.0);
  new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  //DofEdgeVector rotz(edgeMesh, "rotz");
  //rotz.set(new RotZ_d());
  //rotz.writeSharpFile("output/rotzSharp.vtu", &sphere);

  //DofEdgeVector dxyz(edgeMesh, "dxyz");
  //dxyz.set(new DXYZ_d());
  //dxyz.writeSharpFile("output/dxyzSharp.vtu", &sphere);

  DofEdgeVectorPD initSol(edgeMesh, "initSol");
  //Noise_d noiseFun(43);
  //initSol.set(&noiseFun, &noiseFun);
  initSol.interpol(new Michael(0.01));
  //initSol.normalize();
  initSol.writeSharpOnEdgesFile("output/initSolSharp.vtu");

  DecProblemStat decSphere(&sphere, edgeMesh);
  MyInstat sphereInstat(&decSphere, initSol);

  EdgeOperator LaplaceB;
  LaplaceB.addTerm(new LaplaceBeltramiAtEdges());
  decSphere.addMatrixOperator(LaplaceB, 0, 0, sphereInstat.getMinusK3Ptr());
  decSphere.addMatrixOperator(LaplaceB, 1, 1, sphereInstat.getMinusK1Ptr());

  EdgeOperator LaplaceCB;
  LaplaceCB.addTerm(new LaplaceCoBeltramiAtEdges());
  decSphere.addMatrixOperator(LaplaceCB, 0, 0, sphereInstat.getMinusK1Ptr());
  decSphere.addMatrixOperator(LaplaceCB, 1, 1, sphereInstat.getMinusK3Ptr());

  EdgeOperator Scale;
  ValSquaredMinusOne vmo;
  Scale.addTerm(new EdgeVecAtEdges(sphereInstat.getNormPtr(), &vmo));
  decSphere.addMatrixOperator(Scale, 0, 0, sphereInstat.getKnPtr());
  decSphere.addMatrixOperator(Scale, 1, 1, sphereInstat.getKnPtr());

  EdgeOperator DtPrimal;
  DtPrimal.addTerm(new IdentityAtEdges());
  DtPrimal.setUhOld((DofEdgeVector)(initSol));
  decSphere.addMatrixOperator(DtPrimal, 0, 0, sphereInstat.getInvTauPtr());
  decSphere.addVectorOperator(DtPrimal, 0, sphereInstat.getInvTauPtr());

  EdgeOperator DtDual;
  DtDual.addTerm(new IdentityAtEdges());
  DtDual.setUhOld(initSol.getDual());
  decSphere.addMatrixOperator(DtDual, 1, 1, sphereInstat.getInvTauPtr());
  decSphere.addVectorOperator(DtDual, 1, sphereInstat.getInvTauPtr());

  //EdgeOperator K;
  //K.addTerm(new IdentityAtEdges(-1.0));
  //decSphere.addMatrixOperator(K, 0, 0);
  //decSphere.addMatrixOperator(K, 1, 1);

  sphereInstat.solve();
}
