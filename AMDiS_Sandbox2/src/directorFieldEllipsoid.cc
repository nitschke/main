#include <iostream>
#include <fstream>
#include "Dec.h"
#include "ExtremeValueTracker.h"
#include "phiProjection.h"


using namespace std;
using namespace AMDiS;
using namespace dec;



//length of resulting vec depends on the local edge metric
class Noise_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Noise_d(int seed, double fac = 1.0) : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {
    srand(seed); f = fac;
  }

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    //return f;
    return myrand();
  }

private:
  inline double myrand() const{
    return 2.0 * ((double)rand()) / RAND_MAX - 1.0;
  }

  double f;
};

class NoiseFac : public AbstractFunction<double,double> {
public:
  NoiseFac() : AbstractFunction<double,double>() {}

  double operator()(const double &c) const {
    return (2.0 * ((double)rand()) / RAND_MAX - 1.0) * c;
  }
};

// <alpha,[p,q]>
class one_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  one_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return 1.0;
  }
};



class ValSquaredMinusOne : public AbstractFunction<double,double> {
public:
  ValSquaredMinusOne() : AbstractFunction<double,double>() {}

  double operator()(const double &c) const {
    return c * c - 1.0;
  }
};

class ValSquared : public AbstractFunction<double,double> {
public:
  ValSquared() : AbstractFunction<double,double>() {}

  double operator()(const double &c) const {
    return c * c;
  }
};

class Prod2 : public TertiaryAbstractFunction<double,double,double, EdgeElement> {
public:
  Prod2() : TertiaryAbstractFunction<double,double,double, EdgeElement>() {}

  double operator()(const double &a1, const double &a2, const EdgeElement &eel) const {
    double lenP = eel.infoLeft->getEdgeLen(eel.dofEdge);
    return 2.0 * a1 * a2 / (lenP * lenP);
  }
};


class AlphaGrad1 : public TertiaryAbstractFunction<double,double,double, EdgeElement> {
public:
  AlphaGrad1() : TertiaryAbstractFunction<double,double,double, EdgeElement>() {}

  double operator()(const double &a1, const double &a2, const EdgeElement &eel) const {
    double lenP = eel.infoLeft->getEdgeLen(eel.dofEdge);
    return (3.0 * a1 * a1 + a2 * a2) / (lenP * lenP) - 1.0;
  }
};

class AlphaGrad2 : public TertiaryAbstractFunction<double,double,double, EdgeElement> {
public:
  AlphaGrad2() : TertiaryAbstractFunction<double,double,double, EdgeElement>() {}

  double operator()(const double &a1, const double &a2, const EdgeElement &eel) const {
    double lenP = eel.infoLeft->getEdgeLen(eel.dofEdge);
    return 2.0 * a1 * a2 / (lenP * lenP);
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

// for phi-projection
class Phi : public AbstractFunction<double, WorldVector<double> >
{
public:
  Phi() : AbstractFunction<double, WorldVector<double> >(1) {}

  double operator()(const WorldVector<double>& x) const 
  {
    return 0.5 * (x[0]*x[0] + 4.0*x[1]*x[1] + (4.0/9.0)*x[2]*x[2] - 1.0);
  }
};

class GradPhi : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  GradPhi() : AbstractFunction<WorldVector<double>, WorldVector<double> >(1) {}

  WorldVector<double> operator()(const WorldVector<double>& x) const 
  {
    WorldVector<double> rval(x);
    rval[1] *= 4.0;
    rval[2] *= 4.0/9.0;
    return rval;
  }
};


class MyInstat : public DecProblemInstat {
public:
  MyInstat(DecProblemStat *probStat, DofEdgeVectorPD initSol)
      : DecProblemInstat(probStat), 
        normAlpha(initSol.getNormOnEdges()),
        solPrimal((DofEdgeVector)(initSol)),
        solDual(initSol.getDual()),
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
    solPrimal = (DofEdgeVector)(evecPD);
    solDual = evecPD.getDual();

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

  DofEdgeVector* getSolPrimal() {
    return &solPrimal;
  }

  DofEdgeVector* getSolDual() {
    return &solDual;
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

  double* getMinusKnPtr() {
    return &MinusKn;
  }

  ~MyInstat() {csvout.close();}

private:
  DofEdgeVector normAlpha;
  DofEdgeVector solPrimal;
  DofEdgeVector solDual;

  double MinusK0;
  double MinusK1;
  double MinusK3;
  double Kn;
  double MinusKn;

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
  new PhiProject(1, VOLUME_PROJECTION, new Phi(), new GradPhi(), 1.0e-8);
  

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
  Noise_d noiseFun(43);
  //initSol.set(&noiseFun, new Noise_d(43,-1./3.));
  initSol.set(&noiseFun, &noiseFun);
  initSol.normalize(1.E-10);
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

  // explicite -> need little timesteps
  EdgeOperator Scale;
  ValSquaredMinusOne vmo;
  Scale.addTerm(new EdgeVecAtEdges(sphereInstat.getNormPtr(), &vmo));
  decSphere.addMatrixOperator(Scale, 0, 0, sphereInstat.getKnPtr());
  decSphere.addMatrixOperator(Scale, 1, 1, sphereInstat.getKnPtr());

  DofEdgeVector initPrimal = (DofEdgeVector)(initSol);
  DofEdgeVector initDual =  initSol.getDual();



  // time derivation
  EdgeOperator DtPrimal;
  DtPrimal.addTerm(new IdentityAtEdges());
  DtPrimal.setUhOld(initPrimal, 0);
  decSphere.addMatrixOperator(DtPrimal, 0, 0, sphereInstat.getInvTauPtr());
  decSphere.addVectorOperator(DtPrimal, 0, sphereInstat.getInvTauPtr());

  EdgeOperator DtDual;
  DtDual.addTerm(new IdentityAtEdges());
  DtDual.setUhOld(initDual, 1);
  decSphere.addMatrixOperator(DtDual, 1, 1, sphereInstat.getInvTauPtr());
  decSphere.addVectorOperator(DtDual, 1, sphereInstat.getInvTauPtr());

  //shake
  DofEdgeVector one(edgeMesh, "one");
  one.set(new one_d());
  EdgeOperator shake;
  shake.addTerm(new EdgeVecAtEdges(&one, new NoiseFac(), 1.0e-0));
  decSphere.addVectorOperator(shake, 0);
  decSphere.addVectorOperator(shake, 1);


  //decSphere.assembleSystem();
  //using namespace mtl;
  //dense2D<double> mat(60,60);
  //mat = decSphere.getSysMat();
  //edgeMesh->printVolInfos(true,true);
  //cout << sub_matrix(mat,0,30,0,30) << endl << endl;
  //cout << sub_matrix(mat,0,30,30,60) << endl << endl;
  //cout << sub_matrix(mat,30,60,0,30) << endl << endl;
  //cout << sub_matrix(mat,30,60,30,60) << endl << endl;
  //cout << decSphere.getRhs() << endl;

  sphereInstat.solve();
}
