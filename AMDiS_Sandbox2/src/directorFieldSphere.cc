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


// Twist linear form angle 0 to pi/2
class Twist : public AbstractFunction<WorldVector<double>, WorldVector<double> > {

public:
  Twist() : AbstractFunction<WorldVector<double>, WorldVector<double> >() {}

  WorldVector<double> operator()(const WorldVector<double>& coords) const {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldVector<double> r;
    double eps = 1.E-4;
    if (abs(z) < 1.0 - eps) {
      double t1 = z + 1.0;
      double t2 = sqrt(2.0/t1 - 1.0);
      double t3 = sqrt(2.0 - 2.0*z);
      r[0] = (x*z - y*t2) / t3;
      r[1] = (y*z + x*t2) / t3;
      r[2] = -sqrt(1.0 - z)*t1/sqrt(2.0);
    } else {
      r = 0;
    }
    return r;
  }
};

// twist pi/4
class TwistQP : public AbstractFunction<WorldVector<double>, WorldVector<double> > {

public:
  TwistQP() : AbstractFunction<WorldVector<double>, WorldVector<double> >() {}

  WorldVector<double> operator()(const WorldVector<double>& coords) const {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldVector<double> r;
    double eps = 1.E-6;
    if (abs(z) < 1.0 - eps) {
      double t1 = 1.0 / sqrt(2.0 - 2.0*z*z);
      r[0] = (x*z - y) * t1;
      r[1] = (y*z + x) * t1;
      r[2] = (z*z - 1.0) * t1;
    } else {
      r = 0;
    }
    return r;
  }
};



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

// <d(z^2),[p,q]>
class DZ2_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DZ2_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[2]*q[2]*q[0]*q[1] - p[2]*p[2]*p[0]*p[1];
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

    animWriter = new AnimationWriter("output/sphereSharpOnEdges.pvd");
 }


  void closeTimestep() {
    double time = t;
    DecProblemInstat::closeTimestep();
    solPrimal = statProb->getSolution(0);
    solDual =  statProb->getSolution(1);
    DofEdgeVectorPD evecPD(solPrimal, solDual);

    normAlpha = evecPD.getNormOnEdges();

    DofEdgeVector scaledNorm(normAlpha);
    scaledNorm.evalFunction(&delta);
    tracker.trackdownMaxima(scaledNorm, time, 0.1);

    DofEdgeVector normDeviat = normAlpha * normAlpha;
    normDeviat += -1.0;

    double ne =  0.25 * Kn * normDeviat.L2NormSquared(); 
    double dive = evecPD.getDirichletEnergy(-0.5 * MinusK3, 0.0);
    double rote = evecPD.getDirichletEnergy(0.0           , -0.5 * MinusK1);
    double energy = dive + rote + ne;
    csvout << time << "," << dive << "," << rote << "," << ne << ","<< energy << endl;

    //if (t > 0.005) {
    //  double eder = (oldEnergy - energy) / oldEnergy;
    //  if (eder < 1.E-4 && eder > 0.0 && tau < 0.1 ) tau *= 2.0;
    //}
    oldEnergy = energy;

    int prec = 3;
    ostringstream timeoss;
    timeoss << setprecision(prec) << time;
    string fn = "output/sphereSharpOnEdges." + timeoss.str() + ".vtu";
    evecPD.writeSharpOnEdgesFile(fn);
    animWriter->updateAnimationFile(time,fn);
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

  AnimationWriter *animWriter;
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
  Noise_d noiseFun(43);
  //initSol.set(&noiseFun, new Noise_d(43,-1./3.));
  //initSol.set(&noiseFun);
  //initSol.set(new DZ2_d());
  //initSol.bakeDual();
  //initSol.interpol(new Michael(0.01));
  initSol.interpol(new TwistQP());
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

  //// explicite -> need little timesteps
  //EdgeOperator Scale;
  //ValSquaredMinusOne vmo;
  //Scale.addTerm(new EdgeVecAtEdges(sphereInstat.getNormPtr(), &vmo));
  //decSphere.addMatrixOperator(Scale, 0, 0, sphereInstat.getKnPtr());
  //decSphere.addMatrixOperator(Scale, 1, 1, sphereInstat.getKnPtr());

  DofEdgeVector initPrimal = (DofEdgeVector)(initSol);
  DofEdgeVector initDual =  initSol.getDual();


  //Taylor linisarisation
  //(||alpha||^2-1)*alpha
  ValSquaredMinusOne vmo;

  EdgeOperator Norm2Primal;
  Norm2Primal.addTerm(new EdgeVecAtEdges(sphereInstat.getNormPtr(), new ValSquared()));
  Norm2Primal.setUhOld(initPrimal, 0);
  decSphere.addMatrixOperator(Norm2Primal, 0, 0, sphereInstat.getKnPtr());
  decSphere.addVectorOperator(Norm2Primal, 0, sphereInstat.getKnPtr());
  decSphere.addVectorOperator(Norm2Primal, 0, sphereInstat.getKnPtr()); // 2 * ||alpha||^2 *alpha on RHS
  
  EdgeOperator Norm2Dual;
  Norm2Dual.addTerm(new EdgeVecAtEdges(sphereInstat.getNormPtr(), new ValSquared()));
  Norm2Dual.setUhOld(initDual, 1);
  decSphere.addMatrixOperator(Norm2Dual, 1, 1, sphereInstat.getKnPtr());
  decSphere.addVectorOperator(Norm2Dual, 1, sphereInstat.getKnPtr());
  decSphere.addVectorOperator(Norm2Dual, 1, sphereInstat.getKnPtr()); // 2 * ||alpha||^2 *alpha on RHS

  EdgeOperator Id;
  Id.addTerm(new IdentityAtEdges(-1.0));
  decSphere.addMatrixOperator(Id, 0, 0, sphereInstat.getKnPtr());
  decSphere.addMatrixOperator(Id, 1, 1, sphereInstat.getKnPtr());

  EdgeOperator PrimalPrimal;
  PrimalPrimal.addTerm(new EdgeVec2AndEdgeAtEdges(sphereInstat.getSolPrimal(), sphereInstat.getSolPrimal(), new Prod2()));
  decSphere.addMatrixOperator(PrimalPrimal, 0, 0, sphereInstat.getKnPtr());
  
  EdgeOperator DualDual;
  DualDual.addTerm(new EdgeVec2AndEdgeAtEdges(sphereInstat.getSolDual(), sphereInstat.getSolDual(), new Prod2()));
  decSphere.addMatrixOperator(DualDual, 1, 1, sphereInstat.getKnPtr());

  EdgeOperator PrimalDual;
  PrimalDual.addTerm(new EdgeVec2AndEdgeAtEdges(sphereInstat.getSolPrimal(), sphereInstat.getSolDual(), new Prod2()));
  decSphere.addMatrixOperator(PrimalDual, 0, 1, sphereInstat.getKnPtr());
  decSphere.addMatrixOperator(PrimalDual, 1, 0, sphereInstat.getKnPtr());


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

  //EdgeOperator K;
  //K.addTerm(new IdentityAtEdges(1.0));
  //decSphere.addMatrixOperator(K, 0, 0);
  //decSphere.addMatrixOperator(K, 1, 1);

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

  Timer t;
  sphereInstat.solve();
  cout << "prob solved in " << t.elapsed() << " sec" << endl;
}
