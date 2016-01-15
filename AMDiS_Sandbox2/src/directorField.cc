#include <iostream>
#include <fstream>
#include "Dec.h"
#include "ExtremeValueTracker.h"
#include "phiProjection.h"


using namespace std;
using namespace AMDiS;
using namespace dec;

//rotate around y by angle alpha
WorldVector<double> Ry(double alpha, WorldVector<double> X) {
  WorldVector<double> rval;
  double c = cos(alpha);
  double s = sin(alpha);
  rval[0] = c*X[0] + s*X[2];
  rval[1] = X[1];
  rval[2] = -s*X[0] + c*X[2];
  return rval;
}


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

class intValTwoDefectsNonics: public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  intValTwoDefectsNonics() : AbstractFunction<WorldVector<double>, WorldVector<double> >(), mainDir(), extensionDir()
  {
    FUNCNAME("intValTwoDefectsNonics::intValTwoDefectsNonics()");
    mainDir.set(0.0);
    mainDir[0] = 2.0; mainDir[1] = 0.0; mainDir[2] = 1.0;
    
    extensionDir.set(0.0);
    extensionDir[2] = 1.0;
    
    MSG("created director Field with 2 defects located on xy aequator of nonic/octic surface \n");
  }
  WorldVector<double> operator()(const WorldVector<double>& x) const 
  {
    if (x[0] > 1.0 + 1.0e-8)
      return extensionDir;
    
    if (x[0] <= 1.0 + 1.0e-8)
    {
      return mainDir; 
    }
    
    WorldVector<double> nuescht;
    nuescht.set(0.0);
    return nuescht;
  }
private:
  WorldVector<double> mainDir, extensionDir;
};

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


// <dX,[p,q]>
class DX_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DX_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[0] - p[0];
    //return  q[0]*(q[2]+1.0)*(q[2]+1.0) - p[0]*(p[2]+1.0)*(p[2]+1.0);
  }
};

class Df_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Df_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    double a = 0.5*M_PI;
    return  (Ry(a,q))[0] - (Ry(a,p))[0];
  }
};

// <d||X||^2,[p,q]>
class DNorm_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DNorm_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q*q - p*p;
  }
};

// <0,[p,q]>
class Null_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  Null_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  0.0;
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
    string csvNumDefFn = csvfn + "NumberOfDefects.csv";
    csvfn += "Energies.csv";

    csvout.open(csvfn.c_str(), ios::out);
    csvout << "Time,Div,Rot,Norm,Full" << endl;

    csvNumDefout.open(csvNumDefFn.c_str(), ios::out);
    csvNumDefout << "Time,NumberOfDefects" << endl;

    
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

    oldEnergy = 1.e+100;
 }


  void closeTimestep() {
    cout << setprecision(10);
    csvout << setprecision(10);
    double time = t;
    DecProblemInstat::closeTimestep();
    solPrimal = statProb->getSolution(0);
    solDual =  statProb->getSolution(1);
    DofEdgeVectorPD evecPD(solPrimal, solDual);
    normAlpha = evecPD.getNormOnEdges();

    DofEdgeVector scaledNorm(normAlpha);
    scaledNorm.evalFunction(&delta);
    int nmaximas = tracker.trackdownMaxima(scaledNorm, time, 0.1);
    cout << "*** " << nmaximas << " maximas ***" << endl;

    csvNumDefout << time << "," << nmaximas << endl;

    DofEdgeVector normDeviat = normAlpha * normAlpha;
    normDeviat += -1.0;

    double ne =  0.25 * Kn * normDeviat.L2NormSquared(); 
    double dive = evecPD.getDirichletEnergy(-0.5 * MinusK1, 0.0);
    double rote = evecPD.getDirichletEnergy(0.0           , -0.5 * MinusK3);
    double energy = dive + rote + ne;
    csvout << time << "," << dive << "," << rote << "," << ne << ","<< energy << endl;
    
    if (step == 6)  oldEnergy = energy;
    if (step%2 == 1 && step > 6) { 
      cout << "### Energy: " << oldEnergy << " -> " << energy << " ###" << endl;
      double eder = (oldEnergy - energy) / energy;
      cout << "###     rel Diff: " << eder << " ###" << endl;
      cout << "### tau: " << tau << " ###" << endl;

      double eps = 5.E-4; //4.e-5
      double tauMax = 0.5;
      if (eder < eps && tau < tauMax && eder > -1.e-8) {
        t -= tau; // undo in closeTimestep
        tau *= 2.0;
        if (tau > tauMax) tau = tauMax;
        inv_tau = 1. / tau;
        t += tau;
        cout << "### tau -> " << tau << " (coarsening) ###" << endl;
      }

      double eps2 = 5.E-3;
      //double tauMin = 5.e-2;
      double tauMin = 1.e-4;
      if ((eder > eps2 || eder < -1.e-8) && tau > tauMin) {
        t -= tau; // undo in closeTimestep
        tau /= 8.0;
        if (tau < tauMin) tau = tauMin;
        inv_tau = 1. / tau;
        t += tau;
        cout << "### tau -> " << tau << " (refining) ###" << endl;
      }
      oldEnergy = energy;

      if (abs(eder) < 1.0e-14) {
        statProb->writeSolution(t-tau);
        ERROR_EXIT("STAGNATION EXIT\n");
      }
    }
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
  ofstream csvNumDefout;

  Delta delta;
  ExtremeValueTracker tracker;
};



int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  //new PhiProject(1, VOLUME_PROJECTION, new Phi(), new GradPhi(), 1.0e-8);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  int seed = -1;
  Parameters::get("userParameter->seed", seed);
  TEST_EXIT(seed >= 0)("seed must be positive");

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());


  DecProblemStat decSphere(&sphere, edgeMesh);

  DofEdgeVectorPD initSol(edgeMesh, "initSol");
  Noise_d noiseFun(seed);
  //initSol.set(&noiseFun);
  //initSol.set(new DX_d());
  //initSol.set(new Df_d());
  //initSol.set(new DNorm_d());
  //initSol.set(new Null_d());
  //initSol.interpol(new Michael(0.01));
  initSol.interpol(new intValTwoDefectsNonics);
  initSol.normalize(1.E-10);
  initSol.writeSharpOnEdgesFile("output/initSolSharp.vtu");

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

  //EdgeOperator ScalePrimal;
  //ScalePrimal.addTerm(new EdgeVecAtEdges(sphereInstat.getNormPtr(), &vmo, -1.0));
  //ScalePrimal.setUhOld(initPrimal, 0);
  //decSphere.addVectorOperator(ScalePrimal, 0, sphereInstat.getKnPtr());

  //EdgeOperator ScaleDual;
  //ScaleDual.addTerm(new EdgeVecAtEdges(sphereInstat.getNormPtr(), &vmo, -1.0));
  //ScaleDual.setUhOld(initDual, 1);
  //decSphere.addVectorOperator(ScaleDual, 1, sphereInstat.getKnPtr());

  //EdgeOperator Scale;
  //Scale.addTerm(new EdgeVecAtEdges(sphereInstat.getNormPtr(), &vmo));
  //decSphere.addMatrixOperator(Scale, 0, 0,  sphereInstat.getKnPtr());
  //decSphere.addMatrixOperator(Scale, 1, 1, sphereInstat.getKnPtr());

  ////Grad_alpha((||alpha_Old||^2-1)*alpha_Old)*alpha_{|Old}
  //AlphaGrad1 ag1;
  //AlphaGrad2 ag2;

  //EdgeOperator AlphaGradient1Primal;
  //AlphaGradient1Primal.addTerm(new EdgeVec2AndEdgeAtEdges(sphereInstat.getSolPrimal(), sphereInstat.getSolDual(), &ag1));
  //AlphaGradient1Primal.setUhOld(initPrimal, 0);
  //decSphere.addMatrixOperator(AlphaGradient1Primal, 0, 0, sphereInstat.getKnPtr());
  //decSphere.addVectorOperator(AlphaGradient1Primal, 0, sphereInstat.getKnPtr());

  //EdgeOperator AlphaGradient2Primal;
  //AlphaGradient2Primal.addTerm(new EdgeVec2AndEdgeAtEdges(sphereInstat.getSolPrimal(), sphereInstat.getSolDual(), &ag2));
  //AlphaGradient2Primal.setUhOld(initDual, 1);
  //decSphere.addMatrixOperator(AlphaGradient2Primal, 0, 1, sphereInstat.getKnPtr());
  //decSphere.addVectorOperator(AlphaGradient2Primal, 0, sphereInstat.getKnPtr());

  //EdgeOperator AlphaGradient1Dual;
  //AlphaGradient1Dual.addTerm(new EdgeVec2AndEdgeAtEdges(sphereInstat.getSolDual(), sphereInstat.getSolPrimal(), &ag1));
  //AlphaGradient1Dual.setUhOld(initDual, 1);
  //decSphere.addMatrixOperator(AlphaGradient1Dual, 1, 1, sphereInstat.getKnPtr());
  //decSphere.addVectorOperator(AlphaGradient1Dual, 1, sphereInstat.getKnPtr());

  //EdgeOperator AlphaGradient2Dual;
  //AlphaGradient2Dual.addTerm(new EdgeVec2AndEdgeAtEdges(sphereInstat.getSolDual(), sphereInstat.getSolPrimal(), &ag2));
  //AlphaGradient2Dual.setUhOld(initPrimal, 0);
  //decSphere.addMatrixOperator(AlphaGradient2Dual, 1, 0, sphereInstat.getKnPtr());
  //decSphere.addVectorOperator(AlphaGradient2Dual, 1, sphereInstat.getKnPtr());

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
  //decSphere.addMatrixOperator(PrimalPrimal, 1, 1, sphereInstat.getKnPtr());
  
  EdgeOperator DualDual;
  DualDual.addTerm(new EdgeVec2AndEdgeAtEdges(sphereInstat.getSolDual(), sphereInstat.getSolDual(), new Prod2()));
  decSphere.addMatrixOperator(DualDual, 1, 1, sphereInstat.getKnPtr());

  EdgeOperator PrimalDual;
  PrimalDual.addTerm(new EdgeVec2AndEdgeAtEdges(sphereInstat.getSolPrimal(), sphereInstat.getSolDual(), new Prod2()));
  decSphere.addMatrixOperator(PrimalDual, 0, 1, sphereInstat.getKnPtr());
  decSphere.addMatrixOperator(PrimalDual, 1, 0, sphereInstat.getKnPtr());

  //EdgeOperator PrimalDual2;
  //PrimalDual2.addTerm(new EdgeVec2AndEdgeAtEdges(sphereInstat.getSolPrimal(), sphereInstat.getSolDual(), new Prod2()));
  //decSphere.addMatrixOperator(PrimalDual2, 1, 0, sphereInstat.getMinusKnPtr());

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
  //K.addTerm(new IdentityAtEdges(-1.0));
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
