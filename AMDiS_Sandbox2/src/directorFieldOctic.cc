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

    oldEnergy = 1.e+100;
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
    int nmaximas = tracker.trackdownMaxima(scaledNorm, time, 0.1);
    cout << "*** " << nmaximas << " maximas ***" << endl;

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

      double eps = 4.E-4;
      double tauMax = 0.64; // = 0.01*2^6
      if (eder < eps && tau < tauMax && eder > -1.e-8) {
        t -= tau; // undo in closeTimestep
        tau *= 2.0;
        if (tau > tauMax) tau = tauMax;
        inv_tau = 1. / tau;
        t += tau;
        cout << "### tau -> " << tau << " (coarsening) ###" << endl;
      }

      double eps2 = 4.E-3;
      double tauMin = 1.5625e-4; //=0.01*2^(-6)
      if ((eder > eps2 || eder < -1.e-8) && tau > tauMin) {
        t -= tau; // undo in closeTimestep
        tau /= 8.0;
        if (tau < tauMin) tau = tauMin;
        inv_tau = 1. / tau;
        t += tau;
        cout << "### tau -> " << tau << " (refining) ###" << endl;
      }
      oldEnergy = energy;
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

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());


  DecProblemStat decSphere(&sphere, edgeMesh);

  DofEdgeVectorPD initSol(edgeMesh, "initSol");
  Noise_d noiseFun(42);
  //initSol.set(&noiseFun);
  initSol.set(new DX_d());
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
