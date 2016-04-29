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

//Square of Shape 
class S2Torus : public TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> >
{
  public:

  S2Torus() : TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> >()
  {
  }

  double operator()(const WorldVector<double>& v, const WorldVector<double>& w, const WorldVector<double>& coords) const
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldMatrix<double> S2;
    
    if (x*x + z*z > 4) {
      S2[0][0] = (16.*(4.*pow(x,2)*pow(y,2)*pow(4. + sqrt(1. - 4.*pow(y,2)),2)*pow(17. - 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)),2) + pow(-8. + 32.*pow(y,2) + sqrt(1. - 4.*pow(y,2))*(-17. + 4.*pow(y,2)),2)*pow(z,2)))/(pow(4. + sqrt(1. - 4.*pow(y,2)),2)*pow(17. - 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)),3));
      S2[0][1] = (16.*x*y*(-1. + 4.*pow(y,2) - 4.*sqrt(1. - 4.*pow(y,2))))/(17. - 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)));
      S2[0][2] = (16.*x*(-353. + 34512.*pow(y,2) - 70752.*pow(y,4) + 15616.*pow(y,6) - 256.*pow(y,8) + 16.*sqrt(1. - 4.*pow(y,2))*(-17. + 1934.*pow(y,2) - 1344.*pow(y,4) + 96.*pow(y,6)))*z)/(pow(4. + sqrt(1. - 4.*pow(y,2)),2)*pow(17. - 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)),3));
      S2[1][0] = (16.*x*y*(-1. + 4.*pow(y,2) - 4.*sqrt(1. - 4.*pow(y,2))))/(17. - 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)));
      S2[1][1] = 4. - 16.*pow(y,2);
      S2[1][2] = (16.*y*(-1. + 4.*pow(y,2) - 4.*sqrt(1. - 4.*pow(y,2)))*z)/(17. - 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)));
      S2[2][0] = (16.*x*(-353. + 34512.*pow(y,2) - 70752.*pow(y,4) + 15616.*pow(y,6) - 256.*pow(y,8) + 16.*sqrt(1. - 4.*pow(y,2))*(-17. + 1934.*pow(y,2) - 1344.*pow(y,4) + 96.*pow(y,6)))*z)/(pow(4. + sqrt(1. - 4.*pow(y,2)),2)*pow(17. - 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)),3));
      S2[2][1] = (16.*y*(-1. + 4.*pow(y,2) - 4.*sqrt(1. - 4.*pow(y,2)))*z)/(17. - 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)));
      S2[2][2] = (16.*(pow(x,2)*pow(-8. + 32.*pow(y,2) + sqrt(1. - 4.*pow(y,2))*(-17. + 4.*pow(y,2)),2) + 4.*pow(y,2)*pow(4. + sqrt(1. - 4.*pow(y,2)),2)*pow(17. - 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)),2)*pow(z,2)))/(pow(4. + sqrt(1. - 4.*pow(y,2)),2)*pow(17. - 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)),3));
    } else {
      S2[0][0] = (16.*(4.*pow(x,2)*pow(y,2)*(-8177. + 64.*pow(y,6) + 7448.*sqrt(1. - 4.*pow(y,2)) + pow(y,2)*(17292. - 5312.*sqrt(1. - 4.*pow(y,2))) + 48.*pow(y,4)*(-81. + 8.*sqrt(1. - 4.*pow(y,2)))) + (-1. + 4.*pow(y,2))*(353. + 16.*pow(y,4) - 272.*sqrt(1. - 4.*pow(y,2)) + 8.*pow(y,2)*(-49. + 8.*sqrt(1. - 4.*pow(y,2))))*pow(z,2)))/(pow(-4. + sqrt(1. - 4.*pow(y,2)),2)*pow(-17. + 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)),3));
      S2[0][1] = (-16.*x*y*(-1. + 4.*pow(y,2) + 4.*sqrt(1. - 4.*pow(y,2))))/(-17. + 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)));
      S2[0][2] = (16.*x*(-353. + 34512.*pow(y,2) - 70752.*pow(y,4) + 15616.*pow(y,6) - 256.*pow(y,8) - 16.*sqrt(1. - 4.*pow(y,2))*(-17. + 1934.*pow(y,2) - 1344.*pow(y,4) + 96.*pow(y,6)))*z)/(pow(17. - 4.*pow(y,2) - 8.*sqrt(1. - 4.*pow(y,2)),3)*pow(-4. + sqrt(1. - 4.*pow(y,2)),2));
      S2[1][0] = (-16.*x*y*(-1. + 4.*pow(y,2) + 4.*sqrt(1. - 4.*pow(y,2))))/(-17. + 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)));
      S2[1][1] = 4. - 16.*pow(y,2);
      S2[1][2] = (-16.*y*(-1. + 4.*pow(y,2) + 4.*sqrt(1. - 4.*pow(y,2)))*z)/(-17. + 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)));
      S2[2][0] = (16.*x*(-353. + 34512.*pow(y,2) - 70752.*pow(y,4) + 15616.*pow(y,6) - 256.*pow(y,8) - 16.*sqrt(1. - 4.*pow(y,2))*(-17. + 1934.*pow(y,2) - 1344.*pow(y,4) + 96.*pow(y,6)))*z)/(pow(17. - 4.*pow(y,2) - 8.*sqrt(1. - 4.*pow(y,2)),3)*pow(-4. + sqrt(1. - 4.*pow(y,2)),2));
      S2[2][1] = (-16.*y*(-1. + 4.*pow(y,2) + 4.*sqrt(1. - 4.*pow(y,2)))*z)/(-17. + 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)));
      S2[2][2] = (16.*(pow(x,2)*(-1. + 4.*pow(y,2))*(353. + 16.*pow(y,4) - 272.*sqrt(1. - 4.*pow(y,2)) + 8.*pow(y,2)*(-49. + 8.*sqrt(1. - 4.*pow(y,2)))) + 4.*pow(y,2)*(-8177. + 64.*pow(y,6) + 7448.*sqrt(1. - 4.*pow(y,2)) + pow(y,2)*(17292. - 5312.*sqrt(1. - 4.*pow(y,2))) + 48.*pow(y,4)*(-81. + 8.*sqrt(1. - 4.*pow(y,2))))*pow(z,2)))/(pow(-4. + sqrt(1. - 4.*pow(y,2)),2)*pow(-17. + 4.*pow(y,2) + 8.*sqrt(1. - 4.*pow(y,2)),3));
    }

    return (S2*w)*v;
  }

};

class MatPP : public AbstractFunction<double, EdgeElement>
{
  public:

  MatPP(TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> > *evalMat_) 
      : AbstractFunction<double, EdgeElement>(), evalMat(evalMat_) {}

  double operator()(const EdgeElement &eel) const
  {
    WorldVector<double> eP = eel.infoLeft->getEdge(eel.dofEdge);
    double lenP2 = eP*eP;
    WorldVector<double> ce = eel.infoLeft->getEdgeCenter(eel.dofEdge);
    return (*evalMat)(eP,eP,ce)/lenP2;
  }

  private:
   
  TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> > *evalMat;
};

class MatDD : public AbstractFunction<double, EdgeElement>
{
  public:

  MatDD(TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> > *evalMat_) 
      : AbstractFunction<double, EdgeElement>(), evalMat(evalMat_) {}

  double operator()(const EdgeElement &eel) const
  {
    WorldVector<double> eD = eel.infoLeft->getCircumcenter() -  eel.infoRight->getCircumcenter();
    double lenD2 = eD*eD;
    WorldVector<double> ce = eel.infoLeft->getEdgeCenter(eel.dofEdge);
    return (*evalMat)(eD,eD,ce)/lenD2;
  }

  private:
   
  TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> > *evalMat;
};

class MatPD : public AbstractFunction<double, EdgeElement>
{
  public:

  MatPD(TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> > *evalMat_) 
      : AbstractFunction<double, EdgeElement>(), evalMat(evalMat_) {}

  double operator()(const EdgeElement &eel) const
  {
    WorldVector<double> eD = eel.infoLeft->getCircumcenter() -  eel.infoRight->getCircumcenter();
    WorldVector<double> eP = eel.infoLeft->getEdge(eel.dofEdge);
    double lenD = sqrt(eD*eD);
    double lenP = sqrt(eP*eP);
    WorldVector<double> ce = eel.infoLeft->getEdgeCenter(eel.dofEdge);

    return -(*evalMat)(eP,eD,ce)/(lenP*lenD);
  }

  private:
   
  TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> > *evalMat;
};

class MatDP : public AbstractFunction<double, EdgeElement>
{
  public:

  MatDP(TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> > *evalMat_) 
      : AbstractFunction<double, EdgeElement>(), evalMat(evalMat_) {}

  double operator()(const EdgeElement &eel) const
  {
    WorldVector<double> eD = eel.infoLeft->getCircumcenter() -  eel.infoRight->getCircumcenter();
    WorldVector<double> eP = eel.infoLeft->getEdge(eel.dofEdge);
    double lenD = sqrt(eD*eD);
    double lenP = sqrt(eP*eP);
    WorldVector<double> ce = eel.infoLeft->getEdgeCenter(eel.dofEdge);

    return -(*evalMat)(eD,eP,ce)/(lenP*lenD);
  }

  private:
   
  TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> > *evalMat;
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

class P0: public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
public:
  P0(double largeCirc, double smallCirc) : AbstractFunction<WorldVector<double>, WorldVector<double> >(), A(largeCirc), B(smallCirc)
  {
  }

  WorldVector<double> operator()(const WorldVector<double>& coords) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldVector<double> p0;
    p0.set(0.0);

    if (x*x + z*z > 4) {
      p0[0] = (-4.*(2.*B*x*y + A*z))/pow(4. + sqrt(1. - 4.*pow(y,2)),2);
      p0[1] = (2.*B*sqrt(1. - 4.*pow(y,2)))/(4. + sqrt(1. - 4.*pow(y,2)));
      p0[2] = (4.*A*x - 8.*B*y*z)/pow(4. + sqrt(1. - 4.*pow(y,2)),2);
    } else {
      p0[0] = (-4.*(2.*B*x*y + A*z))/pow(-4. + sqrt(1. - 4.*pow(y,2)),2);
      p0[1] = (2.*B*sqrt(1. - 4.*pow(y,2)))/(-4. + sqrt(1. - 4.*pow(y,2)));
      p0[2] = (4.*A*x - 8.*B*y*z)/pow(-4. + sqrt(1. - 4.*pow(y,2)),2);
    }

    return p0;
  }

private:
  
  double A;
  double B;
};


// <dX,[p,q]>
class DX_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DX_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {
    FUNCNAME("DX_d:: DX_d()");
    MSG("created x-unity vector, result in director Field with 4 defects on nonic surface \n");
  }

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
    //double a = 0.5*M_PI;
    //return  (Ry(a,q))[0] - (Ry(a,p))[0];
    return q[1] - p[1];
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
        delta(0.01),
        B2EVs(2,2)
    {
    FUNCNAME("MyInstat::MyInstat(...)");
    string csvfn;
    Parameters::get(probStat->getName() + "->output->filename", csvfn);
    ////forstandalone
    //string divAnim = csvfn + "Divergence.pvd";
    //string rotAnim = csvfn + "Rotation.pvd";
    //divWriter = new AnimationWriter(divAnim);
    //rotWriter = new AnimationWriter(rotAnim);
    //////////////
    string csvNumDefFn = csvfn + "NumberOfDefects.csv";
    csvfn += "Energies.csv";

    csvout.open(csvfn.c_str(), ios::out);
    csvout << "Time,Div,Rot,B2,Norm,Full" << endl;

    csvNumDefout.open(csvNumDefFn.c_str(), ios::out);
    csvNumDefout << "Time,NumberOfDefects" << endl;

    
    K0 = -1.0;
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

    MSG("Init B^2 on all edges: \n");
    double B = 0.0;
    double C = 0.0;
    Parameters::get("userParameter->B", B);
    Parameters::get("userParameter->B", C);
    TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> > *evalMat = new S2Torus();
    B2EVs[0][0] = new DofEdgeVector(probStat->getMesh(),"B2EV PP-Component");
    B2EVs[1][1] = new DofEdgeVector(probStat->getMesh(),"B2EV DD-Component");
    B2EVs[0][1] = new DofEdgeVector(probStat->getMesh(),"B2EV PD-Component");
    B2EVs[1][0] = new DofEdgeVector(probStat->getMesh(),"B2EV DP-Component");
    B2EVs[0][0]->set(new MatPP(evalMat));
    B2EVs[1][1]->set(new MatDD(evalMat));
    B2EVs[0][1]->set(new MatPD(evalMat));
    B2EVs[1][0]->set(new MatDP(evalMat));
    MSG("Done!\n");
 }


  void closeTimestep() {
    cout << setprecision(10);
    csvout << setprecision(10);
    double time = t;
    DecProblemInstat::closeTimestep();
    solPrimal = statProb->getSolution(0);
    solDual =  statProb->getSolution(1);

    ////standalone beginning
    //if (writeSolutions && step%writeEveryithTimestep == 1) {
    //  DOFVector<double> div = solPrimal.divergence();
    //  DOFVector<double> rot = solDual.divergence();
    //  rot *= -1.0;
    //  string fndiv;
    //  string fnrot;
    //  string bn;
    //  Parameters::get(statProb->getName() + "->output->filename", bn);
    //  int prec = 3;
    //  Parameters::get("sphere->output->index precision", prec);
    //  ostringstream timeoss;
    //  timeoss << setprecision(prec) << fixed << time;
    //  fndiv = bn + "Divergence." + timeoss.str() + ".vtu";
    //  fnrot = bn + "Rotation." + timeoss.str() + ".vtu";
    //  io::VtkVectorWriter::writeFile(div,fndiv);
    //  io::VtkVectorWriter::writeFile(rot,fnrot);
    //  divWriter->updateAnimationFile(time,fndiv);
    //  rotWriter->updateAnimationFile(time,fnrot);
    //}
    ////standalone end
    
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

    Vector<DofEdgeVector* > sols(2);
    sols[0] = &solPrimal;
    sols[1] = &solDual;
    DofEdgeVector solB2sol(statProb->getMesh(), "pB2p");
    solB2sol.set(0.0);
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        solB2sol += (*sols[i]).getLocalPDSharp() * (*sols[j]) * (*B2EVs[i][j]);
      }
    }
    double B2e = K0 * solB2sol.surfaceIntegration();


    double energy = dive + rote + B2e + ne ;
    csvout << time << "," << dive << "," << rote << "," << B2e << "," << ne << ","<< energy << endl;
    

    double epsCoarse = -1./0.; // no coarsing
    double tauMax = -1.; // no coarsing
    double facCoarse = 2.0;
    Parameters::get("userParameter->adapt->epsCoarse", epsCoarse);
    Parameters::get("userParameter->adapt->tauMax", tauMax);
    Parameters::get("userParameter->adapt->facCoarse", facCoarse);
    
    
    double epsRefine = 1./0.; // no refining
    double tauMin = 1./0.; // no refining
    double facRefine = 1./8.;
    Parameters::get("userParameter->adapt->epsRefine", epsRefine);
    Parameters::get("userParameter->adapt->tauMin", tauMin);
    Parameters::get("userParameter->adapt->facRefine", facRefine);

    double epsStagnation = 1.e-14;
    Parameters::get("userParameter->trunc->epsStagnation", epsStagnation);

    int firstStep = 6;
    int stepModulo =2;
    Parameters::get("userParameter->adapt->first test at timestep", firstStep);
    Parameters::get("userParameter->adapt->test every i-th timestep", stepModulo);

    if (step == firstStep)  oldEnergy = energy;
    if (step%stepModulo == 1 && step > firstStep) { 
      cout << "### Energy: " << oldEnergy << " -> " << energy << " ###" << endl;
      double eder = (oldEnergy - energy) / energy;
      cout << "###     rel Diff: " << eder << " ###" << endl;
      cout << "### tau: " << tau << " ###" << endl;

      if (eder < epsCoarse && tau < tauMax && eder > -1.e-8) {
        t -= tau; // undo in closeTimestep
        tau *= facCoarse;
        if (tau > tauMax) tau = tauMax;
        inv_tau = 1. / tau;
        t += tau;
        cout << "### tau -> " << tau << " (coarsening) ###" << endl;
      }

      if ((eder > epsRefine || eder < -1.e-8) && tau > tauMin) {
        t -= tau; // undo in closeTimestep
        tau *= facRefine;
        if (tau < tauMin) tau = tauMin;
        inv_tau = 1. / tau;
        t += tau;
        cout << "### tau -> " << tau << " (refining) ###" << endl;
      }
      oldEnergy = energy;

      if (abs(eder) < epsStagnation) {
        statProb->writeSolution(t-tau);
        ERROR_EXIT("STAGNATION EXIT\n");
      }
    }
  }

  DofEdgeVector* getB2PDComponent(int i, int j)
  {
    return B2EVs[i][j];
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

  double* getK0Ptr() {
    return &K0;
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

  double K0;
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

  Matrix<DofEdgeVector* > B2EVs; //Co-Contra-variant PD-Tensor on all edges
  
  ////for standalone part
  //AnimationWriter *divWriter;
  //AnimationWriter *rotWriter;
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

  std::string initFun = "noise";
  Parameters::get("userParameter->initField", initFun);

  double EYRotatedAngle = 0.05;
  Parameters::get("userParameter->rotated_ey->angle", EYRotatedAngle);


  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());


  DecProblemStat decSphere(&sphere, edgeMesh);

  DofEdgeVectorPD initSol(edgeMesh, "initSol");

  //if (initFun == "noise") {
  //  Noise_d noiseFun(seed);
  //  initSol.set(&noiseFun);
  //} 
  //else if (initFun == "ex") initSol.set(new DX_d());
  //else if (initFun == "rotated_ey") initSol.interpol(new EYRotated(EYRotatedAngle));
  ////else if (initFun == "rotated_ey") initSol.interpol(new EYRotated(0.0));
  //else ERROR_EXIT("Don't know this userParameter->initField");
  //initSol.set(new Df_d());
  //initSol.set(new DNorm_d());
  //initSol.set(new Null_d());
  //initSol.interpol(new Michael(0.01));
  initSol.interpol(new P0(0.01,1.0));
  //


  //initSol.normalize(1.E-10);
  //initSol.writeSharpOnEdgesFile("output/initSolSharp.vtu");


  MyInstat sphereInstat(&decSphere, initSol);

  //B^2
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      EdgeOperator *B2 = new EdgeOperator();
      B2->addTerm(new EdgeVecAtEdges(sphereInstat.getB2PDComponent(i,j)));
      decSphere.addMatrixOperator(*B2, i, j, sphereInstat.getK0Ptr());
    }
  }

  //Laplace
  EdgeOperator LaplaceB;
  LaplaceB.addTerm(new LaplaceBeltramiAtEdges());
  decSphere.addMatrixOperator(LaplaceB, 0, 0, sphereInstat.getMinusK3Ptr());
  decSphere.addMatrixOperator(LaplaceB, 1, 1, sphereInstat.getMinusK1Ptr());

  EdgeOperator LaplaceCB;
  LaplaceCB.addTerm(new LaplaceCoBeltramiAtEdges());
  decSphere.addMatrixOperator(LaplaceCB, 0, 0, sphereInstat.getMinusK1Ptr());
  decSphere.addMatrixOperator(LaplaceCB, 1, 1, sphereInstat.getMinusK3Ptr());

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


  Timer t;
  sphereInstat.solve();
  cout << "prob solved in " << t.elapsed() << " sec" << endl;
}
