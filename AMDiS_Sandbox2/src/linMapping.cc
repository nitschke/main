#include "Dec.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

//eval of rotation on Sphere
class EvalRotSphere : public TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> >
{
  public:

  EvalRotSphere(double phi_ = M_PI/2.) : TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> >(), phi(phi_){}

  double operator()(const WorldVector<double>& v, const WorldVector<double>& w, const WorldVector<double>& coords) const
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldMatrix<double> R;
    R.set(0.0);
    //R[0][0] = cos(phi);
    //R[0][1] = z*sin(phi);
    //R[0][2] = -1.*y*sin(phi);
    //R[1][0] = -1.*z*sin(phi);
    //R[1][1] = cos(phi);
    //R[1][2] = x*sin(phi);
    //R[2][0] = y*sin(phi);
    //R[2][1] = -1.*x*sin(phi);
    //R[2][2] = cos(phi);
    
    R[0][0] = cos(z*phi);
    R[0][1] = z*sin(z*phi);
    R[0][2] = -1.*y*sin(z*phi);
    R[1][0] = -1.*z*sin(z*phi);
    R[1][1] = cos(z*phi);
    R[1][2] = x*sin(z*phi);
    R[2][0] = y*sin(z*phi);
    R[2][1] = -1.*x*sin(z*phi);
    R[2][2] = cos(z*phi);
    return (R*w)*v;
  }

  private:

  double phi;
};

//eval of rotation on Ellipsoid(1,0.5,1.5)
class EvalRotEllipsoid : public TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> >
{
  public:

  EvalRotEllipsoid(double phi_ = M_PI/2.) : TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> >(), phi(phi_){}

  double operator()(const WorldVector<double>& v, const WorldVector<double>& w, const WorldVector<double>& coords) const
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldMatrix<double> R;
    R.set(0.0);

    R[0][0] = cos(phi);
    R[0][1] = (8.*z*sqrt(-1.*(-9. + 4.*pow(z,2))*(81.*pow(x,2) + 16.*(81.*pow(y,2) + pow(z,2))))*sin(phi))/((405. - 243.*pow(x,2) + 972.*pow(y,2) - 148.*pow(z,2))*sqrt(9. - 4.*pow(z,2)));
    R[0][2] = (72.*y*sqrt(-1.*(-9. + 4.*pow(z,2))*(81.*pow(x,2) + 16.*(81.*pow(y,2) + pow(z,2))))*sin(phi))/(sqrt(9. - 4.*pow(z,2))*(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2)));
    R[1][0] = (8.*z*sqrt(-1.*(-9. + 4.*pow(z,2))*(81.*pow(x,2) + 16.*(81.*pow(y,2) + pow(z,2))))*sin(phi))/(sqrt(9. - 4.*pow(z,2))*(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2)));
    R[1][1] = cos(phi);
    R[1][2] = (18.*x*sqrt(-1.*(-9. + 4.*pow(z,2))*(81.*pow(x,2) + 16.*(81.*pow(y,2) + pow(z,2))))*sin(phi))/((405. - 243.*pow(x,2) + 972.*pow(y,2) - 148.*pow(z,2))*sqrt(9. - 4.*pow(z,2)));
    R[2][0] = (72.*y*sqrt(-1.*(-9. + 4.*pow(z,2))*(81.*pow(x,2) + 16.*(81.*pow(y,2) + pow(z,2))))*sin(phi))/((405. - 243.*pow(x,2) + 972.*pow(y,2) - 148.*pow(z,2))*sqrt(9. - 4.*pow(z,2)));
    R[2][1] = (18.*x*sqrt(-1.*(-9. + 4.*pow(z,2))*(81.*pow(x,2) + 16.*(81.*pow(y,2) + pow(z,2))))*sin(phi))/(sqrt(9. - 4.*pow(z,2))*(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2)));
    R[2][2] = cos(phi);

    return (R*w)*v;
  }

  private:

  double phi;
};

//Square of Shape (without outer projections Pi -> S2 = JNu.Pi.JNu)
class S2Ellipsoid : public TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> >
{
  public:

  S2Ellipsoid() : TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> >(){}

  double operator()(const WorldVector<double>& v, const WorldVector<double>& w, const WorldVector<double>& coords) const
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldMatrix<double> S2;
    S2.set(0.0);

    S2[0][0] = (-5184.*(6561. - 5148.*pow(z,2) + 1024.*pow(z,4) + 9.*pow(x,2)*(-729. + 256.*pow(z,2)) - 36.*pow(y,2)*(-729. + 256.*pow(z,2))))/pow(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2),3);
    S2[0][1] = (-373248.*x*y*(-729. + 160.*pow(z,2)))/pow(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2),3);
    S2[0][2] = (124416.*x*z*(-15. + 18.*pow(x,2) - 216.*pow(y,2) + 8.*pow(z,2)))/pow(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2),3);
    S2[1][0] = (-93312.*x*y*(-729. + 160.*pow(z,2)))/pow(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2),3);
    S2[1][1] = (5184.*(-6561. + 2664.*pow(z,2) - 400.*pow(z,4) + 9.*pow(x,2)*(-729. + 100.*pow(z,2)) - 36.*pow(y,2)*(-729. + 100.*pow(z,2))))/pow(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2),3);
    S2[1][2] = (31104.*y*z*(237. + 135.*pow(x,2) - 180.*pow(y,2) - 20.*pow(z,2)))/pow(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2),3);
    S2[2][0] = (279936.*x*z*(-15. + 18.*pow(x,2) - 216.*pow(y,2) + 8.*pow(z,2)))/pow(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2),3);
    S2[2][1] = (279936.*y*z*(237. + 135.*pow(x,2) - 180.*pow(y,2) - 20.*pow(z,2)))/pow(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2),3);
    S2[2][2] = (-11664.*(531. + 81.*pow(x,4) + 2160.*pow(y,2) + 1296.*pow(y,4) - 108.*pow(x,2)*(5. + 18.*pow(y,2)) - 200.*pow(z,2) - 16.*pow(z,4)))/pow(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2),3);

    return (S2*w)*v;
  }
};

//S2.Xv on Ellipsoid
class S2Xv : public AbstractFunction<WorldVector<double>, WorldVector<double> >
{
  public:

  S2Xv() : AbstractFunction<WorldVector<double>, WorldVector<double> >() {}

  WorldVector<double> operator()(const WorldVector<double>& coords) const 
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldVector<double> S2Xv;
    S2Xv.set(0.0);

    S2Xv[0] = (20736.*y*(6561. + 216.*(-23. + 3.*pow(x,2) - 4.*pow(y,2))*pow(z,2) + 928.*pow(z,4)))/pow(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2),3);
    S2Xv[1] = (-5184.*x*(6561. + 135.*(-25. + pow(x,2) - 12.*pow(y,2))*pow(z,2) + 460.*pow(z,4)))/pow(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2),3);
    S2Xv[2] = (-279936.*x*y*z*(-207. + 27.*pow(x,2) - 108.*pow(y,2) + 52.*pow(z,2)))/pow(-405. + 243.*pow(x,2) - 972.*pow(y,2) + 148.*pow(z,2),3);

    return S2Xv;
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


class RHS : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  RHS() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return q[2] - p[2]; // -> dz
  }
};

class SqrtFun : public AbstractFunction<double,double>
{
  public:
    SqrtFun() :  AbstractFunction<double,double>(){}

    double operator()(const double &val) const
    {
      return sqrt(val);
    }
};


int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  DecProblemStat decSphere(&sphere, edgeMesh);

  DofEdgeVectorPD rhsEV(edgeMesh, "rhs");
  //rhsEV.set(new RHS());
  rhsEV.interpol(new S2Xv());
  rhsEV.writeSharpFile("output/S2Xv.vtu",&sphere);

  //TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> > *evalMat = new EvalRotSphere(M_PI/2.);
  //TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> > *evalMat = new EvalRotEllipsoid(M_PI/2.);
  TertiaryAbstractFunction<double, WorldVector<double>, WorldVector<double>, WorldVector<double> > *evalMat = new S2Ellipsoid();

  //EdgeOperator funPP;
  //funPP.addTerm(new EdgeFunAtEdges(new MatPP(evalMat)));
  //decSphere.addMatrixOperator(funPP, 0, 0);

  //EdgeOperator funDD;
  //funDD.addTerm(new EdgeFunAtEdges(new MatDD(evalMat)));
  //decSphere.addMatrixOperator(funDD, 1, 1);
  //
  //EdgeOperator funPD;
  //funPD.addTerm(new EdgeFunAtEdges(new MatPD(evalMat)));
  //decSphere.addMatrixOperator(funPD, 0, 1);

  //EdgeOperator funDP;
  //funDP.addTerm(new EdgeFunAtEdges(new MatDP(evalMat)));
  //decSphere.addMatrixOperator(funDP, 1, 0);

  //Alternative
  DofEdgeVector MatPPEV(edgeMesh, "MatPP");
  DofEdgeVector MatDDEV(edgeMesh, "MatDD");
  DofEdgeVector MatPDEV(edgeMesh, "MatPD");
  DofEdgeVector MatDPEV(edgeMesh, "MatDP");
  MatPPEV.set(new MatPP(evalMat));
  MatDDEV.set(new MatDD(evalMat));
  MatPDEV.set(new MatPD(evalMat));
  MatDPEV.set(new MatDP(evalMat));
  Matrix<DofEdgeVector* > MatEVs(2,2);
  MatEVs[0][0] = &MatPPEV;
  MatEVs[1][1] = &MatDDEV;
  MatEVs[0][1] = &MatPDEV;
  MatEVs[1][0] = &MatDPEV;

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      EdgeOperator *fun = new EdgeOperator();
      fun->addTerm(new EdgeVecAtEdges(MatEVs[i][j]));
      decSphere.addMatrixOperator(*fun, i, j);
    }
  }

  //RHS
  EdgeOperator rhsP;
  rhsP.addTerm(new EdgeVecAtEdges(&rhsEV));
  decSphere.addVectorOperator(rhsP, 0);

  DofEdgeVector rhsDualEV = rhsEV.getDual();
  EdgeOperator rhsD;
  rhsD.addTerm(new EdgeVecAtEdges(&rhsDualEV));
  decSphere.addVectorOperator(rhsD, 1);

  decSphere.assembleSystem();

  decSphere.solve();

  decSphere.writeSolution();

  Vector<DofEdgeVector > sols(2);
  sols[0] = decSphere.getSolution(0);
  sols[1] = decSphere.getSolution(1);

  DofEdgeVector solMatsol(edgeMesh, "MatPP");
  solMatsol.set(0.0);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      solMatsol += sols[i].getLocalPDSharp() * sols[j] * (*MatEVs[i][j]);
    }
  }
  solMatsol.writeFile("output/SMS.vtu");
  
  cout << setprecision(10);
  cout << "Int(v.M.v) = " << solMatsol.surfaceIntegration() << endl;

  //using namespace mtl;
  //dense2D<double> mat(60,60);
  //mat = decSphere.getSysMat();
  //edgeMesh->printVolInfos(true,true);
  //cout << sub_matrix(mat,0,30,0,30) << endl << endl;
  //cout << sub_matrix(mat,0,30,30,60) << endl << endl;
  //cout << sub_matrix(mat,30,60,0,30) << endl << endl;
  //cout << sub_matrix(mat,30,60,30,60) << endl << endl;
  //cout << decSphere.getRhs() << endl;


  AMDiS::finalize();
}
