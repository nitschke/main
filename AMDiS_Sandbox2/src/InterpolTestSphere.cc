#include <iostream>
#include <fstream>
#include "Dec.h"

using namespace std;
using namespace AMDiS;
using namespace dec;


class MRotZVec : public AbstractFunction<WorldVector<double>, WorldVector<double> > {

public:
  MRotZVec() : AbstractFunction<WorldVector<double>, WorldVector<double> >() {}

  WorldVector<double> operator()(const WorldVector<double>& coords) const {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldVector<double> r;
    r[0] = -y;
    r[1] = x;
    r[2] = 0.0;
    return r;
  }
};


class DZ_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DZ_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[2] - p[2];
  }
};

class MRotZ2Vec : public AbstractFunction<WorldVector<double>, WorldVector<double> > {

public:
  MRotZ2Vec() : AbstractFunction<WorldVector<double>, WorldVector<double> >() {}

  WorldVector<double> operator()(const WorldVector<double>& coords) const {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    WorldVector<double> r;
    r[0] = -2.0*y*z;
    r[1] = 2.0*x*z;
    r[2] = 0.0;
    return r;
  }
};


class DZ2_d : public BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >
{
public:
  DZ2_d() : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double> >() {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& p, const WorldVector<double>& q) const 
  {
    return  q[2]*q[2] - p[2]*p[2];
  }
};


int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  ///// ===== create projection =====
  ///WorldVector<double> ballCenter;
  ///ballCenter.set(0.0);
  ///new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  EdgeMesh *edgeMesh = new EdgeMesh(sphere.getFeSpace());

  // alpha^# = -RotZ = -(*dZ)^# => *alpha = dZ
  DofEdgeVectorPD mrotz(edgeMesh, "mrotz");
  mrotz.interpol(new MRotZ2Vec());
  mrotz.writeSharpOnEdgesFile("output/mrotzSharp.vtu");

  // dz from interpolated mrotz
  DofEdgeVector dzInterpol = mrotz.getDual();

  // exact dz
  DofEdgeVector dz(edgeMesh, "dz");
  dz.set(new DZ2_d());

  DofEdgeVector errform = dz - dzInterpol;
  errform.writeFile("output/errform.vtu");

  double errMax = errform.absMax();
  double errl2 = errform.l2Norm();

  cout << "ErrMax: " << errMax << endl;
  cout << "Errl2:  " << errl2 << endl;
}
