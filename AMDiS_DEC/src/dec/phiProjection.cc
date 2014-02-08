#include "AMDiS.h"
#include "phiProjection.h"
#include "WorldVectorHelper.h"

namespace AMDiS {

  void PhiProject::project(WorldVector<double> &x) {
    double evalPhi = (*phi)(x);
    double evalPhiOld;
    WorldVector<double> xWedge;
    int n = 0;
    int nMax = 20;
    while(abs(evalPhi) > eps && n < nMax) {
      WorldVector<double> evalGradPhi = (*gradPhi)(x);
      xWedge = x - (evalPhi / dot(evalGradPhi, evalGradPhi)) * evalGradPhi;
      //evalPhiOld = evalPhi;

      x = xWedge;
      evalPhi = (*phi)(x);

      n++;
      if (n == nMax) cout << evalPhi << endl;
      //cout << evalPhi << endl;
    }
  }

  //void PhiProject::project(WorldVector<double> &x) {
  //  double evalPhi = (*phi)(x);
  //  WorldVector<double> evalGradPhi = (*gradPhi)(x);
  //  double evalPhiOld;
  //  double f = evalPhi;
  //  double df;
  //  double h = 0;
  //  int n = 0;
  //  WorldVector<double> xWedge;
  //  while(abs(f) > eps && n < 10) {
  //    xWedge = x - h * evalGradPhi;
  //    WorldVector<double> evalGradPhiWedge = (*gradPhi)(xWedge);
  //    f = (*phi)(xWedge);
  //    df = dot(evalGradPhi, evalGradPhiWedge);
  //    h -= f / df;
  //    n++;

  //    //cout << f << endl;

  //  }
  //  x += h * evalGradPhi;
  //}

}
