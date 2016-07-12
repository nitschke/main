#include "AMDiS.h"
#include "phiProjection.h"

using namespace std;
namespace AMDiS {

  void PhiProject::project(WorldVector<double> &x) {
    double evalPhi = (*phi)(x);
    double evalPhiOld;
    WorldVector<double> xWedge;
    int n = 0;
    double c = 1.0;
    while(abs(evalPhi) > eps && n < nMax) {
      WorldVector<double> evalGradPhi = (*gradPhi)(x);
      xWedge = x - c * (evalPhi / (evalGradPhi * evalGradPhi)) * evalGradPhi;
      //evalPhiOld = evalPhi;

      x = xWedge;
      evalPhi = (*phi)(x);

      n++;
      if (n == nMax) cout << evalPhi << endl;
      //cout << evalPhi << endl;
    }
  }

  void PhiProject::project(const WorldVector<double> &x, WorldVector<double> &v) {
    WorldVector<double> normal = getNormal(x);
    v -= (v * normal) * normal;
  }


  WorldVector<double> PhiProject::getNormal(const WorldVector<double> &x) {
     WorldVector<double> evalGradPhi = (*gradPhi)(x);
     return (1./std::sqrt(evalGradPhi*evalGradPhi)) * evalGradPhi;
  }

}
