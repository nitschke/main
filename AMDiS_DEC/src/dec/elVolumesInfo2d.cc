#include "elVolumesInfo2d.h"
#include "WorldVectorHelper.h"

namespace AMDiS {

  ElVolumesInfo2d::ElVolumesInfo2d(const ElInfo *el) {
    dualVertexVol = WorldVector<double>();
    oppEdgeLen = WorldVector<double>();
    dualOppEdgeLen = WorldVector<double>();

    WorldVector<double> p0 = el->getCoord(0);
    WorldVector<double> p1 = el->getCoord(1);
    WorldVector<double> p2 = el->getCoord(2);

    WorldVector<double> e0 = p2 - p1;
    WorldVector<double> e1 = p0 - p2;
    WorldVector<double> e2 = p1 - p0;

    double d = 2.0 * el->getDet() * el->getDet();

    // bary coords of the circumcenter
    double a1 = dot(e1, e1) * dot(e2, e0) / d;
    double a2 = dot(e2, e2) * dot(e1, e0) / d;
    //double a0 = - dot(e0, e0) * dot(e1, e2) / d;
    //double a1 = - dot(e1, e1) * dot(e0, e2) / d;
    //double a2 = - dot(e2, e2) * dot(e0, e1) / d;

    // circumcenter
    WorldVector<double> cc = p0 - a1*e2 + a2*e1;
    //WorldVector<double> cc = a0*p0 + a1*p1 + a2*p2; 
    //appendToCSV(cc, "circumcenter.csv");

    // edge circumcenters
    WorldVector<double> cc0 = p1 + 0.5*e0;
    WorldVector<double> cc1 = p2 + 0.5*e1;
    WorldVector<double> cc2 = p0 + 0.5*e2;

    // dual edges (in T)
    WorldVector<double> starE0 = cc - cc0; 
    WorldVector<double> starE1 = cc - cc1; 
    WorldVector<double> starE2 = cc - cc2; 

    // edge len of the opp vertex
    oppEdgeLen[0] = sqrt(dot(e0, e0));
    oppEdgeLen[1] = sqrt(dot(e1, e1));
    oppEdgeLen[2] = sqrt(dot(e2, e2));

    // dual edge len of the opp vertex
    dualOppEdgeLen[0] = sqrt(dot(starE0, starE0));
    dualOppEdgeLen[1] = sqrt(dot(starE1, starE1));
    dualOppEdgeLen[2] = sqrt(dot(starE2, starE2));

    // Vol of the dual vertex (voronoi cell) (in T)
    dualVertexVol[0] = 0.25 * (dualOppEdgeLen[1]*oppEdgeLen[1] + dualOppEdgeLen[2]*oppEdgeLen[2]);
    dualVertexVol[1] = 0.25 * (dualOppEdgeLen[0]*oppEdgeLen[0] + dualOppEdgeLen[2]*oppEdgeLen[2]);
    dualVertexVol[2] = 0.25 * (dualOppEdgeLen[1]*oppEdgeLen[1] + dualOppEdgeLen[0]*oppEdgeLen[0]);

    //double vol = dualVertexVol[1] + dualVertexVol[2] +  dualVertexVol[0];
    //cout << "\n*** VolErr ***" << endl;
    //cout << abs(2.0*vol - el->getDet()) << endl;
    //cout << "*** Coords ***" << endl;
    //cout << p0 << endl;
    //cout << p1 << endl;
    //cout << p2 << endl;
    //cout << "*** Circumc ***" << endl;
    //cout << cc << endl;
  }
}