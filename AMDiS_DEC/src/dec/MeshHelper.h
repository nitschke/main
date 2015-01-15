#ifndef MESH_HELPER
#define MESH_HELPER

#include "AMDiS.h"
#include <iostream>
#include <fstream>


namespace AMDiS {

DOFVector<int> getConnections(const FiniteElemSpace *feSpace);

DOFVector<double> get1RingVols(const FiniteElemSpace *feSpace);

DOFVector<double> getDualVols(const FiniteElemSpace *feSpace);

// 1-Ring approximation, don't need wellcentered
DOFVector<double> getVoronoiRadii(const FiniteElemSpace *feSpace);

// need wellcentered
DOFVector<double> getVoronoiRadiiDualApprox(const FiniteElemSpace *feSpace);

// Area weighted
DOFVector<WorldVector<double> > getNormals(const FiniteElemSpace *feSpace, bool norma = false);

DOFVector<WorldVector<double> > getNormalsVoronoiAverage(const FiniteElemSpace *feSpace, bool norma = false);

DOFVector<WorldVector<double> > getNormalsConnectionAverage(const FiniteElemSpace *feSpace, bool normalize = false);

//normalized
DOFVector<WorldVector<double> > getNormalsAngleAverage(const FiniteElemSpace *feSpace);

//normalized
DOFVector<WorldVector<double> > getNormalsAngleEdgeReciprocalAverage(const FiniteElemSpace *feSpace);

//normalized
DOFVector<WorldVector<double> > getNormalsEdgeReciprocalAverage(const FiniteElemSpace *feSpace);

//normalized; only for positive mean curvature
DOFVector<WorldVector<double> > getNormalsBeltramiAverage(const FiniteElemSpace *feSpace);

DOFVector<WorldVector<double> > getConnectionForces(const FiniteElemSpace *feSpace, bool constantRadii = false, double k = 0.75);

double getMaxMagnitude(DOFVector<WorldVector<double> > F);

//need wellcentered, voronoiarea average
DOFVector<double> getAverage(DOFVector<double> f);

// all neigbours get weight 1 and vertex self 1-Ring*midPointWeight, not scaled! 
DOFVector<WorldVector<double> > getAverage(DOFVector<WorldVector<double> > v, double midPointWeight = 1.0);

// all neigbours get weight 1 and vertex self 1-Ring*midPointWeight, not scaled!; only on defects 
DOFVector<WorldVector<double> > getAverageDefects(DOFVector<WorldVector<double> > v, double midPointWeight = 1.0);

class MeshInfoCSVWriter {

public:
  MeshInfoCSVWriter(string file) : out(file.data(), ios::out), n(0) {
    out << "Nr,AvDiameter,MaxDiameter,AvArea,MinArea,MaxArea,AvMaxAngle,MaxMaxAngle,AvAngleRatio,MaxAngleRatio"<< endl;
  }

  ~MeshInfoCSVWriter() {out.close();}

  void appendData(const FiniteElemSpace *feSpace, bool verbose = false);

private:
  ofstream out;
  int n;

};
  
}

#endif
