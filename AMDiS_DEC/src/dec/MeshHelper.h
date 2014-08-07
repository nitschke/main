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

// weighted, not nomalized
DOFVector<WorldVector<double> > getNormals(const FiniteElemSpace *feSpace);

DOFVector<WorldVector<double> > getNormalsConnectionAverage(const FiniteElemSpace *feSpace);

DOFVector<WorldVector<double> > getConnectionForces(const FiniteElemSpace *feSpace, bool constantRadii = false, double k = 0.75);

double getMaxMagnitude(DOFVector<WorldVector<double> > F);

//need wellcentered
DOFVector<double> getAverage(DOFVector<double> f);


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
