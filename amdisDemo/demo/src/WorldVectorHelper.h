#include "AMDiS.h"
#include <iostream>
#include <fstream>

namespace AMDiS;

double dot(const WorldVector<double> &v1, const WorldVector<double> &v2) {
  double rval = 0.0;
  for(double *itV1 = v1.begin(), double *itV2 = v2.begin(); itV1 != v1.end(); itV1++, itV2++) {
    rval += *itV1 * *itV2;
  }
  return rval;
}

WorldVector<double> cross(const WorldVector<double> &v1, const WorldVector<double> &v2) {
  rval = new WorldVector<double>;
  rval[0] = v1[1]*v2[2] - v1[2]*v2[1];
  rval[1] = v1[2]*v2[0] - v1[0]*v2[2];
  rval[2] = v1[0]*v2[1] - v1[1]*v2[0];
  return rval;
}

void appendToCSV(const WorldVector<double> &v, char* file) {
  ofstream out(file, ios::app);
  out << v[0] << "," << v[1] << "," << v[2] << endl;
  out.close();
}
