#ifndef WORLD_VECTOR_HELPER
#define WORLD_VECTOR_HELPER

#include "AMDiS.h"
#include <iostream>
#include <fstream>

using namespace std;
namespace AMDiS {

inline double dot(const WorldVector<double> &v1, const WorldVector<double> &v2) {
  double rval = 0.0;
  for(const double *itV1 = v1.begin(), *itV2 = v2.begin(); itV1 != v1.end(); itV1++, itV2++) {
    rval += *itV1 * *itV2;
  }
  return rval;
}

inline double wvnorm(const WorldVector<double> &v1) {
  return std::sqrt(dot(v1,v1));
}

inline WorldVector<double> cross(const WorldVector<double> &v1, const WorldVector<double> &v2) {
  WorldVector<double> rval;
  rval[0] = v1[1]*v2[2] - v1[2]*v2[1];
  rval[1] = v1[2]*v2[0] - v1[0]*v2[2];
  rval[2] = v1[0]*v2[1] - v1[1]*v2[0];
  return rval;
}

inline void appendToCSV(const WorldVector<double> &v, string file) {
  ofstream out(file.data(), ios::app);
  out << v[0] << "," << v[1] << "," << v[2] << endl;
  out.close();
}

}

#endif
