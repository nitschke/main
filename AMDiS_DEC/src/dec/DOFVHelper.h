#include "AMDiS.h"

namespace AMDiS {

  const DOFVector<double>& minus(const DOFVector<double>& x, const DOFVector<double>& y);

  DOFVector<double> halfMag(const DOFVector<double>& x, const DOFVector<double>& y, const DOFVector<double>& z);

  DOFVector<double> prod01(const DOFVector<WorldVector<double> >& v);

  DOFVector<double> halfSum01(const DOFVector<WorldVector<double> >& v);

  void printError(const DOFVector<double> &dofv,const DOFVector<double> &sol, string name);

}
