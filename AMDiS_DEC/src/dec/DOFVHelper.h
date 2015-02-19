#include "AMDiS.h"

using namespace std;
namespace AMDiS {

  const DOFVector<double>& minus(const DOFVector<double>& x, const DOFVector<double>& y);

  const DOFVector<double>& minus(const DOFVector<WorldVector<double> >& x, const DOFVector<WorldVector<double> >& y);

  DOFVector<double> sqrtScale(const DOFVector<double>& x);

  DOFVector<double> mag(const DOFVector<WorldVector<double> >& v);

  DOFVector<double> halfMag(const DOFVector<double>& x, const DOFVector<double>& y, const DOFVector<double>& z,
                          std::string resName = "HalfOfMagXYZ",
                            bool respOrientation = false);

  DOFVector<double> prod01(const DOFVector<WorldVector<double> >& v);

  DOFVector<double> halfSum01(const DOFVector<WorldVector<double> >& v);

  DOFVector<WorldVector<double> > toWorld(const DOFVector<double>& x, const DOFVector<double>& y, const DOFVector<double>& z);

  DOFVector<double> getComp(int i, const DOFVector<WorldVector<double> >& v, std::string resName = "v_i");

  DOFVector<WorldVector<double> > normalize(const DOFVector<WorldVector<double> >& v);

  void printError(const DOFVector<double> &dofv,const DOFVector<double> &sol, string name);

  void printError(SystemVector &sysv, int i0 ,int i1, int i2, const DOFVector<WorldVector<double> > &sol, string name);

  void printError(const  DOFVector<WorldVector<double> > dofv, const DOFVector<WorldVector<double> > &sol, string name);

}
