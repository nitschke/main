#include "AMDiS.h"

namespace AMDiS {

inline int matIndex(int z, int s) {return 3 * z + s;}

inline int matIndexZ(int pos) {return pos / 3;}
inline int matIndexS(int pos) {return pos % 3;}

DOFVector<WorldVector<double> > getEigenVals(WorldMatrix< DOFVector<double> * > II);

}
