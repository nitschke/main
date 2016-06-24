#ifndef EXTREMEVALUETRACKER_H
#define EXTREMEVALUETRACKER_H

#include "Dec_fwd.h"
#include "DecProblemStat.h"

using namespace std;
namespace AMDiS { namespace dec {
class DofVertexVector;
class ExtremeValueTracker {
public:

  ExtremeValueTracker(DecProblemStat *ps) {
    string csvfn;
    Parameters::get(ps->getName() + "->output->filename", csvfn);
    csvfn += "ExtremeValues.csv";
    csvout.open(csvfn.c_str(), ios::out);
    csvout << "t,x,y,z" << endl;
  }

  // for scalar values stored at edge centers-> return number of maximas
  int trackdownMaxima(const DofEdgeVector &dofe, double time, double minVal); 

  ~ExtremeValueTracker() {
    csvout.close();
  }


private:

  bool isMaximum(const DofEdgeVector &dofe, const EdgeElement &eel) const;

  ofstream csvout;

};

}}

#endif
