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
  int trackdownMaxima(const DofEdgeVector &dofe, double time, double minVal = 0.0); 

  ~ExtremeValueTracker() {
    csvout.close();
  }


protected:

  ExtremeValueTracker(){}

  bool isMaximum(const DofEdgeVector &dofe, const EdgeElement &eel) const;
  bool isMinimum(const DofEdgeVector &dofe, const EdgeElement &eel) const;

  ofstream csvout;

};

//find the overall minimum of an edgevector for Z > 0 and write out the Z coord
class OneMinValTrackerInPositiveZ : public ExtremeValueTracker {
public:
  
  OneMinValTrackerInPositiveZ(DecProblemStat *ps, bool writeHeader = false) {
    string csvfn;
    Parameters::get(ps->getName() + "->output->filename", csvfn);
    csvfn += "MinValues.csv";
    csvout.open(csvfn.c_str(), ios::out);
    if (writeHeader) csvout << "t,z" << endl;
  }

  int trackdownMinima(const DofEdgeVector &dofe, double time);

};

}}

#endif
