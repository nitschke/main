#ifndef DOFCONDITION
#define DOFCONDITION

#include "AMDiS.h"

using namespace AMDiS;

class ProblemStatExt : public ProblemStat
{
  public:
    ProblemStatExt(ProblemStat prob) : ProblemStat(prob) {};

    void setDOFCondition(int iFESpace, 
                         int iDOF, 
                         double val,
                         bool verbose);


    void updateSolverMatrix()
    {
      solverMatrix.setMatrix(*systemMatrix);
    };
};

#endif
