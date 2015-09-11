#ifndef DECPROBLEMINSTAT_H
#define DECPROBLEMINSTAT_H

#include "Dec_fwd.h"
#include "DecProblemStat.h"

using namespace std;
namespace AMDiS { namespace dec {

//TODO: make it more AMDiS-like
class DecProblemInstat {
public:

  DecProblemInstat(DecProblemStat *probStat, DecProblemStat *initStatProb = NULL);

  virtual void oneIteration();

  void solve();

  virtual void initTimestep();

  virtual void closeTimestep();

  double *getTimePtr() {
    return &t;
  }

  double *getTauPtr() {
    return &tau;
  }

  double *getInvTauPtr() {
    return &inv_tau;
  }

  unsigned int getTimestep() {
    return step;
  }

private:

void updateUhOlds();
void updateUhOlds_EdgeOperators(list<pair<DecOperator*, double*> > &ops, int i);
void updateUhOlds_VertexOperators(list<pair<DecOperator*, double*> > &ops, int i);

protected:

DecProblemStat *statProb;
DecProblemStat *initProb;

double t0; // start time
double t1; // end time
double t; //current time
double tau; //time step size
unsigned int step; // current time step

Timer timer;

double inv_tau; // 1 / tau

bool writeSolutions;
unsigned int writeEveryithTimestep;
};

}}
#endif
