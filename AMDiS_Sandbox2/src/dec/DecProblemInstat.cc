#include "DecProblemInstat.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

DecProblemInstat::DecProblemInstat(DecProblemStat *probStat, DecProblemStat *initStatProb) 
        : statProb(probStat),  initProb(initStatProb) {
  FUNCNAME("DecProblemInstat::DecProblemInstat");
  string keyprefix = statProb->getProblemStat()->getName() + "->adapt";
  
  t0 = 0.0;
  Parameters::get(keyprefix + "->start time", t0);

  t1 = 1.0;
  Parameters::get(keyprefix + "->end time", t1);
  TEST_EXIT(t0 <= t1)("Start time must be lower then end time.");

  tau = 0.1;
  Parameters::get(keyprefix + "->timeStep", tau);
  TEST_EXIT(tau > 0.0)("Time step size must be positiv. Only Chuck Norris can solve PDEs correctly backward in time.")
  
  step = 0;
  t = t0;
}

void DecProblemInstat::closeTimestep() {
  updateUhOlds(); 
  t += tau;
  step++;
}

void DecProblemInstat::oneIteration() {
  statProb->assembleSystem();
  statProb->solveSystem();
  closeTimestep();
}

void DecProblemInstat::solve() {
  FUNCNAME("DecProblemInstat::solve");
  
  if (initProb) {
    ERROR_EXIT("solving init problem is not implemented") //TODO: implement
  }

  while (t < t1) oneIteration();
}

void DecProblemInstat::updateUhOlds() {
  FUNCNAME("DecProblemInstat::updateUhOlds");
  for (int i = 0; i < statProb->nComponents; i++) {
    switch(statProb->spaceTypes[i]) {
      case EDGESPACE:
          updateUhOlds_EdgeOperators((statProb->vectorOperators)[i], i);
          break;
      default:
        ERROR_EXIT("Das haette nicht passieren duerfen!");
    }
  }
}

void DecProblemInstat::updateUhOlds_EdgeOperators(list<DecOperator*> &ops, int i) {
  for (list<DecOperator*>::const_iterator opIter = ops.begin(); opIter != ops.end(); ++opIter) {
    EdgeOperator *eop = dynamic_cast<EdgeOperator*>(*opIter);
    eop->setUhOld(&(statProb->getSolution(i)));
  }
}
