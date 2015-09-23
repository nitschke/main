#include "DecProblemInstat.h"
#include "EdgeOperator.h"
#include "VertexOperator.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

DecProblemInstat::DecProblemInstat(DecProblemStat *probStat, DecProblemStat *initStatProb) 
        : statProb(probStat),  initProb(initStatProb) {
  FUNCNAME("DecProblemInstat::DecProblemInstat");
  string keyprefix = const_cast<ProblemStat*>(statProb->getProblemStat())->getName();
  
  t0 = 0.0;
  Parameters::get(keyprefix + "->adapt->start time", t0);

  t1 = 1.0;
  Parameters::get(keyprefix + "->adapt->end time", t1);
  TEST_EXIT(t0 <= t1)("Start time must be lower then end time.");

  tau = 0.1;
  Parameters::get(keyprefix + "->adapt->timestep", tau);
  TEST_EXIT(tau > 0.0)("Time step size must be positiv. Only Chuck Norris can solve PDEs correctly backward in time.");
  
  step = 0;
  t = t0 + tau;
  inv_tau = 1. / tau;

  writeSolutions = true;
  Parameters::get(keyprefix + "->output->write every i-th timestep", writeEveryithTimestep);
}

void DecProblemInstat::initTimestep() {
  FUNCNAME("DecProblemInstat::initTimestep");
  MSG("*** Start time iteration step %d at t = %g ***\n", step, t);
  timer.reset();
}

void DecProblemInstat::closeTimestep() {
  FUNCNAME("DecProblemInstat::closeTimestep");
  updateUhOlds(); 
  MSG("*** End time iteration step, needed %.5f seconds ***\n\n", timer.elapsed());

  if (writeSolutions && step%writeEveryithTimestep == 0) {
    //statProb->writeSolution("." + boost::lexical_cast<std::string>(step));
    statProb->writeSolution(t);
  }
  t += tau; // move to init timestep?
  step++;
}

void DecProblemInstat::oneIteration() {
  initTimestep();
  statProb->assembleSystem();
  //cout << statProb->getSysMat() << endl;
  //cout << statProb->getRhs() << endl;
  statProb->solve();
  closeTimestep();
}

void DecProblemInstat::solve() {
  FUNCNAME("DecProblemInstat::solve");
  
  if (initProb) {
    ERROR_EXIT("solving init problem is not implemented"); //TODO: implement
  }

  while (t <= t1) oneIteration();
}

//TODO: implement general uhOld-concept
void DecProblemInstat::updateUhOlds() {
  FUNCNAME("DecProblemInstat::updateUhOlds");
  for (int i = 0; i < statProb->nComponents; i++) {
    switch(statProb->spaceTypes[i]) {
      case VERTEXSPACE:
          updateUhOlds_VertexOperators((statProb->vectorOperators)[i], i);
          break;
      case EDGESPACE:
          updateUhOlds_EdgeOperators((statProb->vectorOperators)[i], i);
          break;
      default:
        ERROR_EXIT("Das haette nicht passieren duerfen!");
    }
  }
}

void DecProblemInstat::updateUhOlds_EdgeOperators(list<pair<DecOperator*, double*> > &ops, int i) {
  FUNCNAME("DecProblemInstat::updateUhOlds_EdgeOperators");
  for (list<pair<DecOperator*, double*> >::const_iterator opIter = ops.begin(); opIter != ops.end(); ++opIter) {
    EdgeOperator *eop = dynamic_cast<EdgeOperator*>(opIter->first);
    if (eop->isUhOldSet()) {
      short rProbNum = i;
      if (eop->rProbNum > -1) {
        rProbNum = eop->rProbNum;
        TEST_EXIT(statProb->spaceTypes[i] == statProb->spaceTypes[rProbNum])("SpaceType mismatch");
      }
      DofEdgeVector *soli = new DofEdgeVector(statProb->getSolution(rProbNum));
      eop->setUhOld(soli);
    }
  }
}

void DecProblemInstat::updateUhOlds_VertexOperators(list<pair<DecOperator*, double*> > &ops, int i) {
  for (list<pair<DecOperator*, double*> >::const_iterator opIter = ops.begin(); opIter != ops.end(); ++opIter) {
    VertexOperator *vop = dynamic_cast<VertexOperator*>(opIter->first);
    if (vop->isUhOldSet()) {
      DOFVector<double> *soli = new DOFVector<double>(statProb->getVertexSolution(i));
      vop->setUhOld(soli);
    }
  }
}
