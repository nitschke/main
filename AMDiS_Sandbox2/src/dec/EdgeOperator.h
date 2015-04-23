#ifndef EDGEOPERATOR_H
#define EDGEOPERATOR_H

#include "Dec_fwd.h"
#include "DecOperator.h"
#include "EdgeOperatorTerm.h"

namespace AMDiS { namespace dec {


class EdgeOperator : public DecOperator {
public:
  
  EdgeOperator() : DecOperator(EDGESPACE) {
    uhold = NULL;
  }

  // work on uhold, so their will change in time
  void setUhOld(EdgeVector *oldSolution) {
    if (uhold) delete uhold;
    uhold = oldSolution;
  }

  void addTerm(EdgeOperatorTerm *term);

  list< EdgeOperatorTerm* >::const_iterator begin() {return opTs.begin();}
  list< EdgeOperatorTerm* >::const_iterator end() {return opTs.end();}

  ~EdgeOperator() {
    for (list< EdgeOperatorTerm* >::iterator opIter = opTs.begin();
         opIter != opTs.end(); ++opIter) {
      delete (*opIter);      
    }
  }

private:

  list< EdgeOperatorTerm* > opTs;

  EdgeVector *uhold;

  friend class DecProblemStat;
};

}}
#endif
