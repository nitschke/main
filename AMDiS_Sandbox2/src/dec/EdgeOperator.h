#ifndef EDGEOPERATOR_H
#define EDGEOPERATOR_H

#include "Dec_fwd.h"
#include "DecOperator.h"
#include "EdgeOperatorTerm.h"
#include "DofEdgeVector.h"

namespace AMDiS { namespace dec {


class EdgeOperator : public DecOperator {
public:
  
  EdgeOperator() : DecOperator(EDGESPACE) {
    uhold = NULL;
  }

  // work on uhold, so their will change in time
  void setUhOld(DofEdgeVector *oldSolution) {
    //if (uhold) delete uhold;
    uhold = oldSolution;
  }

  bool isUhOldSet() {
    return uhold;
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

  DofEdgeVector *uhold;

  friend class DecProblemStat;
};

}}
#endif
