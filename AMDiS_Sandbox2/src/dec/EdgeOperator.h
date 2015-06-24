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
    rProbNum = -1;
  }

  bool isUhOldSet() {
    return uhold;
  }

  //TODO: revise the uhOld-concept!
  void setUhOld(const DofEdgeVector &oldSolution, short rowProblemNumber) {
    uhold = new DofEdgeVector(oldSolution);
    rProbNum = rowProblemNumber;
  }

  void addTerm(EdgeOperatorTerm *term) {
    if (!colType) colType = term->getColType();
    if (colType == term->getColType()) {
      opTs.push_back(term);
    } else {
      ERROR_EXIT("EdgeOperatorTerm--EdgeOperator col space type mismatch!");
    }
  }

  list< EdgeOperatorTerm* >::const_iterator begin() {return opTs.begin();}
  list< EdgeOperatorTerm* >::const_iterator end() {return opTs.end();}

  ~EdgeOperator() {
    for (list< EdgeOperatorTerm* >::iterator opIter = opTs.begin();
         opIter != opTs.end(); ++opIter) {
      delete (*opIter);      
    }
  }

private:

  void setUhOld(DofEdgeVector *oldSolution) {
    if (uhold) delete uhold;
    uhold = oldSolution;
  }

  list< EdgeOperatorTerm* > opTs;

  DofEdgeVector *uhold;
  short rProbNum;

  friend class DecProblemStat;
  friend class DecProblemInstat;
};

}}
#endif
