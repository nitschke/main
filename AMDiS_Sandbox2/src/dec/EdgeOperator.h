#ifndef EDGEOPERATOR_H
#define EDGEOPERATOR_H

#include "EdgeOperator.h"
#include "DecOperator"

class EdgeOperator : public DecOperator {
  
  EdgeOperator() : colType(EDGESPACE) {}

  void addTerm(EdgeOperatorTerm *term);

  list< EdgeOperatorTerm* >::const_iterator begin() {return opTs.cbegin()}
  list< EdgeOperatorTerm* >::const_iterator end() {return opTs.cend()}


private:

  list< EdgeOperatorTerm* > opTs;
};

#endif
