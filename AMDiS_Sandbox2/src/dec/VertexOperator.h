#ifndef VERTEXOPERATOR_H
#define VERTEXOPERATOR_H

#include "Dec_fwd.h"
#include "DecOperator.h"
#include "VertexOperatorTerm.h"

namespace AMDiS { namespace dec {


class VertexOperator : public DecOperator {
public:
  
  VertexOperator() : DecOperator(VERTEXSPACE) {
    uhold = NULL;
  }

  bool isUhOldSet() {
    return uhold;
  }

  void setUhOld(const DOFVector<double> &oldSolution) {
    uhold = new DOFVector<double>(oldSolution);
  }

  void addTerm(VertexOperatorTerm *term) {
    if (!colType) colType = term->getColType();
    if (colType == term->getColType()) {
      opTs.push_back(term);
    } else {
      ERROR_EXIT("VertexOperatorTerm--VertexOperator col space type mismatch!");
    }
  }

  list< VertexOperatorTerm* >::const_iterator begin() {return opTs.begin();}
  list< VertexOperatorTerm* >::const_iterator end() {return opTs.end();}

  ~VertexOperator() {
    for (list< VertexOperatorTerm* >::iterator opIter = opTs.begin();
         opIter != opTs.end(); ++opIter) {
      delete (*opIter);      
    }
  }

private:

  void setUhOld(DOFVector<double> *oldSolution) {
    if (uhold) delete uhold;
    uhold = oldSolution;
  }

  list< VertexOperatorTerm* > opTs;

  DOFVector<double> *uhold;

  friend class DecProblemStat;
  friend class DecProblemInstat;
};

}}
#endif
