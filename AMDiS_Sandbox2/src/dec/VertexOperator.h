#ifndef VERTEXOPERATOR_H
#define VERTEXOPERATOR_H

#include "Dec_fwd.h"
#include "DecOperator.h"
#include "VertexOperatorTerm.h"

namespace AMDiS { namespace dec {


class VertexOperator : public DecOperator {
public:
  
  VertexOperator() : DecOperator(VERTEXSPACE) {
    uholdAtVs = NULL;
    uholdAtEs = NULL;
    compNum = -1;
  }

  bool isUhOldSet() {
    return uholdAtVs || uholdAtEs;
  }

  void setUhOldAtVertices(const DOFVector<double> &oldSolution, short componentNumber) {
    FUNCNAME("void VertexOperator::setUhOldAtVertices(const DOFVector<double> &oldSolution)");
    TEST_EXIT(colType == VERTEXSPACE)("VertexOperator col space type mismatch: You can only use DofVertexVector for old Solution!");
    uholdAtVs = new DOFVector<double>(oldSolution);
    compNum = componentNumber;
  }

  void setUhOldAtEdges(const DofEdgeVector &oldSolution, short componentNumber) {
    FUNCNAME("void VertexOperator::setUhOldAtEdges(const DOFVector<double> &oldSolution)");
    TEST_EXIT(colType == EDGESPACE)("VertexOperator col space type mismatch: You can only use DofEdgeVector for old Solution!");
    uholdAtEs = new DofEdgeVector(oldSolution);
    compNum = componentNumber;
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

  void setUhOldAtVertices(DOFVector<double> *oldSolution) {
    if (uholdAtVs) delete uholdAtVs;
    uholdAtVs = oldSolution;
  }

  void setUhOldAtEdges(DofEdgeVector *oldSolution) {
    if (uholdAtEs) delete uholdAtEs;
    uholdAtEs = oldSolution;
  }

  list< VertexOperatorTerm* > opTs;

  DOFVector<double> *uholdAtVs;
  DofEdgeVector *uholdAtEs;

  short compNum;

  friend class DecProblemStat;
  friend class DecProblemInstat;
};

}}
#endif
