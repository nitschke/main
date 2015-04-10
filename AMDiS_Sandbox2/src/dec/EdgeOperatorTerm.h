#ifndef EDGEOPERATORTERM_H
#define EDGEOPERATORTERM_H

#include "DecOperatorTerm.h"
#include "EdgeMesh.h"
#include "DofEdgeVector.h"

typedef map<DegreeOfFreedom, double> edgeRowValMapper;

class EdgeOperatorTerm : public DecOperatorTerm {
public:
  EdgeOperatorTerm() : colType(EDGESPACE) {}

  virtual edgeRowValMapper evalRow(EdgeElement *eel, double factor);
};


class  Discrete1FormAtEdge: public EdgeOperatorTerm {
public:
  Discrete1FormAtEdge(DofEdgeVector *edgeVector, double f = 1.0) 
      : EdgeOperatorTerm(), rowType(EDGESPACE), evec(edgeVector), factor(f) {};
  
  edgeRowValMapper evalRow(EdgeElement *eel, double factor);

private:
  DofEdgeVector *evec;
  double factor;
};

#endif
