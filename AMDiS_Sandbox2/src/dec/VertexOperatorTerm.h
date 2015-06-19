#ifndef VERTEXOPERATORTERM_H
#define VERTEXOPERATORTERM_H

#include "Dec_fwd.h"
#include "DecOperatorTerm.h"
#include "EdgeMesh.h"

namespace AMDiS { namespace dec {

typedef map<DegreeOfFreedom, double> vertexRowValMapper; // vertex -> value

typedef enum {
  FIRSTVERTEX = 1,
  SECONDVERTEX = 2
} VertexPosition; //EdgeRingIteratorType conform

class VertexOperatorTerm : puplic DecOperatorTerm {
public:
  VertexOperatorTerm(SpaceType colType = UNDEFINEDSPACE) : DecOperatorTerm(VERTEXSPACE, colType) {}

  virtual vertexRowValMapper evalRow(const EdgeElement &eel, VertexPosition pos, double factor) {};
};



// < Delta_B(f) , vertex >
class  LaplaceBeltramiAtVertices: public VertexOperatorTerm {
public:
  LaplaceBeltramiAtVertices(double f = 1.0) 
      : VertexOperatorTerm(VERTEXSPACE), fac(f) {};
  
  vertexRowValMapper evalRow(const EdgeElement &eel, VertexPosition pos, double factor);

private:
  double fac;
};



}}

#endif
