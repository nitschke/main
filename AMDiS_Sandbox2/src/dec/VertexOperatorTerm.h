#ifndef VERTEXOPERATORTERM_H
#define VERTEXOPERATORTERM_H

#include "Dec_fwd.h"
#include "DecOperatorTerm.h"

using namespace std;

namespace AMDiS { namespace dec {

typedef map<DegreeOfFreedom, double> vertexRowValMapper; // vertex -> value

//typedef enum {
//  FIRSTVERTEX = 1,
//  SECONDVERTEX = 2
//} VertexPosition; //EdgeRingIteratorType conform

//TODO: throw exception for *FACE types or make it friends to problemstat
typedef EdgeRingIteratorType VertexPosition;

class VertexOperatorTerm : public DecOperatorTerm {
public:
  VertexOperatorTerm(SpaceType colType = UNDEFINEDSPACE) : DecOperatorTerm(VERTEXSPACE, colType) {}

  virtual vertexRowValMapper evalRow(const EdgeElement &eel, VertexPosition pos, double factor) {};
};



// < LB(f) , vertex >
class  LaplaceBeltramiAtVertices: public VertexOperatorTerm {
public:
  LaplaceBeltramiAtVertices(double f = 1.0) 
      : VertexOperatorTerm(VERTEXSPACE), fac(f) {};
  
  vertexRowValMapper evalRow(const EdgeElement &eel, VertexPosition pos, double factor);

protected:
  double fac;
};

// < g(vec[vertex])*LB(f) , vertex >
class  VertexVecLaplaceBeltramiAtVertices: public LaplaceBeltramiAtVertices {
public:
  VertexVecLaplaceBeltramiAtVertices(DOFVector<double> *dofv, 
                  AbstractFunction<double,double> *function = NULL, double f = 1.0) 
      : LaplaceBeltramiAtVertices(f), dv(dofv), func(function) {};
  
  vertexRowValMapper evalRow(const EdgeElement &eel, VertexPosition pos, double factor);

private:
  DOFVector<double> *dv;
  AbstractFunction<double,double> *func;
};


class IdentityAtVertices: public VertexOperatorTerm {
public:
  IdentityAtVertices(double f = 1.0)
      : VertexOperatorTerm(VERTEXSPACE), fac(f) {};

  vertexRowValMapper evalRow(const EdgeElement &eel, VertexPosition pos, double factor);

private:
  double fac;
};

class VertexVecAtVertices: public VertexOperatorTerm {
public:
  VertexVecAtVertices(DOFVector<double> *dofv, double f = 1.0)
      : VertexOperatorTerm(VERTEXSPACE), fac(f), dv(dofv) {};

  vertexRowValMapper evalRow(const EdgeElement &eel, VertexPosition pos, double factor);

private:
  double fac;
  DOFVector<double> *dv;
};



// < div(alpha) , vertex >
class  DivAtVertices: public VertexOperatorTerm {
public:
  DivAtVertices(double f = 1.0) 
      : VertexOperatorTerm(EDGESPACE), fac(f) {name="DivAtVertices";};
  
  vertexRowValMapper evalRow(const EdgeElement &eel, VertexPosition pos, double factor);

protected:
  double fac;
};

}}

#endif
