#ifndef EDGEOPERATORTERM_H
#define EDGEOPERATORTERM_H

#include "Dec_fwd.h"
#include "DecOperatorTerm.h"

using namespace std;
namespace AMDiS { namespace dec {


typedef map<DegreeOfFreedom, double> edgeRowValMapper;

class EdgeOperatorTerm : public DecOperatorTerm {
public:
  EdgeOperatorTerm(SpaceType colType = UNDEFINEDSPACE) : DecOperatorTerm(EDGESPACE, colType) {}

  virtual edgeRowValMapper evalRow(const EdgeElement &eel, double factor) {};
};


// < alpha, edge >
class  EdgeVecAtEdges: public EdgeOperatorTerm {
public:
  EdgeVecAtEdges(DofEdgeVector *edgeVector, double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), evec(edgeVector), fac(f) {};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  DofEdgeVector *evec;
  double fac;
};


// < alpha , edge >
class  IdentityAtEdges: public EdgeOperatorTerm {
public:
  IdentityAtEdges(double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), fac(f) {};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
};


// < Delta_B(alpha) , edge >
class  LaplaceBeltramiAtEdges: public EdgeOperatorTerm {
public:
  LaplaceBeltramiAtEdges(double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), fac(f) {};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
};



// < Delta_CB(alpha) , edge >
class  LaplaceCoBeltramiAtEdges: public EdgeOperatorTerm {
public:
  LaplaceCoBeltramiAtEdges(double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), fac(f) {};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
};

// < ||alpha||, edge >; averaging over all linear local spaces Span{edge, adjacent-edge} 
class  NormOfEdgeVecAtEdges: public EdgeOperatorTerm {
public:
  NormOfEdgeVecAtEdges(DofEdgeVector *edgeVector, double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), fac(f), evec(edgeVector) {};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
  DofEdgeVector *evec;
};



}}
#endif
