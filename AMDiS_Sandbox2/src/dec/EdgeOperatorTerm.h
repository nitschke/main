#ifndef EDGEOPERATORTERM_H
#define EDGEOPERATORTERM_H

#include "Dec_fwd.h"
#include "DecOperatorTerm.h"

using namespace std;
namespace AMDiS { namespace dec {


typedef map<DegreeOfFreedom, double> edgeRowValMapper; //edge -> value

class EdgeOperatorTerm : public DecOperatorTerm {
public:
  EdgeOperatorTerm(SpaceType colType = UNDEFINEDSPACE) : DecOperatorTerm(EDGESPACE, colType) {}

  virtual edgeRowValMapper evalRow(const EdgeElement &eel, double factor) {};
};


// < alpha, edge >
class  EdgeVecAtEdges: public EdgeOperatorTerm {
public:
  EdgeVecAtEdges(DofEdgeVector *edgeVector, 
                 AbstractFunction<double,double> *function = NULL, 
                 double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), evec(edgeVector), fac(f), func(function) {};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  DofEdgeVector *evec;
  double fac;
  AbstractFunction<double,double> *func;
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

// < f(||alpha||), edge >; averaging over all linear local spaces Span{edge, adjacent-edge} 
class  NormSquaredEdgeVecAtEdges: public EdgeOperatorTerm {
public:
  NormSquaredEdgeVecAtEdges(DofEdgeVector *edgeVector, 
                            AbstractFunction<double,double> *function = NULL, 
                            double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), fac(f), evec(edgeVector), func(function) {};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
  DofEdgeVector *evec;
  AbstractFunction<double,double> *func;
};



}}
#endif
