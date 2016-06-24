#ifndef EDGEOPERATORTERM_H
#define EDGEOPERATORTERM_H

#include "Dec_fwd.h"
#include "DecOperatorTerm.h"

using namespace std;
namespace AMDiS { namespace dec {


typedef map<DegreeOfFreedom, double> edgeRowValMapper; //edge |-> value

class EdgeOperatorTerm : public DecOperatorTerm {
public:
  EdgeOperatorTerm(SpaceType colType = UNDEFINEDSPACE) : DecOperatorTerm(EDGESPACE, colType) {}

  virtual edgeRowValMapper evalRow(const EdgeElement &eel, double factor) {};

};


// < f(<beta,edge>) * alpha, edge >
class  EdgeVecAtEdges: public EdgeOperatorTerm {
public:
  EdgeVecAtEdges(DofEdgeVector *edgeVector, 
                 AbstractFunction<double,double> *function = NULL, 
                 double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), evec(edgeVector), fac(f), func(function) {name = "EdgeVecAtEdges";};

  EdgeVecAtEdges(DofEdgeVector *edgeVector,  double f)
      : EdgeOperatorTerm(EDGESPACE), evec(edgeVector), fac(f), func(NULL) {name = "EdgeVecAtEdges";};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  DofEdgeVector *evec;
  double fac;
  AbstractFunction<double,double> *func;
};

// < f(edge) * alpha, edge >
class  EdgeFunAtEdges: public EdgeOperatorTerm {
public:
  EdgeFunAtEdges(AbstractFunction<double, EdgeElement> *function, 
                 double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), fac(f), func(function) {name = "EdgeFunAtEdges";};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
  AbstractFunction<double,EdgeElement> *func;
};

// < f(<beta1,edge>, <beta2,edge>, edge) * alpha, edge >
class  EdgeVec2AndEdgeAtEdges: public EdgeOperatorTerm {
public:
  EdgeVec2AndEdgeAtEdges(DofEdgeVector *edgeVector1, DofEdgeVector *edgeVector2,
                 TertiaryAbstractFunction<double,double,double, EdgeElement> *function, 
                 double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), evec1(edgeVector1), evec2(edgeVector2), fac(f), func(function) {name = "EdgeVec2AndEdgeAtEdges";};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  DofEdgeVector *evec1;
  DofEdgeVector *evec2;
  double fac;
  TertiaryAbstractFunction<double,double,double,EdgeElement> *func;
};


// < alpha , edge >
class  IdentityAtEdges: public EdgeOperatorTerm {
public:
  IdentityAtEdges(double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), fac(f) {name = "IdentityAtEdges";};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
};


// < Delta_B(alpha) , edge >
class  LaplaceBeltramiAtEdges: public EdgeOperatorTerm {
public:
  LaplaceBeltramiAtEdges(double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), fac(f) {name = "LaplaceBeltramiAtEdges";};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
};



// < Delta_CB(alpha) , edge >
class  LaplaceCoBeltramiAtEdges: public EdgeOperatorTerm {
public:
  LaplaceCoBeltramiAtEdges(double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), fac(f) {name = "LaplaceCoBeltramiAtEdges";};
  
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

// < df , edge >
class  ExteriorDerivativeAtEdges: public EdgeOperatorTerm {
public:
  ExteriorDerivativeAtEdges(double f = 1.0) 
      : EdgeOperatorTerm(VERTEXSPACE), fac(f) {name = "ExteriorDerivativeAtEdges";};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
};


}}
#endif
