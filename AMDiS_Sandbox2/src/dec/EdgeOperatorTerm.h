#ifndef EDGEOPERATORTERM_H
#define EDGEOPERATORTERM_H

#include "Dec_fwd.h"
#include "DecOperatorTerm.h"
#include "VertexOperatorTerm.h"

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
// DEPRECATED
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

// < *df , edge >
class  RotAtEdges: public EdgeOperatorTerm {
public:
  RotAtEdges(double f = 1.0) 
      : EdgeOperatorTerm(VERTEXSPACE), fac(f) {name = "RotAtEdges";};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
};

// < f alpha , edge > 
// with given alpha
// f will be middled on edge
class  AverageVertexAndEdgeVecAtEdges: public EdgeOperatorTerm {
public:
  AverageVertexAndEdgeVecAtEdges(DofEdgeVector *edgeVector, double f = 1.0) 
      : EdgeOperatorTerm(VERTEXSPACE), evec(edgeVector), fac(f) {name = "AverageVertexAndEdgeVecAtEdges";};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
  DofEdgeVector *evec;
};

// < rot(alpha) gamma , e > = <rot(alpha),c(e)><gamma,e>
// with given gamma, rot alpha will be calculated around left and right face of edge e
class  RotAtEdgeCenterAndEdgeVecAtEdges: public EdgeOperatorTerm {
public:
  RotAtEdgeCenterAndEdgeVecAtEdges(DofEdgeVector *edgeVector, double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), evec(edgeVector), fac(f) {name = " RotAtEdgeCenterAndEdgeVecAtEdge";};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
  DofEdgeVector *evec;
};


//TODO: fails (no convergence) in directionalDerivativeTest2, hence do a standalone test,
//      maybe, the considered area must be convex -> way out: also take opposite voronoi areas in calculation
// < div(alpha) gamma , e > = <div(alpha),c(e)><gamma,e>
// with given gamma, div alpha will be calculated in the voronoi arias around first and second vertex of edge e
class  DivAtEdgeCenterAndEdgeVecAtEdges: public EdgeOperatorTerm {
public:
  DivAtEdgeCenterAndEdgeVecAtEdges(DofEdgeVector *edgeVector, double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE),
        evec(edgeVector), fac(f) {
          name = " DivAtEdgeCenterAndEdgeVecAtEdge";
          divOp = new DivAtVertices();
        }
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
  DofEdgeVector *evec;

  DivAtVertices *divOp;
};


// < *alpha , edge >
class  HodgeAtEdges: public EdgeOperatorTerm {
public:
  HodgeAtEdges(double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), fac(f) {name = "LaplaceCoBeltramiAtEdges";};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
};


}}
#endif
