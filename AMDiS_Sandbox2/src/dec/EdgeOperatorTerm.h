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


class  Discrete1FormAtEdges: public EdgeOperatorTerm {
public:
  Discrete1FormAtEdges(DofEdgeVector *edgeVector, double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), evec(edgeVector), fac(f) {};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  DofEdgeVector *evec;
  double fac;
};

class  IdentityAtEdges: public EdgeOperatorTerm {
public:
  IdentityAtEdges(double f = 1.0) 
      : EdgeOperatorTerm(EDGESPACE), fac(f) {};
  
  edgeRowValMapper evalRow(const EdgeElement &eel, double factor);

private:
  double fac;
};


}}
#endif
