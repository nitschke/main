#ifndef DECOPERATOR_H
#define DECOPERATOR_H

#include "DecOperatorTerm.h"

class DecOperator {
public:
  
  SpaceType getColType() {return colType;}
  SpaceType getRowType() {return rowType;}

  void setFactorRef(double *factorRef) {factor = factorRef;}

  double getFactor() {
    return (!factor) ? 1.0 : *factor; 
  }
  
protected:
 SpaceType colType = 0;
 SpaceType rowType = 0;

 double *factor = NULL;
}

#endif
