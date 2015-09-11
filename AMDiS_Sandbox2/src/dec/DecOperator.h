#ifndef DECOPERATOR_H
#define DECOPERATOR_H

#include "Dec_fwd.h"

namespace AMDiS { namespace dec {


class DecOperator {
public:

  DecOperator() 
    : colType(UNDEFINEDSPACE), 
      rowType(UNDEFINEDSPACE) {}
      //{factor = NULL;}

  DecOperator(SpaceType rowType_)
    : colType(UNDEFINEDSPACE), 
      rowType(rowType_) {}
      //{factor = NULL;}
  
  SpaceType getColType() {return colType;}
  SpaceType getRowType() {return rowType;}

  //void setFactorRef(double *factorRef) {factor = factorRef;}

  //double getFactor() {
  //  return (!factor) ? 1.0 : *factor; 
  //}

virtual ~DecOperator() {};
  
protected:
 SpaceType colType;
 SpaceType rowType;

 //double *factor;
};

}}
#endif
