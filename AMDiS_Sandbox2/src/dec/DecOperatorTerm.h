#ifndef DECOPERATORTERM_H
#define DECOPERATORTERM_H

namespace AMDiS { namespace dec {


class DecOperatorTerm {
public:

  DecOperatorTerm(SpaceType rowType_ = UNDEFINEDSPACE, SpaceType colType_ = UNDEFINEDSPACE)
    : colType(colType_),
      rowType(rowType_) {}
  
  SpaceType getColType() {return colType;}
  SpaceType getRowType() {return rowType;}

  std::string getName() {return name;}
  

protected:
 SpaceType colType; //maps to SpaceType
 SpaceType rowType; //maps from SpaceType

 std::string name;
};

}}
#endif
