#ifndef DECOPERATORTERM_h
#define DECOPERATORTERM_H

typedef enum {
  UNDEFINEDSPACE = 0,
  VERTEXSPACE= 1,
  EDGESPACE = 2,
  FACESPACE = 3
} SpaceType;


class DecOperatorTerm {
public:
  
  SpaceType getColType() {return colType;}
  SpaceType getRowType() {return rowType;}
  
protected:
 SpaceType colType = 0; //maps from
 SpaceType rowType = 0; //maps to
}

#endif
