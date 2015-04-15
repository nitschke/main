#include "EdgeOperator.h"

using namespace AMDiS;
using namespace dec;


void EdgeOperator::addTerm(EdgeOperatorTerm *term) {
    if (!colType) colType = term->getColType();
    if (colType == term->getColType()) {
      opTs.push_back(term);
    } else {
      ERROR_EXIT("EdgeOperatorTerm--EdgeOperator col space type mismatch!");
    }
  }

