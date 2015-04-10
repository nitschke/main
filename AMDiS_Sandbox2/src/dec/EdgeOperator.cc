#include "EdgeOperator.h"

void EdgeOperator::addTerm(EdgeOperatorTerm *term) {
    if (!rowType) rowType = term->getRowType();
    if (rowType == term->getRowType()) {
      opTs.push_back(term);
    } else {
      ERROR_EXIT("EdgeOperatorTerm--EdgeOperator row space type mismatch!");
    }
  }

