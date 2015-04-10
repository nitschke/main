#ifndef EDGEPROBLEMSTAT_H
#define EDGEPROBLEMSTAT_H

#include "AMDiS.h"
#include "EdgeMesh.h"
#include "EdgeOperator.h"

typedef mtl::compressed2D<double> SparseMatrix;
typedef mtl::dense_vector<double> DenseVector;

class DecProblemStat {
public:
  
  DecProblemStat(ProblemStat *problem);

  const EdgeMesh* getMesh() {
    return emesh;
  }

  void addMatrixOperator(DecOperator *op, int row, int col, double *factor = NULL);
  void addMatrixOperator(DecOperator &op, int row, int col, double *factor = NULL) {
    addMatrixOperator(&op, row, col, factor);
  }

  void addVectorOperator(DecOperator *op, int row, double *factor = NULL);
  void addVectorOperator(DecOperator &op, int row, double *factor = NULL) {
    addVectorOperator(&op, row, factor);
  }
  
  void assembleSystem();

private:

  void assembleMatrixBlock(list<DecOperator*> &ops, SpaceType colType, int ohrow, int ohcol);

//TODO: DESTructur (emesh, sysmat, rhs)

private:



  ProblemStat *ps;

  EdgeMesh *emesh;

  SparseMatrix *sysMat;
  DenseVector *rhs;

  Matrix< list<DecOperator*> > matrixOperators;
  Vector< list<DecOperator*> > vectorOperators;
  Vector< SpaceType > spaceTypes;

  int nComponents;
};

#endif
