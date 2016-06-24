#ifndef DECPROBLEMSTAT_H
#define DECPROBLEMSTAT_H

#include "Dec_fwd.h"
#include "SolverInterface.h"
#include "DofEdgeVector.h"

using namespace std;
namespace AMDiS { namespace dec {


class DecProblemStat {
public:
  
  DecProblemStat(ProblemStat *problem, EdgeMesh *edgeMesh = NULL);

  const EdgeMesh* getMesh() {
    return emesh;
  }

  const ProblemStat* getProblemStat() {
    return ps;
  }

  string getName() {
    return ps->getName();
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


  const SparseMatrix getSysMat() {
    TEST_EXIT(sysMat)("sysMat is not initialized\n");
    return *sysMat;
  }

  const DenseVector getRhs() {
    TEST_EXIT(rhs)("rhs is not initialized\n");
    return *rhs;
  }

  const DenseVector getFullSolution() {
    TEST_EXIT(fullSolution)("there is no fullSolution\n");
    return *fullSolution;
  }

  void solve();
  void solveDeprecated();

  DofEdgeVector getSolution(int i = 0);

  DofEdgeVector getEdgeSolution(int i) {return getSolution(i);}

  Vector<DofEdgeVector> getSolutions()
  {
    Vector<DofEdgeVector> sols(nComponents);
    for (int i = 0; i < nComponents; ++i) {
      sols[i] = getSolution(i);
    }
    return sols;
  }

  DofVertexVector getVertexSolution(int i);

  void writeSolution(string nameAddition = "");

  void writeSolution(double time, string nameAddition = "");

private:

  inline void assembleMatrixBlock_EdgeEdge(list<pair<DecOperator*,double*> > &ops, int ohrow, int ohcol);
  inline void assembleMatrixBlock_VertexVertex(list<pair<DecOperator*,double*> > &ops, int ohrow, int ohcol);

  inline void assembleVectorBlock_Edge(list<pair<DecOperator*,double*> > &ops, int ohrow);
  inline void assembleVectorBlock_Vertex(list<pair<DecOperator*, double*> > &ops, int ohrow);


//TODO: DESTructur (emesh, sysmat, rhs)

private:

  int nComponents;

  int n;
  vector<int> ns;

  ProblemStat *ps;

  EdgeMesh *emesh;

  SparseMatrix *sysMat;
  DenseVector *rhs;

  Matrix< list<pair<DecOperator*,double*> > > matrixOperators;
  Vector< list<pair<DecOperator*,double*> > > vectorOperators;
  Vector< SpaceType > spaceTypes;

  DenseVector *fullSolution;

  SolverInterface solver;

  bool writeSharps;
  bool writeFlats;

  bool writeAnimation;
  Vector < AnimationWriter* > *animWriterFlat;
  Vector < AnimationWriter* > *animWriterSharp;

  friend class DecProblemInstat; 
};

}}
#endif
