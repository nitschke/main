#ifndef DECPROBLEMSTAT_H
#define DECPROBLEMSTAT_H

#include "Dec_fwd.h"
#include "DofEdgeVector.h"
#include "SolverInterface.h"
#include "AnimationWriter.h"

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

  DofEdgeVector getSolution(int i) {
    TEST_EXIT(i < nComponents)("The stationary problem has only %d components!\n", nComponents);
    TEST_EXIT(fullSolution)("there is no solution ");
    TEST_EXIT(spaceTypes[i] == EDGESPACE)("Wrong SpaceType\n");
    
    DofEdgeVector soli(emesh, "Sol_" + boost::lexical_cast<std::string>(i));
    int oh = 0;
    for (int k = 0; k < i; ++k) oh += ns[k];
    mtl::irange range(oh, oh + ns[i]);
    soli.set((*fullSolution)[range]);
    return soli;
  }

  void writeSolution(string nameAddition = "");

  void writeSolution(double time, string nameAddition = "");

private:

  inline void assembleMatrixBlock_EdgeEdge(list<DecOperator*> &ops, int ohrow, int ohcol);

  inline void assembleVectorBlock_Edge(list<DecOperator*> &ops, int ohrow);

//TODO: DESTructur (emesh, sysmat, rhs)

private:


  int nComponents;

  int n;
  vector<int> ns;

  ProblemStat *ps;

  EdgeMesh *emesh;

  SparseMatrix *sysMat;
  DenseVector *rhs;

  Matrix< list<DecOperator*> > matrixOperators;
  Vector< list<DecOperator*> > vectorOperators;
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
