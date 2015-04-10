#include "DecProblemStat.h"
#include "AMDiS.h"

DecProblemStat::DecProblemStat(ProblemStat *problem) 
        : emesh(NULL),
          sysMat(NULL),
          rhs(NULL),
          ps(problem),
          nComponents(problem->getNumComponents())
          matrixOperators(nComponents, nComponents),
          vectorOperators(nComponents),
          spaceTypes(nComponents) {
  spaceTypes.fill(UNDEFINEDSPACE);
}

void DecProblemStat::addMatrixOperator(DecOperator *op, int row, int col, double *factor) {
  TEST_EXIT(row < nComponents && col < nComponents)
      ("Cannot add matrix operator at position %d/%d. The stationary problem has only %d components!\n",
       row, col, nComponents);

  SpaceType rowType = op->getRowType();
  SpaceType colType = op->getRowType();
  
  // if there is no EdgeMesh but needed
  if ((!eMesh) && (rowType == EDGESPACE || colType == EDGESPACE) emesh = new EdgeMesh(ps->getFeSpace());
  
  if (!spaceTypes[row]) spaceTypes[row] = rowType; 
  if (!spaceTypes[col]) spaceTypes[col] = colType; 

  TEST_EXIT(spaceTypes[row] == rowType)("Wrong SpaceType for the row at adding matrix operator");
  TEST_EXIT(spaceTypes[col] == colType)("Wrong SpaceType for the col at adding matrix operator");

  op->setFactorRef(factor);
  matrixOperators[row][col].push_back(op);
}

void DecProblemStat::addVectorOperator(DecOperator *op, int row, double *factor) {
  TEST_EXIT(row < nComponents)
      ("Cannot add vector operator at position %d. The stationary problem has only %d components!\n",
       row, nComponents);

  SpaceType rowType = op->getRowType();
  
  // if there is no EdgeMesh but needed
  if ((!eMesh) && (rowType == EDGESPACE) emesh = new EdgeMesh(ps->getFeSpace());
  
  if (!spaceTypes[row]) spaceTypes[row] = rowType; 

  TEST_EXIT(spaceTypes[row] == rowType)("Wrong SpaceType for the row at adding vector operator");

  op->setFactorRef(factor);
  vectorOperators[row].push_back(op);
}

void DecProblemStat::assembleSystem() {
  int n = 0;
  vector<int> ns(nComponents);
  for (int i = 0; i < nComponents; ++i) {
     TEST_EXIT(!spaceTypes[i])("Space %d is not set!", i);
     switch (spaceTypes[i]) {
       case EDGESPACE:
                ns[i] = emesh->getNumberOfEdges();
                n += ns[i]; 
                break;
       default: 
                ERROR_EXIT("SpaceType %d is unknown or not implemented", spaceTypes[i]);
     }
  }

  sysMat = new SparseMatrix(n, n);
  rhs = new denseVector(n);

  int ohrow = 0; //overhead
  for (int r = 0; r < nComponents; ++r) {
    int ohcol = 0; //overhead
    for( int c = 0; c < nComponents; ++c) {
      assembleMatrixBlock(matrixOperators[r][c], spaceType[c], ohrow, ohcol);
      ohcol += ns[c];
    }
    //assembleVectorBlock(vectorOperators[r], ???, ohrow);
    ohrow += ns[r];
  }
}

void DecProblemStat::assembleMatrixBlock(list<DecOperator*> &ops, SpaceType colType, int ohrow, int ohcol) {
  list<DecOperator*>::const_iterator opIter.begin(); 
  switch (colType) {
    case EDGESPACE:
          //TODO: nochmal durchdenken, ich bin raus!
          for (; opIter != ops.end(); ++opIter) {
            EdgeOperator *eop = dynamic_cast<EdgeOperator*>(*opIter);
            double factor = eop->getFactor();
            list< EdgeOperatorTerm* >::const_iterator termIter = eop->begin();
            for (; termIter != eop->end(); ++termIter) {
              
            }
          }
          break;
    default:
          ERROR_EXIT("Das haette nicht passieren duerfen!");
  }
}
