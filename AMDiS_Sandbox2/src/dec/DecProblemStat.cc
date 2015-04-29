#include "DecProblemStat.h"
#include "DecOperator.h"
#include "EdgeOperator.h"


using namespace AMDiS;
using namespace dec;

DecProblemStat::DecProblemStat(ProblemStat *problem, EdgeMesh *edgeMesh) 
        : emesh(edgeMesh),
          sysMat(NULL),
          rhs(NULL),
          ps(problem),
          n(-1),
          nComponents(problem->getNumComponents()),
          ns(nComponents),
          matrixOperators(nComponents, nComponents),
          vectorOperators(nComponents),
          spaceTypes(nComponents),
          fullSolution(NULL),
          solver(ps->getName()){
  spaceTypes.fill(UNDEFINEDSPACE);
}

void DecProblemStat::addMatrixOperator(DecOperator *op, int row, int col, double *factor) {
  TEST_EXIT(row < nComponents && col < nComponents)
      ("Cannot add matrix operator at position %d/%d. The stationary problem has only %d components!\n",
       row, col, nComponents);

  SpaceType rowType = op->getRowType();
  SpaceType colType = op->getRowType();
  
  // if there is no EdgeMesh but needed
  if ((!emesh) && (rowType == EDGESPACE || colType == EDGESPACE)) emesh = new EdgeMesh(ps->getFeSpace());
  
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
  if ((!emesh) && (rowType == EDGESPACE)) emesh = new EdgeMesh(ps->getFeSpace());
  
  if (!spaceTypes[row]) spaceTypes[row] = rowType; 

  TEST_EXIT(spaceTypes[row] == rowType)("Wrong SpaceType for the row at adding vector operator");

  op->setFactorRef(factor);
  vectorOperators[row].push_back(op);
}

void DecProblemStat::assembleSystem() {
  FUNCNAME("DecProblemStat::assembleSystem()");
  MSG("Assemble System ... \n");
  Timer t;

  if ( n < 0 ) {   // not precalc 
    n = 0;
    for (int i = 0; i < nComponents; ++i) {
       TEST_EXIT(spaceTypes[i])("Space %d is not set!", i);
       switch (spaceTypes[i]) {
         case EDGESPACE:
                  ns[i] = emesh->getNumberOfEdges();
                  n += ns[i]; 
                  break;
         default: 
                  ERROR_EXIT("SpaceType %d is unknown or not implemented", spaceTypes[i]);
       }
    }
  }

  if (!sysMat) sysMat = new SparseMatrix(n, n);
  if (!rhs) rhs = new DenseVector(n);

  int ohrow = 0; //overhead
  for (int r = 0; r < nComponents; ++r) {
    int ohcol = 0; //overhead
    for( int c = 0; c < nComponents; ++c) {
      switch(spaceTypes[r] + 3*spaceTypes[c]) {
        case 8: //EDGESPACE x EDGESPACE
            assembleMatrixBlock_EdgeEdge(matrixOperators[r][c], ohrow, ohcol); 
            break;
        default:
          ERROR_EXIT("Das haette nicht passieren duerfen!");
      }
      ohcol += ns[c];
    }
    switch(spaceTypes[r]) {
      case EDGESPACE:
          assembleVectorBlock_Edge(vectorOperators[r], ohrow);
          break;
      default:
        ERROR_EXIT("Das haette nicht passieren duerfen!");
    }
    ohrow += ns[r];
  }

  MSG("done (needed %.5f seconds)\n", t.elapsed());
  int nnz = sysMat->nnz();
  int numRows = sysMat->num_rows();
  int numCols = sysMat->num_cols();
  MSG("Fill-in of assembled (%ix%i)-Matrix: %i (approx. %.1f per row)\n", numRows, numCols, nnz, ((double)(nnz)/numRows)); 
}

inline void DecProblemStat::assembleMatrixBlock_EdgeEdge(list<DecOperator*> &ops, int ohrow, int ohcol) {
  typedef typename mtl::Collection<SparseMatrix>::value_type vtype;
  mtl::matrix::inserter<SparseMatrix, update_plus<vtype> > insSysMat(*sysMat);
  list<DecOperator*>::const_iterator opIter = ops.begin(); 
  for (; opIter != ops.end(); ++opIter) {
    EdgeOperator *eop = dynamic_cast<EdgeOperator*>(*opIter);
    double factor = eop->getFactor();
    list< EdgeOperatorTerm* >::const_iterator termIter = eop->begin();
    for (; termIter != eop->end(); ++termIter) {
      //TODO: implement edgeMesh::iterator
      vector<EdgeElement>::const_iterator edgeIter = emesh->getEdges()->begin();
      for (int r = ohrow; edgeIter != emesh->getEdges()->end(); ++r, ++edgeIter) {
        edgeRowValMapper mapper = (*termIter)->evalRow(*edgeIter, factor);
        edgeRowValMapper::iterator mapIter = mapper.begin();
        for (; mapIter != mapper.end(); ++mapIter) {
          int c = ohcol + mapIter->first;
          insSysMat[r][c] << mapIter->second;
        }
      }
    }
  }
}

//TODO: act on uhold if is set in operator
inline void DecProblemStat::assembleVectorBlock_Edge(list<DecOperator*> &ops, int ohrow) {
  list<DecOperator*>::const_iterator opIter = ops.begin(); 
  for (; opIter != ops.end(); ++opIter) {
    EdgeOperator *eop = dynamic_cast<EdgeOperator*>(*opIter);
    double factor = eop->getFactor();
    list< EdgeOperatorTerm* >::const_iterator termIter = eop->begin();
    for (; termIter != eop->end(); ++termIter) {
      //TODO: implement edgeMesh::iterator
      vector<EdgeElement>::const_iterator edgeIter = emesh->getEdges()->begin();
      for (int r = ohrow; edgeIter != emesh->getEdges()->end(); ++r, ++edgeIter) {
        double val = 0.0;
        edgeRowValMapper mapper = (*termIter)->evalRow(*edgeIter, factor);
        edgeRowValMapper::iterator mapIter = mapper.begin();
        for (; mapIter != mapper.end(); ++mapIter) {
          if (eop->uhold) {
            val += mapIter->second * ((*(eop->uhold))[mapIter->first]);
          } else {
            val += mapIter->second;
          }
        }
        (*rhs)[r] += val;
      }
    }
  }
}


void DecProblemStat::solve() {
  using namespace mtl;
  using namespace itl;
  FUNCNAME("DecProblemStat::solve()");

  if (!fullSolution) fullSolution = new DenseVector(n);

  solver.init(sysMat);
  solver.solve(*rhs, *fullSolution);
}

void DecProblemStat::writeSolution(string nameAddition) {
  string basename = "output/" + ps->getName();
  Parameters::get(ps->getName() + "->output->filename", basename);
  
  bool writeSharps = false;
  Parameters::get(ps->getName() + "->output->edgeForms sharp", writeSharps);

  bool writeFlats = false;
  Parameters::get(ps->getName() + "->output->edgeForms flat", writeFlats);

  for (int i = 0; i < nComponents; ++i) {
    string bni = basename + boost::lexical_cast<std::string>(i);
    switch(spaceTypes[i]) {
      case EDGESPACE: {
            DofEdgeVector soli = getSolution(i);
            if (writeFlats) soli.writeFile(bni + nameAddition + ".vtu");
            if (writeSharps) {
              DOFVector< WorldVector<double> > soliSharp = soli.getSharpFaceAverage();
              io::VtkVectorWriter::writeFile(soliSharp, bni + "Sharp" + nameAddition + ".vtu");
            }
          }
          break;
      default:
        ERROR_EXIT("Das haette nicht passieren duerfen!");
    }
  }
}
