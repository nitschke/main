#include "DecProblemStat.h"
#include "DecOperator.h"
#include "EdgeOperator.h"
#include "VertexOperator.h"
#include "SolverInterface.h"
#include "AnimationWriter.h"
#include "DofVertexVector.h"
#include "DofEdgeVector.h"


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
          solver(ps->getName()),
          animWriterFlat(NULL),
          animWriterSharp(NULL)
{
  spaceTypes.set(UNDEFINEDSPACE);

  writeSharps = false;
  Parameters::get(ps->getName() + "->output->edgeForms sharp", writeSharps);

  writeFlats = false;
  Parameters::get(ps->getName() + "->output->edgeForms flat", writeFlats);

  writeAnimation = false;
  Parameters::get(ps->getName() + "->output->ParaView animation", writeAnimation);
  //TODO: only EDGESPACEs -> temp.fix: VertexSolutions in FlatWriter
  if (writeAnimation) {
    string basename = "output/" + ps->getName();
    Parameters::get(ps->getName() + "->output->filename", basename);
    if (writeSharps) {
      animWriterSharp = new Vector<AnimationWriter*>(nComponents);
      for (int i = 0; i < nComponents; ++i) {
        (*animWriterSharp)[i] = new AnimationWriter(basename + boost::lexical_cast<std::string>(i) + "Sharp.pvd");
      }
    }
    //if (writeFlats) {
      animWriterFlat = new Vector<AnimationWriter*>(nComponents);
      for (int i = 0; i < nComponents; ++i) {
        (*animWriterFlat)[i] = new AnimationWriter(basename + boost::lexical_cast<std::string>(i) + ".pvd");
      }
    //}
  }
}

void DecProblemStat::addMatrixOperator(DecOperator *op, int row, int col, double *factor) {
  TEST_EXIT(row < nComponents && col < nComponents)
      ("Cannot add matrix operator at position %d/%d. The stationary problem has only %d components!\n",
       row, col, nComponents);
  TEST_EXIT(emesh)("No EdgeMesh!");

  SpaceType rowType = op->getRowType();
  SpaceType colType = op->getColType();
  
  if (!spaceTypes[row]) spaceTypes[row] = rowType; 
  if (!spaceTypes[col]) spaceTypes[col] = colType; 

  TEST_EXIT(spaceTypes[row] == rowType)("Wrong SpaceType for the row at adding matrix operator in (%d,%d) \n", row, col);
  TEST_EXIT(spaceTypes[col] == colType)("Wrong SpaceType for the col at adding matrix operator in (%d,%d) \n", row, col);

  matrixOperators[row][col].push_back(make_pair(op,factor));
}

void DecProblemStat::addVectorOperator(DecOperator *op, int row, double *factor) {
  TEST_EXIT(row < nComponents)
      ("Cannot add vector operator at position %d. The stationary problem has only %d components!\n",
       row, nComponents);
  TEST_EXIT(emesh)("No EdgeMesh!");

  SpaceType rowType = op->getRowType();
  
  if (!spaceTypes[row]) spaceTypes[row] = rowType; 

  TEST_EXIT(spaceTypes[row] == rowType)("Wrong SpaceType for the row at adding vector operator in (%d) \n", row);

  vectorOperators[row].push_back(make_pair(op,factor));
}

void DecProblemStat::assembleSystem() {
  FUNCNAME("DecProblemStat::assembleSystem()");
  MSG("Assemble System ... \n");
  Timer t;

  //TODO: dofCompress here, if the mesh change (-> change datastruct dofEdge in EdgeElement too)

  if ( n < 0 ) {   // not precalc 
    n = 0;
    for (int i = 0; i < nComponents; ++i) {
       TEST_EXIT(spaceTypes[i])("Space %d is not set!", i);
       switch (spaceTypes[i]) {
         case VERTEXSPACE:
                  ns[i] = emesh->getFeSpace()->getMesh()->getNumberOfVertices();
                  break;
         case EDGESPACE:
                  ns[i] = emesh->getNumberOfEdges();
                  break;
         default: 
                  ERROR_EXIT("SpaceType %d is unknown or not implemented", spaceTypes[i]);
       }
       n += ns[i];
    }
  }

  if (!sysMat) sysMat = new SparseMatrix(n, n);
  if (!rhs) rhs = new DenseVector(n);

  *rhs = 0.0;
  sysMat->make_empty();

  int ohrow = 0; //overhead
  for (int r = 0; r < nComponents; ++r) {
    int ohcol = 0; //overhead
    for( int c = 0; c < nComponents; ++c) {
      //TODO: only the rowspacetype looks that is neccesarry for the choice of assembling, but proof this in detail
      switch(spaceTypes[r] + 3*spaceTypes[c]) {
        case 4: //VERTEXSPACE x VERTEXSPACE
            assembleMatrixBlock_VertexVertex(matrixOperators[r][c], ohrow, ohcol);
            break;
        case 8: //EDGESPACE x EDGESPACE
            //cout << "\nAssemble Matrixblock (" << r << "," << c << "):" <<  endl;
            assembleMatrixBlock_EdgeEdge(matrixOperators[r][c], ohrow, ohcol); 
            break;
        case 7: //VERTEXSPACE x EDGESPACE
            assembleMatrixBlock_VertexVertex(matrixOperators[r][c], ohrow, ohcol);
            break;
        case 5: //EDGESPACE x VERTEXSPACE
            //cout << "\nAssemble Matrixblock (" << r << "," << c << "):" <<  endl;
            assembleMatrixBlock_EdgeEdge(matrixOperators[r][c], ohrow, ohcol); 
            break;
        default:
          TEST_EXIT(matrixOperators[r][c].empty())
                   ("Unknown SpaceType combination; RowSpace[%d] = %d; ColSpace[%d] = %d \n", r, spaceTypes[r], c, spaceTypes[c]);
      }
      ohcol += ns[c];
    }
    switch(spaceTypes[r]) {
      case VERTEXSPACE:
          assembleVectorBlock_Vertex(vectorOperators[r], ohrow);
          break;
      case EDGESPACE:
          //cout << "\nAssemble Vectorblock (" << r << "):" <<  endl;
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

  for (list< DofCondition >::const_iterator dofCondIter = oneDofConditions.begin(); dofCondIter != oneDofConditions.end(); ++dofCondIter) {
    int row = dofCondIter->rowDof;
    int col = dofCondIter->colDof;
    for (int k = 0; k < dofCondIter->rowComp; ++k) row += ns[k];
    for (int k = 0; k < dofCondIter->colComp; ++k) col += ns[k];
    assembleOneDofCondition(row, col, dofCondIter->val);
  }
}

inline void DecProblemStat::assembleMatrixBlock_EdgeEdge(list<pair<DecOperator*, double*> > &ops, int ohrow, int ohcol) {
  typedef typename mtl::Collection<SparseMatrix>::value_type vtype;
  mtl::matrix::inserter<SparseMatrix, update_plus<vtype> > insSysMat(*sysMat);
  list<pair<DecOperator*,double*> >::const_iterator opIter = ops.begin(); 
  for (; opIter != ops.end(); ++opIter) {
    EdgeOperator *eop = dynamic_cast<EdgeOperator*>(opIter->first);
    double* fref =  opIter->second;
    double factor = (fref) ? (*fref) : 1.0;
    list< EdgeOperatorTerm* >::const_iterator termIter = eop->begin();
    for (; termIter != eop->end(); ++termIter) {
      //TODO: implement edgeMesh::iterator
      Timer tot;
      vector<EdgeElement>::const_iterator edgeIter = emesh->getEdges()->begin();
      for (int r = ohrow; edgeIter != emesh->getEdges()->end(); ++r, ++edgeIter) {
        edgeRowValMapper mapper = (*termIter)->evalRow(*edgeIter, factor);
        edgeRowValMapper::iterator mapIter = mapper.begin();
        for (; mapIter != mapper.end(); ++mapIter) {
          int c = ohcol + mapIter->first;
          insSysMat[r][c] << mapIter->second;
        }
      }
      //cout << "   " << (*termIter)->getName() << ": " << tot.elapsed() << " sec" << endl;
    }
  }
}

//TODO: If there's a vertexmesh (with connected edge relations) implemented, then use vertex iterations 
//      and change also vertex operators in the same manner like the edge operators
inline void DecProblemStat::assembleMatrixBlock_VertexVertex(list<pair<DecOperator*,double*> > &ops, int ohrow, int ohcol) {
  typedef typename mtl::Collection<SparseMatrix>::value_type vtype;
  mtl::matrix::inserter<SparseMatrix, update_plus<vtype> > insSysMat(*sysMat);
  list<pair<DecOperator*,double*> >::const_iterator opIter = ops.begin(); 
  for (; opIter != ops.end(); ++opIter) {
    VertexOperator *vop = dynamic_cast<VertexOperator*>(opIter->first);
    double* fref =  opIter->second;
    double factor = (fref) ? (*fref) : 1.0;
    list< VertexOperatorTerm* >::const_iterator termIter = vop->begin();
    for (; termIter != vop->end(); ++termIter) {
      DOFVector<bool> visited(emesh->getFeSpace(), "visited");
      visited = false;
      vector<EdgeElement>::const_iterator edgeIter = emesh->getEdges()->begin();
      for (; edgeIter != emesh->getEdges()->end(); ++edgeIter) {
        // first vertex
        DegreeOfFreedom dof = edgeIter->dofEdge.first;
        if (!visited[dof]) {
          vertexRowValMapper mapper = (*termIter)->evalRow(*edgeIter, FIRSTVERTEX, factor);
          vertexRowValMapper::iterator mapIter = mapper.begin();
          for (; mapIter != mapper.end(); ++mapIter) {
            insSysMat[ohrow + dof][ohcol +  mapIter->first] << mapIter->second;
          }
          visited[dof] = true;
        }
        // second vertex
        dof = edgeIter->dofEdge.second;
        if (!visited[dof]) {
          vertexRowValMapper mapper = (*termIter)->evalRow(*edgeIter, SECONDVERTEX, factor);
          vertexRowValMapper::iterator mapIter = mapper.begin();
          for (; mapIter != mapper.end(); ++mapIter) {
            insSysMat[ohrow + dof][ohcol + mapIter->first] << mapIter->second;
          }
          visited[dof] = true;
        }
      }
    }
  }
}

//act on uhold if is set in operator
//TODO: see assembleMatrixBlock_VertexVertex
inline void DecProblemStat::assembleVectorBlock_Vertex(list<pair<DecOperator*, double*> > &ops, int ohrow) {
  list<pair<DecOperator*, double*> >::const_iterator opIter = ops.begin(); 
  for (; opIter != ops.end(); ++opIter) {
    VertexOperator *vop = dynamic_cast<VertexOperator*>(opIter->first);
    double* fref =  opIter->second;
    double factor = (fref) ? (*fref) : 1.0;
    list< VertexOperatorTerm* >::const_iterator termIter = vop->begin();
    for (; termIter != vop->end(); ++termIter) {
      DOFVector<bool> visited(emesh->getFeSpace(), "visited");
      visited = false;
      vector<EdgeElement>::const_iterator edgeIter = emesh->getEdges()->begin();
      for (; edgeIter != emesh->getEdges()->end(); ++edgeIter) {
        // first vertex
        DegreeOfFreedom dof = edgeIter->dofEdge.first;
        if (!visited[dof]) {
          double val = 0.0;
          vertexRowValMapper mapper = (*termIter)->evalRow(*edgeIter, FIRSTVERTEX, factor);
          vertexRowValMapper::iterator mapIter = mapper.begin();
          for (; mapIter != mapper.end(); ++mapIter) {
            if (vop->uholdAtVs) {
              val += mapIter->second * ((*(vop->uholdAtVs))[mapIter->first]);
            } else if (vop->uholdAtEs) {
              val += mapIter->second * ((*(vop->uholdAtEs))[mapIter->first]);
            } else {
              val += mapIter->second;
            }
          }
          (*rhs)[ohrow + dof] += val;
          visited[dof] = true;
        }
        // second vertex
        dof = edgeIter->dofEdge.second;
        if (!visited[dof]) {
          double val = 0.0;
          vertexRowValMapper mapper = (*termIter)->evalRow(*edgeIter, SECONDVERTEX, factor);
          vertexRowValMapper::iterator mapIter = mapper.begin();
          for (; mapIter != mapper.end(); ++mapIter) {
            if (vop->uholdAtVs) {
              val += mapIter->second * ((*(vop->uholdAtVs))[mapIter->first]);
            } else if (vop->uholdAtEs) {
              val += mapIter->second * ((*(vop->uholdAtEs))[mapIter->first]);
            } else {
              val += mapIter->second;
            }
          }
          (*rhs)[ohrow + dof] += val;
          visited[dof] = true;
        }
      }
    }
  }
}


//act on uhold if is set in operator
inline void DecProblemStat::assembleVectorBlock_Edge(list<pair<DecOperator*, double*> > &ops, int ohrow) {
  list<pair<DecOperator*,double*> >::const_iterator opIter = ops.begin(); 
  for (; opIter != ops.end(); ++opIter) {
    EdgeOperator *eop = dynamic_cast<EdgeOperator*>(opIter->first);
    double* fref =  opIter->second;
    double factor = (fref) ? (*fref) : 1.0;
    list< EdgeOperatorTerm* >::const_iterator termIter = eop->begin();
    for (; termIter != eop->end(); ++termIter) {
      //TODO: implement edgeMesh::iterator
      Timer tot;
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
      //cout << "   " << (*termIter)->getName() << ": " << tot.elapsed() << " sec" << endl;
    }
  }
}

void DecProblemStat::setValAtDof(int rowComp, int colComp, DegreeOfFreedom rowDof, DegreeOfFreedom colDof, double val) {
  FUNCNAME("DecProblemStat::setValAtDof(int comp, DegreeOfFreedom dof, double val)");

  DofCondition dofCond = {rowComp, colComp, rowDof, colDof, val};
  oneDofConditions.push_back(dofCond);
}

void DecProblemStat::assembleOneDofCondition(int row, int col, double val) {
  using namespace mtl; 
  FUNCNAME("DecProblemStat::assembleOneDofCondition(const pair<unsigned int, double> &dofCond)");

  { 
    matrix::inserter<SparseMatrix> ins(*sysMat);
    ins.make_empty(row);
    ins[row][col] << 1.0;
  }
  sysMat->crop();
  
  (*rhs)[row] = val;

  int nnz = sysMat->nnz();
  int numRows = sysMat->num_rows();
  int numCols = sysMat->num_cols();
  MSG("New Fill-in of assembled (%ix%i)-Matrix: %i (approx. %.1f per row)\n", numRows, numCols, nnz, ((double)(nnz)/numRows));
}


void DecProblemStat::solve() {
  using namespace mtl;
  using namespace itl;
  FUNCNAME("DecProblemStat::solve()");

  TEST_EXIT(n > 0)("System is not assembled\n");

  if (!fullSolution) fullSolution = new DenseVector(n, 0.0);

  MSG("Solve system ... (with %s)\n", solver.getSolverName().c_str());
  Timer t;

  solver.init(sysMat);
  solver.solve(*rhs, *fullSolution);

  MSG("solving needed %.5f seconds\n", t.elapsed());
}

void DecProblemStat::writeSolution(string nameAddition) {
  string basename = "output/" + ps->getName();
  Parameters::get(ps->getName() + "->output->filename", basename);

  for (int i = 0; i < nComponents; ++i) {
    string bni = basename + boost::lexical_cast<std::string>(i);
    switch(spaceTypes[i]) {
      case VERTEXSPACE: {
            DOFVector<double> soli = getVertexSolution(i);
            io::VtkVectorWriter::writeFile(soli, bni + nameAddition + ".vtu");
          }
          break;
      case EDGESPACE: {
            DofEdgeVector soli = getSolution(i);
            if (writeFlats) soli.writeFile(bni + nameAddition + ".vtu");
            if (writeSharps) {
              string fn = bni + "Sharp" + nameAddition + ".vtu";
              soli.writeSharpFile(fn, ps);
            }
          }
          break;
      default:
        ERROR_EXIT("Das haette nicht passieren duerfen!");
    }
  }
}

void DecProblemStat::writeSolution(double time, string nameAddition) {
  string basename = "output/" + ps->getName();
  Parameters::get(ps->getName() + "->output->filename", basename);
  
  int prec = 3;
  Parameters::get("sphere->output->index precision", prec);
  TEST_EXIT(prec > 0)("Precision must be greater than 0");
  ostringstream timeoss;
  timeoss << setprecision(prec) << fixed << time;
  nameAddition += "." + timeoss.str();

  for (int i = 0; i < nComponents; ++i) {
    string bni = basename + boost::lexical_cast<std::string>(i);
    switch(spaceTypes[i]) {
      case VERTEXSPACE: {
            DOFVector<double> soli = getVertexSolution(i);
            string fn = bni + nameAddition + ".vtu";
            io::VtkVectorWriter::writeFile(soli, fn);
            (*animWriterFlat)[i]->updateAnimationFile(time, fn);
          }
          break;
      case EDGESPACE: {
            DofEdgeVector soli = getSolution(i);
            if (writeFlats) {
              string fn = bni + nameAddition + ".vtu";
              soli.writeFile(fn);
              (*animWriterFlat)[i]->updateAnimationFile(time, fn);
            }
            if (writeSharps) {
              string fn = bni + "Sharp" + nameAddition + ".vtu";
              soli.writeSharpFile(fn, ps);
              (*animWriterSharp)[i]->updateAnimationFile(time, fn);
            }
          }
          break;
      default:
        ERROR_EXIT("Das haette nicht passieren duerfen!");
    }
  }
}

//void DecProblemStat::solveDeprecated() {
//  using namespace mtl;
//  using namespace itl;
//  FUNCNAME("DecProblemStat::solveDeprecated()");
//
//  string solverName = "cgs";
//  Parameters::get(ps->getName() + "->solver", solverName);
//  double tol = 1.e-6;
//  Parameters::get(ps->getName() + "->solver->tolerance", tol);
//  int maxIter = 1000;
//  Parameters::get(ps->getName() + "->solver->max iteration", maxIter);
//
//  MSG("Solve system ... (with %s)\n", solverName.c_str());
//  Timer t;
//
//  if (!fullSolution) fullSolution = new DenseVector(n);
//  pc::identity<SparseMatrix> L(*sysMat);
//  pc::identity<SparseMatrix> R(*sysMat);
//  cyclic_iteration<double> iter(*rhs, maxIter, tol);
//  while (true) {
//    if (solverName == "cgs") {cgs(*sysMat, *fullSolution, *rhs, L, iter); break;}
//    if (solverName == "umfpack") {umfpack_solve(*sysMat, *fullSolution, *rhs); break;}
//    if (solverName == "bicgstab2") {bicgstab_ell(*sysMat, *fullSolution, *rhs, L, R, iter, 2); break;}
//    if (solverName == "bicgstab_ell") { int ell = 3;
//                          Parameters::get(ps->getName() + "->solver->ell", ell);
//                          bicgstab_ell(*sysMat, *fullSolution, *rhs, L, R, iter, ell);
//                          break;}
//    if (solverName == "tfqmr") {tfqmr(*sysMat, *fullSolution, *rhs, L, R, iter); break;}
//    ERROR_EXIT("Solver %s is not known\n", solverName.c_str());
//  }
//
//  MSG("solving needed %.5f seconds\n", t.elapsed());
//}


DofEdgeVector DecProblemStat::getSolution(int i) {
    FUNCNAME("DecProblemStat::getSolution(int i) aka DecProblemStat::getEdgeSolution(int i)");
    TEST_EXIT(i < nComponents)("The stationary problem has only %d components!\n", nComponents);
    TEST_EXIT(fullSolution && i >= 0)("there is no solution ");
    TEST_EXIT(spaceTypes[i] == EDGESPACE)("Wrong SpaceType\n");
    
    DofEdgeVector soli(emesh, "Sol_" + boost::lexical_cast<std::string>(i));
    int oh = 0;
    for (int k = 0; k < i; ++k) oh += ns[k];
    mtl::irange range(oh, oh + ns[i]);
    soli.set((*fullSolution)[range]);
    return soli;
}

DofVertexVector  DecProblemStat::getVertexSolution(int i) {
    FUNCNAME("DecProblemStat::getVertexSolution(int i)");
    TEST_EXIT(i < nComponents)("The stationary problem has only %d components!\n", nComponents);
    TEST_EXIT(fullSolution && i >= 0)("there is no solution ");
    TEST_EXIT(spaceTypes[i] == VERTEXSPACE)("Wrong SpaceType\n");
    
    DofVertexVector soli(emesh, "Sol_" + boost::lexical_cast<std::string>(i));
    int oh = 0;
    for (int k = 0; k < i; ++k) oh += ns[k];
    for (int k = 0; k < ns[i]; k++) {
      soli[k] = (*fullSolution)[oh + k];
    }
    return soli;
}

