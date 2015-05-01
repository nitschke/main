#ifndef SOLVERINTERFACE_H
#define SOLVERINTERFACE_H

#include "Dec_fwd.h"
#include "solver/LinearSolverInterface.h"


namespace AMDiS { namespace dec {

class SolverInterface {
public:

  SolverInterface(std::string name)
    : solver(NULL), runner(NULL)
  {
    FUNCNAME("SolverInterface::SolverInterface");
    std::string initFileStr = name + "->solver";
    std::string solverName = "H2SO4";

    Parameters::get(initFileStr, solverName);
    TEST_EXIT(solverName != "H2SO4")("You have to specify the parameter %s in the initfile!\n", initFileStr.c_str());
    
    solverType = "mtl_" + solverName;

    LinearSolverCreator *solverCreator = 
	    dynamic_cast<LinearSolverCreator*>(CreatorMap<LinearSolverInterface>::getCreator(solverType, initFileStr));
    TEST_EXIT(solverCreator)("No valid solver type found in parameter \"%s\"\n", initFileStr.c_str());
    
    solverCreator->setName(initFileStr);

    solver = dynamic_cast<LinearSolverInterface*>(solverCreator->create());
    TEST_EXIT(solver)("Can't create Solver correctly\n");

    runner = dynamic_cast< RunnerBase<SparseMatrix, DenseVector>* >(solver->getRunner());
    TEST_EXIT(runner)("Can't create Runner from %s correctly\n", solver->getName().c_str());
  }

  void init(SparseMatrix *matrix) {
    mat = matrix;
    SolverMatrix<Matrix<DOFMatrix*> > dummy;
    runner->init(dummy, *mat);
  }

  void solve(const DenseVector &b, DenseVector &x) const {
    runner->solve(*mat, x, b);
  }

  std::string getSolverName() {
    return solverType;
  }


private:
  SparseMatrix *mat;

  LinearSolverInterface* solver;
  
  RunnerBase<SparseMatrix, DenseVector> *runner;

  std::string solverType;
};

}}

#endif
