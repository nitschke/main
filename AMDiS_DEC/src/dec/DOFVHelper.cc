#include "AMDiS.h"

namespace AMDiS {

  DOFVector<double> minus(const DOFVector<double>& x, const DOFVector<double>& y)
  {
    FUNCNAME("DOFVector<double>::minus(const DOFVector<double>& x, const DOFVector<double>& y)");

    DOFVector<double>::Iterator xIterator(const_cast<DOFVector<double>*>(&x), USED_DOFS);
    DOFVector<double>::Iterator yIterator(const_cast<DOFVector<double>*>(&y), USED_DOFS);

    DOFVector<double> rval(x.getFeSpace(), "x-y");
    DOFVector<double>::Iterator rvalIterator(&rval, USED_DOFS);
    
    for (xIterator.reset(), yIterator.reset(), rvalIterator.reset(); !xIterator.end();
	 ++xIterator, ++yIterator, ++rvalIterator) {
      *rvalIterator = *xIterator - *yIterator; 
    }

    return rval;
  }

  DOFVector<double> halfMag(const DOFVector<double>& x, const DOFVector<double>& y, const DOFVector<double>& z)
  {
    FUNCNAME("DOFVector<double>::halfMag(const DOFVector<double>& x, const DOFVector<double>& y, const DOFVector<double>& z)");

    DOFVector<double>::Iterator xIterator(const_cast<DOFVector<double>*>(&x), USED_DOFS);
    DOFVector<double>::Iterator yIterator(const_cast<DOFVector<double>*>(&y), USED_DOFS);
    DOFVector<double>::Iterator zIterator(const_cast<DOFVector<double>*>(&z), USED_DOFS);

    DOFVector<double> rval(x.getFeSpace(), "0.5*mag(x,y,z)");
    DOFVector<double>::Iterator rvalIterator(&rval, USED_DOFS);
    
    for (xIterator.reset(), yIterator.reset(), zIterator.reset(), rvalIterator.reset(); !xIterator.end();
	 ++xIterator, ++yIterator, ++zIterator, ++rvalIterator) {
      *rvalIterator = 0.5 * sqrt((*xIterator)*(*xIterator) + (*yIterator)*(*yIterator) + (*zIterator)*(*zIterator)); 
    }

    return rval;
  }


  DOFVector<double> prod01(const DOFVector<WorldVector<double> >& v)
  {
    FUNCNAME("DOFVector<double>::prod01(const DOFVector<WorldVector<double> >& v");

    DOFVector<WorldVector<double> >::Iterator vIterator(const_cast<DOFVector<WorldVector<double> >*>(&v), USED_DOFS);

    DOFVector<double> rval(v.getFeSpace(), "v[0]*v[1]");
    DOFVector<double>::Iterator rvalIterator(&rval, USED_DOFS);
    
    for (vIterator.reset(), rvalIterator.reset(); !vIterator.end();
	 ++vIterator, ++rvalIterator) {
      *rvalIterator = (*vIterator)[0] * (*vIterator)[1]; 
    }

    return rval;
  }

  DOFVector<double> halfSum01(const DOFVector<WorldVector<double> >& v)
  {
    FUNCNAME("DOFVector<double>::halfSum(const DOFVector<WorldVector<double> >& v");

    DOFVector<WorldVector<double> >::Iterator vIterator(const_cast<DOFVector<WorldVector<double> >*>(&v), USED_DOFS);

    DOFVector<double> rval(v.getFeSpace(), "0.5*(v[0]+v[1])");
    DOFVector<double>::Iterator rvalIterator(&rval, USED_DOFS);
    
    for (vIterator.reset(), rvalIterator.reset(); !vIterator.end();
	 ++vIterator, ++rvalIterator) {
      *rvalIterator = 0.5 * ((*vIterator)[0] + (*vIterator)[1]); 
    }

    return rval;
  }


// alte berechnungen faelschlicherweise mit l2norm -_-
  void printError(const DOFVector<double> &dofv,const DOFVector<double> &sol, string name) {
    DOFVector<double> err = minus(dofv,sol);
    cout << "************ " << name << " ************" << endl;
    cout << "L2Err_" << name << "....... " << err.L2Norm() << endl;
    cout << "MaxErr_" << name << "...... " << err.absMax() << endl;
    cout << "L2ErrRel_" << name << ".... " << (err.L2Norm() / sol.L2Norm()) << endl;
    cout << "MaxErrRel_" << name << "... " << (err.absMax() / sol.absMax()) << endl << endl;
  }


}
