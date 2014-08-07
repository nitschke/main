#include "AMDiS.h"
#include "WorldVectorHelper.h"

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

  DOFVector<double> minus(const DOFVector<WorldVector<double> >& x, const DOFVector<WorldVector<double> >& y)
  {
    FUNCNAME("DOFVector<WorldVector<double> >::minus(const DOFVector<WorldVector<double> >& x, const DOFVector<WorldVector<double> >& y)");

    DOFVector<WorldVector<double> >::Iterator xIterator(const_cast<DOFVector<WorldVector<double> >*>(&x), USED_DOFS);
    DOFVector<WorldVector<double> >::Iterator yIterator(const_cast<DOFVector<WorldVector<double> >*>(&y), USED_DOFS);

    DOFVector<double> rval(x.getFeSpace(), "mag(x-y)");
    DOFVector<double>::Iterator rvalIterator(&rval, USED_DOFS);
    
    for (xIterator.reset(), yIterator.reset(), rvalIterator.reset(); !xIterator.end();
	 ++xIterator, ++yIterator, ++rvalIterator) {
      *rvalIterator = wvnorm(*xIterator - *yIterator); 
    }

    return rval;
  }

  DOFVector<double> mag(const DOFVector<WorldVector<double> >& v)
  {
    FUNCNAME("DOFVector<double>::halfMag(const DOFVector<double>& x, const DOFVector<double>& y, const DOFVector<double>& z)");

    DOFVector<WorldVector<double> >::Iterator vIterator(const_cast<DOFVector<WorldVector<double> >*>(&v), USED_DOFS);

    DOFVector<double> rval(v.getFeSpace(), "mag(v)");
    DOFVector<double>::Iterator rvalIterator(&rval, USED_DOFS);
    
    for (vIterator.reset(), rvalIterator.reset(); !vIterator.end(); ++vIterator, ++rvalIterator) {
      *rvalIterator = wvnorm(*vIterator); 
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

    DOFVector<double> rval(v.getFeSpace(), "v0_Times_v1");
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

    DOFVector<double> rval(v.getFeSpace(), "HalfOf_v0_Plus_v1");
    DOFVector<double>::Iterator rvalIterator(&rval, USED_DOFS);
    
    for (vIterator.reset(), rvalIterator.reset(); !vIterator.end();
	 ++vIterator, ++rvalIterator) {
      *rvalIterator = 0.5 * ((*vIterator)[0] + (*vIterator)[1]); 
    }

    return rval;
  }

  DOFVector<WorldVector<double> > toWorld(const DOFVector<double>& x, const DOFVector<double>& y, const DOFVector<double>& z)
  {
    FUNCNAME("DOFVector<WorldVector<double> > toWorld(const DOFVector<double>& x, const DOFVector<double>& y, const DOFVector<double>& z)");

    DOFVector<double>::Iterator xIterator(const_cast<DOFVector<double>*>(&x), USED_DOFS);
    DOFVector<double>::Iterator yIterator(const_cast<DOFVector<double>*>(&y), USED_DOFS);
    DOFVector<double>::Iterator zIterator(const_cast<DOFVector<double>*>(&z), USED_DOFS);

    DOFVector<WorldVector<double> > rval(x.getFeSpace(), "[x,y,z]");
    DOFVector<WorldVector<double> >::Iterator rvalIterator(&rval, USED_DOFS);
    
    for (xIterator.reset(), yIterator.reset(), zIterator.reset(), rvalIterator.reset(); !xIterator.end();
	 ++xIterator, ++yIterator, ++zIterator, ++rvalIterator) {
      (*rvalIterator)[0] = *xIterator;
      (*rvalIterator)[1] = *yIterator;
      (*rvalIterator)[2] = *zIterator;
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

  void printError(SystemVector &sysv, int i0 ,int i1, int i2, const DOFVector<WorldVector<double> > &sol, string name) {
    DOFVector<WorldVector<double> > dofv = toWorld(*sysv.getDOFVector(i0), *sysv.getDOFVector(i1), *sysv.getDOFVector(i2));
    DOFVector<double> err = minus(dofv,sol);
    DOFVector<double> solMag = mag(sol);
    cout << "************ " << name << " ************" << endl;
    cout << "L2Err_" << name << "....... " << err.L2Norm() << endl;
    cout << "MaxErr_" << name << "...... " << err.absMax() << endl;
    cout << "L2ErrRel_" << name << ".... " << (err.L2Norm() / solMag.L2Norm()) << endl;
    cout << "MaxErrRel_" << name << "... " << (err.absMax() / solMag.absMax()) << endl << endl;
  }


}
