#include "AMDiS.h"

namespace AMDiS {

  template<typename T>
  const DOFVector<T>& minus(const DOFVector<T>& x, const DOFVector<T>& y)
  {
    FUNCNAME_DBG("DOFVector<T>::operator-=(DOFVector<T>& x, const DOFVector<T>& y)");

    TEST_EXIT_DBG(x.getFeSpace() && y.getFeSpace())
      ("feSpace is NULL: %8X, %8X\n", x.getFeSpace(), y.getFeSpace());
    TEST_EXIT_DBG(x.getFeSpace()->getAdmin() &&
	      (x.getFeSpace()->getAdmin() == y.getFeSpace()->getAdmin()))
      ("no admin or different admins: %8X, %8X\n",
       x.getFeSpace()->getAdmin(), y.getFeSpace()->getAdmin());
    TEST_EXIT_DBG(x.getSize() == y.getSize())("different sizes\n");
    
    typename DOFVector<T>::Iterator xIterator(dynamic_cast<DOFIndexed<T>*>(const_cast<DOFVector<T>*>(&x)), USED_DOFS);
    typename DOFVector<T>::Iterator yIterator(dynamic_cast<DOFIndexed<T>*>(const_cast<DOFVector<T>*>(&y)), USED_DOFS);
    for (xIterator.reset(), yIterator.reset(); !xIterator.end();
	 ++xIterator, ++yIterator)
      *xIterator -= *yIterator; 

    return x;
  }


  inline void printError(const DOFVector<double> &dofv,const DOFVector<double> &sol, string name) {
    DOFVector<double> m = minus(dofv,sol);
    cout << "MaxErr_" << name << ": " << m.absMax() << endl;
    cout << "L2Err_" << name << ":  " << m.l2norm() << endl;
  }
}
