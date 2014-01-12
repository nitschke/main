#include "decOperator.h"
#include "elVolumesInfo2d.h"

namespace AMDiS {

  using namespace std;

  inline void DecOperator::updateUserMat(ElementMatrix& userMat, ElementMatrix& opMat) {
    userMat += opFactor * opMat;
  }

  inline void DecOperator::updateUserVec(ElementVector& userVec, ElementVector& opVec) {
    userVec += opFactor * opVec;
  }


  void LBeltramiDEC::getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor) {
    
    ElVolumesInfo2d volInfo(elInfo);
    opMat = 0.0;
    
    for (int i = 0; i < 3; i++) {
      for (int j = (i+1)%3; j != i; j = (j+1)%3) {
        int k = (2*(i+j))%3;
        double c = volInfo.getDualOppEdgeLen(k) / (volInfo.getOppEdgeLen(k));
        opMat(i, i) = opMat(i, i) - c;
        opMat(i, j) = opMat(i, j) + c;
      }
    }

    opMat *= factor;
    //cout << "Beltrami: \n" << opMat << endl;
    updateUserMat(userMat, opMat);
  }

  void JacobianDEC::getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor) {
    //cout << "*** TIC ***\n";
    
    ElVolumesInfo2d volInfo(elInfo);
    opMat = 0.0;

    double psiCC = (*psi)(volInfo.getCircumcenter());
    
    for (int i = 0; i < 3; i++) {
      // next point
      int iPNext = (i+1)%3;
      // edge [p_i, p_iPNext]
      int kENext = (2*(i+iPNext))%3;
      double c = volInfo.getDualOppEdgeLen(kENext) / (volInfo.getOppEdgeLen(kENext));
      opMat(i, i)      = opMat(i, i) + c * psiCC;
      opMat(i, iPNext) = opMat(i, iPNext) + c * psiCC;

      // prev point
      int iPPrev = (i+2)%3;
      // edge [p_i, p_iPPrev]
      int kEPrev = (2*(i+iPPrev))%3;
      c = volInfo.getDualOppEdgeLen(kEPrev) / (volInfo.getOppEdgeLen(kEPrev));
      opMat(i, i)      = opMat(i, i) - c * psiCC;
      opMat(i, iPPrev) = opMat(i, iPPrev) - c * psiCC;
    }
    
    opMat *= factor;
    updateUserMat(userMat, opMat);
    //cout << "*** TOC ***\n";
  }

  void FunctionDEC::getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor) {

    ElVolumesInfo2d volInfo(elInfo);

    for (int i = 0; i < 3; i++) {
      opVec(i) = volInfo.getDualVertexVol(i)*(*f)(elInfo->getCoord(i));
    }

    opVec *= factor;
    updateUserVec(userVec, opVec);
  }

  // TODO: untested
  void FunctionDEC::getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor) {
    
    ElVolumesInfo2d volInfo(elInfo);
    
    opMat = 0.0;

    for (int i = 0; i < 3; i++) {
      opMat(i, i) = volInfo.getDualVertexVol(i)*(*f)(elInfo->getCoord(i));
    }
    
    opMat *= factor;
    updateUserMat(userMat, opMat);
  }

  
  void SimpleDEC::getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor) {
    
    ElVolumesInfo2d volInfo(elInfo);
    
    opMat = 0.0;

    for (int i = 0; i < 3; i++) {
      opMat(i, i) = volInfo.getDualVertexVol(i);
    }
    
    opMat *= factor;
    //cout << "SimpleDEC:\n" << opMat << endl;
    updateUserMat(userMat, opMat);
  }


  void SimpleDEC::getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor) {

    ElVolumesInfo2d volInfo(elInfo);

    if (uhOld) {
      //cout << "*";
      uhOld->getLocalVector(elInfo->getElement(), opVec);
      //cout << opVec << endl;
    } else {
      opVec = 1.0;
    }

    for (int i = 0; i < 3; i++) {
      opVec(i) *= volInfo.getDualVertexVol(i);
    }


    opVec *= factor;
    updateUserVec(userVec, opVec);
  }

}
