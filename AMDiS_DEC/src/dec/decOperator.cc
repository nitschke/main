#include "decOperator.h"
#include "elVolumesInfo2d.h"

namespace AMDiS {

  using namespace std;

  void LBeltramiDEC::getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor) {
    
    ElVolumesInfo2d volInfo(elInfo);
    
    userMat = 0.0;

    for (int i = 0; i < 3; i++) {
      for (int j = (i+1)%3; j != i; j = (j+1)%3) {
        int k = (2*(i+j))%3;
        double c = volInfo.getDualOppEdgeLen(k) / (volInfo.getOppEdgeLen(k));
        userMat(i, i) = userMat(i, i) - c;
        userMat(i, j) = userMat(i, j) + c;
      }
    }

    userMat *= opFactor * factor;
    cout << "Beltrami: \n" << userMat << endl;
  }

  void FunctionDEC::getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor) {

    ElVolumesInfo2d volInfo(elInfo);

    for (int i = 0; i < 3; i++) {
      userVec(i) = volInfo.getDualVertexVol(i)*(*f)(elInfo->getCoord(i));
    }

    userVec *= opFactor * factor;

  }

  // TODO: untested
  void FunctionDEC::getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor) {
    
    ElVolumesInfo2d volInfo(elInfo);
    
    userMat = 0.0;

    for (int i = 0; i < 3; i++) {
      userMat(i, i) = volInfo.getDualVertexVol(i)*(*f)(elInfo->getCoord(i));
    }
    
    userMat *= opFactor * factor;
  }

  
  void SimpleDEC::getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor) {
    
    ElVolumesInfo2d volInfo(elInfo);
    
    userMat = 0.0;

    for (int i = 0; i < 3; i++) {
      userMat(i, i) = volInfo.getDualVertexVol(i);
    }
    
    //cout << factor << endl;
    //cout << opFactor << endl;
    userMat *= opFactor * factor;
    cout << "SimpleDEC:\n" << userMat << endl;
  }


  void SimpleDEC::getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor) {

    ElVolumesInfo2d volInfo(elInfo);

    if (uhOld) {
      //cout << "*";
      uhOld->getLocalVector(elInfo->getElement(), userVec);
      //cout << userVec << endl;
    } else {
      userVec = 0.0;
    }

    for (int i = 0; i < 3; i++) {
      userVec(i) *= volInfo.getDualVertexVol(i);
    }


    userVec *= opFactor * factor;

  }

}
