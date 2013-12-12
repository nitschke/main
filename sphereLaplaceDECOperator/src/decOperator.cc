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
      double volStarP = volInfo.getDualVertexVol(i);
      for (int j = (i+1)%3; j != i; j = (j+1)%3) {
        int k = (2*(i+j))%3;
        double c = volInfo.getDualOppEdgeLen(k) / (volInfo.getOppEdgeLen(k) * volStarP);
        userMat(i, i) = userMat(i, i) - c;
        userMat(i, j) = userMat(i, j) + c;
      }
    }
    //cout << userMat << endl;
  }

  void FunctionDEC::getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor) {

    userVec = 0.0;

    for (int i = 0; i < 3; i++) {
      userVec(i) = (*f)(elInfo->getCoord(i));
    }

  }

}
