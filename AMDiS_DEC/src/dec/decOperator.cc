#include "decOperator.h"
#include "elVolumesInfo2d.h"
#include "WorldVectorHelper.h"

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

  void LBeltramiInteriorFunctionDEC::getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor) {
    
    ElVolumesInfo2d volInfo(elInfo);
    opVec = 0.0;
    
    for (int i = 0; i < 3; i++) {
      for (int j = (i+1)%3; j != i; j = (j+1)%3) {
        int k = (2*(i+j))%3;
        double c = volInfo.getDualOppEdgeLen(k) / (volInfo.getOppEdgeLen(k)) * (*q)(elInfo->getCoord(j));
        opVec(i) -= c;
        opVec(j) += c;
      }
    }

    opVec *= factor;
    updateUserVec(userVec, opVec);
  }


  void LBeltramiWeightedDEC::getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor) {
    
    ElVolumesInfo2d volInfo(elInfo);
    opMat = 0.0;
    
    for (int i = 0; i < 3; i++) {
      for (int j = (i+1)%3; j != i; j = (j+1)%3) {
        int k = (2*(i+j))%3;
        double c = (volInfo.getDualOppEdgeLen(k) / (volInfo.getOppEdgeLen(k)))
                      * 0.5 * ((*kappa)(elInfo->getCoord(i)) + (*kappa)(elInfo->getCoord(j)));
        opMat(i, i) = opMat(i, i) - c;
        opMat(i, j) = opMat(i, j) + c;
      }
    }

    opMat *= factor;
    //cout << "Beltrami: \n" << opMat << endl;
    updateUserMat(userMat, opMat);
  }
  
  void PrimalPrimalGradFunctionDEC::getElementVector(const ElInfo *elInfo, 
		    ElementVector& userVec, 
				double factor) {

    ElVolumesInfo2d volInfo(elInfo);
    opVec = 0.0;

    for (int i = 0; i < 3; i++) {
      double c = volInfo.getDualVertexVol(i);
      for (int j = (i+1)%3; j != i; j = (j+1)%3) {
        opVec[i] += c * (elInfo->getGrdLambda()[j])[l] * ((*f)(elInfo->getCoord(j)) - (*f)(elInfo->getCoord(i)));
      }
    }

    //cout << opVec << endl;
    opVec *= factor;
    updateUserVec(userVec, opVec);
    }

  void PrimalPrimalGradDEC::getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor) {

    ElVolumesInfo2d volInfo(elInfo);
    opMat = 0.0;

    for (int i = 0; i < 3; i++) {
      double c = volInfo.getDualVertexVol(i);
      for (int j = (i+1)%3; j != i; j = (j+1)%3) {
        double Cgrd_jl = c * (elInfo->getGrdLambda()[j])[l];
        opMat[i][j] += Cgrd_jl;
        opMat[i][i] -= Cgrd_jl;
      }
    }

    //cout << opVec << endl;
    opMat *= factor;
    updateUserMat(userMat, opMat);
    }

  void DualPrimalNormalDEC::getElementVector(const ElInfo *elInfo, 
		    ElementVector& userVec, 
				double factor) {

    ElVolumesInfo2d volInfo(elInfo);
    opVec = 0.0;
    
    WorldVector<double> elNormal;
    elInfo->getElementNormal(elNormal);
    for (int i = 0; i < 3; i++) {
      opVec[i] += volInfo.getDualVertexVol(i) * elNormal[l];
    }

    opVec *= factor;
    updateUserVec(userVec, opVec);
  }

  void JacobianDEC::getElementMatrix(const ElInfo *elInfo, 
		    ElementMatrix& userMat, 
				double factor) {
    //cout << "*** TIC ***\n";
    
    ElVolumesInfo2d volInfo(elInfo);
    opMat = 0.0;

    ElementVector phiVec(3);
    for (int i = 0; i < 3; i++) phiVec(i) = (*phi)(elInfo->getCoord(i));
    
    for (int i = 0; i < 3; i++) {
      double c = volInfo.getDualVertexVol(i) / elInfo->getDet();
      for (int k = 0; k < 3; k++) {
        int kk = (k+1)%3;
        opMat(i, k) -= c * phiVec(kk);
        opMat(i, kk) += c * phiVec(k);
      }
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
      //opMat(i, i) = 1.0;
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

  void GaussCurvatureDEC::getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor) {

    ElVolumesInfo2d volInfo(elInfo);
    opVec = 0.0;
    
    //TODO: improve!
    for (int i = 0; i < 3; i++) {
      for (int j = (i+1)%3; j != i; j = (j+1)%3) {
        opVec[i] += 2.0 * M_PI / 12.0 - atan2(volInfo.getDualOppEdgeLen(j), 0.5 * volInfo.getOppEdgeLen(j));
        //cout << atan2(volInfo.getDualOppEdgeLen(j), 0.5 * volInfo.getOppEdgeLen(j)) << endl;
      }
    }

    opVec *= factor;
    updateUserVec(userVec, opVec);


  }

}
