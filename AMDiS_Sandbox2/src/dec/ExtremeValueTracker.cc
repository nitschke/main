#include "ExtremeValueTracker.h"
#include "EdgeMesh.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

int ExtremeValueTracker::trackdownMaxima(const DofEdgeVector &dofe, double time, double minVal) {
  int counter = 0;
  vector<EdgeElement>::const_iterator edgeIter = dofe.getEdgeMesh()->getEdges()->begin();
  for (; edgeIter != dofe.getEdgeMesh()->getEdges()->end(); ++edgeIter) {
    if (dofe[*edgeIter] > minVal && isMaximum(dofe, *edgeIter)) {
      WorldVector<double> ec = edgeIter->infoLeft->getEdgeCenter(edgeIter->dofEdge);
      csvout << time << "," << ec[0] << "," << ec[1] << "," << ec[2] << endl;
      counter++;
    } 
  }
  return counter;
}

bool ExtremeValueTracker::isMaximum(const DofEdgeVector &dofe, const EdgeElement &eel) const {
  bool isExtreme = true;
    
  EdgeElement::EdgeRingIterator ringIter1(&eel, FIRSTVERTEX);
  for (; !ringIter1.isEnd() && isExtreme; ++ringIter1) {
    if ( dofe[eel] < dofe[*ringIter1] ) isExtreme = false;
  }
  
  EdgeElement::EdgeRingIterator ringIter2(&eel, SECONDVERTEX);
  for (; !ringIter2.isEnd() && isExtreme; ++ringIter2) {
    if ( dofe[eel] < dofe[*ringIter2] ) isExtreme = false;
  }

  return isExtreme;
}
