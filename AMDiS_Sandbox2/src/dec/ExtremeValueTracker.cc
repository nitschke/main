#include "ExtremeValueTracker.h"
#include "EdgeMesh.h"
#include <limits>

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

int OneMinValTrackerInPositiveZ::trackdownMinima(const DofEdgeVector &dofe, double time) {
  std::list< std::pair< double , double > > ZsVals;
  vector<EdgeElement>::const_iterator edgeIter = dofe.getEdgeMesh()->getEdges()->begin();
  for (; edgeIter != dofe.getEdgeMesh()->getEdges()->end(); ++edgeIter) {
    if(isMinimum(dofe, *edgeIter)) {
      WorldVector<double> ec = edgeIter->infoLeft->getEdgeCenter(edgeIter->dofEdge);
      if (ec[2] >= 0.0) ZsVals.push_back(std::make_pair(ec[2], dofe[*edgeIter]));
    }
  }

  const double max_double = std::numeric_limits<double>::max();
  pair< double , double > overallMin(0. , max_double);
  for (std::list< std::pair< double , double > >::iterator it = ZsVals.begin();
       it != ZsVals.end(); ++it) {
    if (it->second < overallMin.second) overallMin = *it;
  }

  csvout << time << "," << overallMin.first << endl;

  return (overallMin.second == max_double) ? 0 : 1;
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

bool ExtremeValueTracker::isMinimum(const DofEdgeVector &dofe, const EdgeElement &eel) const {
  bool isExtreme = true;
    
  EdgeElement::EdgeRingIterator ringIter1(&eel, FIRSTVERTEX);
  for (; !ringIter1.isEnd() && isExtreme; ++ringIter1) {
    if ( dofe[eel] > dofe[*ringIter1] ) isExtreme = false;
  }
  
  EdgeElement::EdgeRingIterator ringIter2(&eel, SECONDVERTEX);
  for (; !ringIter2.isEnd() && isExtreme; ++ringIter2) {
    if ( dofe[eel] > dofe[*ringIter2] ) isExtreme = false;
  }

  return isExtreme;
}
