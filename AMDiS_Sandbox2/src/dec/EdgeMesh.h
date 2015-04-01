#ifndef EDGEMESH_H
#define EDGEMESH_H

#include "AMDiS.h"
#include "elVolumesInfo2d.h"


using namespace AMDiS;
using namespace std;


/*                  /\
 *                 /  \
 *     left.first /    \ left.second
 *               /      \
 *              /  left  \
 *        first ---------> second
 *              \  right /
 *               \      /
 *    right.first \    / right.second
 *                 \  /
 *                  \/
 */
struct EdgeElement {
  DegreeOfFreedom edgeDof; // aka int  -> edge index
  DofEdge dofEdge;         // aka pair<int,int> -> (first,second) vertex index

  ElVolumesInfo2d *infoLeft;
  pair<EdgeElement*, EdgeElement*> edgesLeft;

  ElVolumesInfo2d *infoRight;
  pair<EdgeElement*, EdgeElement*> edgesRight;
};


class EdgeMesh {

public:
  // works only with no boundaries
  EdgeMesh(const FiniteElemSpace *feSpace_);

  const FiniteElemSpace* getFeSpace() const {return feSpace;}

  int getNumberOfEdges() const {return nEdges;} 
  

  const vector<EdgeElement>* getEdges() const {return &edges;}

  DOFVector< list<EdgeElement> > getEdgeRings() const;

  map<int, pair<ElVolumesInfo2d*, vector<EdgeElement> > > getFaceEdges() const;

  // volume of the edge skeleton (Vol(|K^(1)|)...length of all edges)
  double getVol() const;

 



private:
  const FiniteElemSpace *feSpace;
  vector<EdgeElement> edges;

  int nEdges;
};

inline ostream &operator <<(ostream &out, const EdgeElement &eel) {
  out << eel.edgeDof << "(" << eel.dofEdge.first << "," << eel.dofEdge.second << "): " << endl;
  out << "left  face: " << *(eel.infoLeft) << endl;
  out << "right face: " << *(eel.infoRight) << endl;
  out << "left  edges (first, second): " << eel.edgesLeft.first->edgeDof << " " <<  eel.edgesLeft.second->edgeDof << endl;
  out << "right edges (first, second): " << eel.edgesRight.first->edgeDof << " " <<  eel.edgesRight.second->edgeDof << endl;
  return out;
}

inline ostream &operator <<(ostream &out, const EdgeMesh &eMesh) {
  out << eMesh.getNumberOfEdges() << " edges:" << endl;
  vector<EdgeElement>::const_iterator eIter = eMesh.getEdges()->begin();
  for (; eIter != eMesh.getEdges()->end(); ++eIter) {
    out << *eIter;
  }
  return out;
}

#endif
