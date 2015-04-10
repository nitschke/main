#ifndef EDGEMESH_H
#define EDGEMESH_H

#include "AMDiS.h"
#include "elVolumesInfo2d.h"


using namespace AMDiS;
using namespace std;

typedef enum {
  FIRSTVERTEX = 1,
  SECONDVERTEX = 2,
  LEFTFACE = 3,
  RIGHTFACE = 4
} EdgeRingIteratorType;

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

  // const_iterator over edges (always counterclockwise)
  // FIRSTVERTEX : iterator over the edge ring of edge.first vertex
  // SECONDVERTEX: iterator over the edge ring of edge.second vertex
  // LEFTFACE    : iterator over the edge ring of the left face (i.e. 3 edges)
  // RIGHTFACE   : iterator over the edge ring of the right face (i.e. 3 edges)
  class EdgeRingIterator {
    public:
      EdgeRingIterator(const EdgeElement *eel, EdgeRingIteratorType type) : t(type), start(eel), position(eel) {
        switch(t) {
          case FIRSTVERTEX: refVert = &(eel->dofEdge.first); break;
          case SECONDVERTEX: refVert = &(eel->dofEdge.second); break;
          case LEFTFACE: refFace = eel->infoLeft; break;
          case RIGHTFACE: refFace = eel->infoRight; break;
          default: ERROR_EXIT("unknow iterator type for EdgeRingIterator");
        }
        counter = 0;
      }

      // Prefix operator++
      const EdgeRingIterator& operator++() {
        switch(t) {
          case FIRSTVERTEX: 
                    position = position->edgesLeft.first; 
                    if (position->dofEdge.second == *refVert) t = SECONDVERTEX;
                    break;
          case SECONDVERTEX: 
                    position = position->edgesRight.second; 
                    if (position->dofEdge.first == *refVert) t = FIRSTVERTEX;
                    break;
          case LEFTFACE:
                    position = position->edgesLeft.second; 
                    if (position->infoRight == refFace) t = RIGHTFACE;
                    break;
          case RIGHTFACE:
                    position = position->edgesRight.first; 
                    if (position->infoLeft == refFace) t = LEFTFACE;
                    break;
          default: ERROR_EXIT("unknow iterator type for EdgeRingIterator");
        }
        counter++;
        return *this;
      }

      // dereference operator*
      const EdgeElement operator*() {
        return *position;
      }

      // dereference operator->
      const EdgeElement* operator->() {
        return position;
      }

      // only true after every cycle
      bool isEnd() {
        return start == position && counter > 0;
      }

      // after one cycle: counter gives the number of edges in the ring
      unsigned short getCounter() {
        return counter;
      }

    private:
      const EdgeElement *start;
      const EdgeElement *position;

      EdgeRingIteratorType t;

      const DegreeOfFreedom *refVert = NULL;
      const ElVolumesInfo2d *refFace = NULL;

      unsigned short counter;
  };
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
