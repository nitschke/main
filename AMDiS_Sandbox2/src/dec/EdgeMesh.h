#ifndef EDGEMESH_H
#define EDGEMESH_H

#include "Dec_fwd.h"
#include "ElVolumesInfo2d.h" // <<(ElVolumesInfo2d&)


using namespace std;
namespace AMDiS { namespace dec {


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
  //TODO: pointer, so that dofCompress is consistent
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
      EdgeRingIterator(const EdgeElement *eel, EdgeRingIteratorType type) : t(type), reft(type), start(eel), position(eel) {
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

      // return iteratorType, this can be useful, because the type change iff the orientation change.
      EdgeRingIteratorType getType() {
        return t;
      }

      // return true if the edge sign induced from vertex resp. face differ from the beginning sign.
      // This mean for ...
      // FIRSTVERTEX : true if edge point to reference vertex
      // SECONDVERTEX: false if edge point to reference vertex
      // LEFTFACE    : false if edge point in counterclockwise direction
      // RIGHTFACE   : true if edge point in counterclockwise direction
      bool changedSign() {
        return t != reft;
      }

      // true if the edge point outward from the reference Vertex
      bool pointOutward() {
        if (reft == FIRSTVERTEX) return (t == reft);
        if (reft == SECONDVERTEX) return (t != reft);
        ERROR_EXIT("Wrong iterator type for EdgeElement::EdgeRingIterator::pointOutward. Only *VERTEX iteration are permitted.");
      }

      bool pointCounterclockwise() {
        if (reft == LEFTFACE) return t == reft;
        if (reft == RIGHTFACE) return t != reft;
        ERROR_EXIT("Wrong iterator type for EdgeElement::EdgeRingIterator::pointCounterclockwise. Only *FACE iteration are permitted.");
      }

      // return the reference face for a *FACE iteration
      // or the next face in iteration direction (counterclockwise) for a *VERTEX iteration 
      const ElVolumesInfo2d* getFace() {
        if (t == RIGHTFACE || t == LEFTFACE) {
          return refFace;
        } else {
          if (t == FIRSTVERTEX) return position->infoLeft;
          if (t == SECONDVERTEX) return position->infoRight;
        }
      }

    private:
      const EdgeElement *start;
      const EdgeElement *position;

      EdgeRingIteratorType t;
      EdgeRingIteratorType reft;

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

  // maximum (outer) diameter of all faces
  double getMaxFaceDiameter() const;

  // maximum diameter(=length) of all edges
  double getMaxEdgeDiameter() const;

  //only for debuging
  void printVolInfos(bool ePLen = false, bool eDLen = false) {
  vector<EdgeElement>::const_iterator eIter = edges.begin();
    for (; eIter != edges.end(); ++eIter) {
      cout << eIter->edgeDof << "(" << eIter->dofEdge.first << "," << eIter->dofEdge.second << "): " << endl;
      if (ePLen) cout << " ePLen: " << eIter->infoLeft->getEdgeLen(eIter->dofEdge) << endl;
      if (eDLen) cout << " eDLen: " << (eIter->infoLeft->getDualEdgeLen(eIter->dofEdge) 
                                       +eIter->infoRight->getDualEdgeLen(eIter->dofEdge)) << endl;
    }
  }

  ~EdgeMesh() {
    vector< ElVolumesInfo2d* > elVols(feSpace->getMesh()->getMacroElements().size());
    vector<EdgeElement>::iterator eIter = edges.begin();
    for (; eIter != edges.end(); ++eIter) {
      elVols[eIter->infoLeft->getElInfo()->getMacroElement()->getIndex()] = eIter->infoLeft;
      elVols[eIter->infoRight->getElInfo()->getMacroElement()->getIndex()] = eIter->infoRight;
    }
    vector< ElVolumesInfo2d* >::iterator infoIter = elVols.begin();
    for (; infoIter != elVols.end(); ++infoIter) delete (*infoIter);
  }
 



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

}}
#endif
