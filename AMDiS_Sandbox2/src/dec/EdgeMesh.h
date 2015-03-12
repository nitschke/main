#ifndef EDGEMESH_H
#define EDGEMESH_H

#include "AMDiS.h"
#include "elVolumesInfo2d.h"


using namespace AMDiS;
using namespace std;


/*          /\
 *         /  \
 *        /    \
 *       /      \
 *      /  left  \
 *first ---------> second
 *      \  right /
 *       \      /
 *        \    /
 *         \  /
 *          \/
 */
struct EdgeElement {
  DegreeOfFreedom edgeDof; // aka int  -> edge index
  DofEdge dofEdge;         // aka pair<int,int> -> (first,second) vertex index

  ElVolumesInfo2d *infoLeft;
  ElVolumesInfo2d *infoRight;
};


class EdgeMesh {

public:
  // works only with no boundaries
  EdgeMesh(const FiniteElemSpace *feSpace_);

  const FiniteElemSpace* getFeSpace() const {return feSpace;}

  int getNumberOfEdges() const {return nEdges;} 
  

  const vector<EdgeElement>* getEdges() const {return &edges;}

  DOFVector< list<EdgeElement> > getEdgeRings() const;

 



private:
  const FiniteElemSpace *feSpace;
  vector<EdgeElement> edges;

  int nEdges;
};

#endif
