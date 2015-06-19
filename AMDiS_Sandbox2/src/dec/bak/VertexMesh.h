#ifndef EDGEMESH_H
#define EDGEMESH_H

#include "Dec_fwd.h"
#include "EdgeMesh.h"

using namespace std;
namespace AMDiS { namespace dec {

struct VertexElement {
  DegreeOfFreedom vertexDof; // aka int -> vertex index (same as AMDiS dofs)
  vector<EdgeElement*> eels; // all edges containing this vertex
};

class VertexMesh {

public:
  VertexMesh(const EdgeMesh *edgeMesh);

private:

  const EdgeMesh *eMesh;
  vector<VertexElement> vertices;

};

}}
