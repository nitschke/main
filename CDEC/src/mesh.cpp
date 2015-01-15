#include "mesh.h"

Mesh::Mesh(const char* meshFile) 
{
  TiXmlDocument *doc = new TiXmlDocument(meshFile);
  bool loadOk = doc->LoadFile();
  if (loadOk) 
  {
    cout << doc->FirstChild(PrimalMesh)->FirstChild(Vertices)->Attribute(NumberOf);
  } 
  else 
  {
    cout << "askjdfsajkdh" << endl;
  }
}
