#include "AMDiS.h"
#include "decOperator.h"
#include "phiProjection.h"
#include "dummyProjection.h"
#include "torusProjection.h"
#include "meshCorrector.h"
#include "MeshHelper.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  new DummyProject(1, VOLUME_PROJECTION);
  //new TorusProject(1, VOLUME_PROJECTION);
  //WorldVector<double> ballCenter;
  //ballCenter.set(0.0);
  //new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);
  

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);


  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
					       &sphere,
					       adaptInfo);
  
  double h;
  Parameters::get("meshCorrector->h", h);
  int nMax;
  Parameters::get("meshCorrector->nMax", nMax);

  DOFVector<int> conn = getConnections(sphere.getFeSpace());
  VtkVectorWriter::writeFile(conn, "output/Connections.vtu");
  DOFVector<int>::Iterator connIter(const_cast<DOFVector<int>*>(&conn), USED_DOFS);
  int counter= 0;
  cout << endl;
  bool allesImLotAufmBoot = true;
  for (connIter.reset(); !connIter.end(); ++connIter, ++counter)
    if (*connIter < 5) {
      cout << "Only " << *connIter << " elements on vertex " << counter << endl;
      allesImLotAufmBoot = false;
    }

  MeshCorrector mc(sphere.getFeSpace());
  if (allesImLotAufmBoot) mc.iterate(nMax, h, "bunny");


  // ===== start adaption loop =====
  adapt->adapt();


  string meshOut;
  Parameters::get("meshCorrector->outName", meshOut);
  //MacroWriter::writeMacro(new DataCollector<double>(sphere.getFeSpace()), meshOut.c_str());
  //DOFVector<double> *phi = new DOFVector<double>(sphere.getFeSpace(),"phi");
  //phi->interpol(new Phi());
  //VtkVectorWriter::writeFile(phi, meshOut + string("_phi.vtu"));

  AMDiS::finalize();
}


