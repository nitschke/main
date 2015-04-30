#include "AnimationWriter.h"

using namespace std;
using namespace AMDiS;
using namespace dec;

void AnimationWriter::updateAnimationFile(double time, string frame) {
  animFiles.push_back(make_pair(time, frame));

  boost::iostreams::filtering_ostream file;
  file.push(boost::iostreams::file_descriptor_sink(fname, std::ios::trunc));

  file << "<?xml version=\"1.0\"?>\n";
	file << "<VTKFile type=\"Collection\" version=\"0.1\" >"  << "\n";
	file << "<Collection>\n";

  for (list< pair<double,string> >::iterator it = animFiles.begin();
       it != animFiles.end(); ++it) {
    file << "<DataSet timestep=\"" << it->first
		<< "\" part=\"0\" file=\"" << it->second << "\"/>\n";      
  }

  file << "</Collection>\n";
	file << "</VTKFile>\n";
}
