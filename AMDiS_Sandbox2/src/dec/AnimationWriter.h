#ifndef ANIMATIONWRITER_H
#define ANIMATIONWRITER_H

#include "Dec_fwd.h"

using namespace std;
namespace AMDiS { namespace dec {


class AnimationWriter {
public:

AnimationWriter(string fileName) : fname(fileName) {};

void updateAnimationFile(double time, string frame);

private:

string fname;

list< pair<double,string>  > animFiles; 
};

}}
#endif
