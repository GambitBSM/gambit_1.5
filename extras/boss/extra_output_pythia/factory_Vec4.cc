#include "Pythia8/Basics.h"

Pythia8::Abstract__Vec4* Factory_Vec4(double xIn, double yIn, double zIn, double tIn)
{
    return new Pythia8::Vec4( xIn,  yIn,  zIn,  tIn);
}

Pythia8::Abstract__Vec4* Factory_Vec4(const Pythia8::Vec4& v)
{
    return new Pythia8::Vec4( v);
}


