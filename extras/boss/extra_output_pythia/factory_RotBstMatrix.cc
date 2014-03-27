#include "Pythia8/Basics.h"

Pythia8::Abstract__RotBstMatrix* Factory_RotBstMatrix()
{
    return new Pythia8::RotBstMatrix();
}

Pythia8::Abstract__RotBstMatrix* Factory_RotBstMatrix(const Pythia8::RotBstMatrix& Min)
{
    return new Pythia8::RotBstMatrix( Min);
}


