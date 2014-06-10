#include "Pythia8/Event.h"

Pythia8::Abstract__Particle* Factory_Particle()
{
    return new Pythia8::Particle();
}

Pythia8::Abstract__Particle* Factory_Particle(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, double pxIn, double pyIn, double pzIn, double eIn, double mIn, double scaleIn, double polIn)
{
    return new Pythia8::Particle( idIn,  statusIn,  mother1In,  mother2In,  daughter1In,  daughter2In,  colIn,  acolIn,  pxIn,  pyIn,  pzIn,  eIn,  mIn,  scaleIn,  polIn);
}

Pythia8::Abstract__Particle* Factory_Particle(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, Pythia8::Vec4 pIn, double mIn, double scaleIn, double polIn)
{
    return new Pythia8::Particle( idIn,  statusIn,  mother1In,  mother2In,  daughter1In,  daughter2In,  colIn,  acolIn,  pIn,  mIn,  scaleIn,  polIn);
}

Pythia8::Abstract__Particle* Factory_Particle(const Pythia8::Particle& pt)
{
    return new Pythia8::Particle( pt);
}


