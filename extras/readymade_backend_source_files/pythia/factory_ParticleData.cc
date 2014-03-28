#include "Pythia8/ParticleData.h"

Pythia8::Abstract__ParticleData* Factory_ParticleData(const Pythia8::ParticleData& arg_1)
{
    return new Pythia8::ParticleData( arg_1);
}

Pythia8::Abstract__ParticleData* Factory_ParticleData()
{
    return new Pythia8::ParticleData();
}


