#include "wrapperdeleter.h"

void wrapper_deleter(Pythia8::Particle_GAMBIT* wptr)
{
    delete wptr;
}

void wrapper_deleter(Pythia8::Info_GAMBIT* wptr)
{
    delete wptr;
}

void wrapper_deleter(Pythia8::Vec4_GAMBIT* wptr)
{
    delete wptr;
}

void wrapper_deleter(Pythia8::Hist_GAMBIT* wptr)
{
    delete wptr;
}

void wrapper_deleter(Pythia8::Event_GAMBIT* wptr)
{
    delete wptr;
}

void wrapper_deleter(Pythia8::Pythia_GAMBIT* wptr)
{
    delete wptr;
}
