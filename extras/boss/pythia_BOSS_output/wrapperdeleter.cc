#include "wrapperdeleter.h"

void wrapper_deleter(Pythia8::Particle__BOSS* wptr)
{
    delete wptr;
}

void wrapper_deleter(Pythia8::Info__BOSS* wptr)
{
    delete wptr;
}

void wrapper_deleter(Pythia8::Vec4__BOSS* wptr)
{
    delete wptr;
}

void wrapper_deleter(Pythia8::Hist__BOSS* wptr)
{
    delete wptr;
}

void wrapper_deleter(Pythia8::Event__BOSS* wptr)
{
    delete wptr;
}

void wrapper_deleter(Pythia8::Pythia__BOSS* wptr)
{
    delete wptr;
}
