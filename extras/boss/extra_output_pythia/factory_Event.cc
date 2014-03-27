#include "Pythia8/Event.h"

Pythia8::Abstract__Event* Factory_Event(int capacity)
{
    return new Pythia8::Event( capacity);
}

Pythia8::Abstract__Event* Factory_Event(const Pythia8::Event& oldEvent)
{
    return new Pythia8::Event( oldEvent);
}


