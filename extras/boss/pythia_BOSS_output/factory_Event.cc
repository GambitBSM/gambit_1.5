#include "Pythia8/Event.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Event_decl.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Event_def.h"
#include "abstracttypedefs.h"
#include "wrappertypedefs.h"

namespace Pythia8
{
    Abstract_Event* Factory_Event(int capacity)
    {
        return new Event(capacity);
    }
    
    Abstract_Event* Factory_Event()
    {
        return new Event();
    }
    
}

