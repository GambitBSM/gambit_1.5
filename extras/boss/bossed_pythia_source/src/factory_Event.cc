#include "Pythia8/Event.h"
#include "backend_types/Pythia_8_209/wrapper_Event.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_Event* Factory_Event_0(int capacity)
    {
        return new Event(capacity);
    }
    
    Abstract_Event* Factory_Event_1()
    {
        return new Event();
    }
    
}

