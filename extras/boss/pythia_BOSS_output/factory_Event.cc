#include "Pythia8/Event.h"
#include "backend_types/Pythia_8_186/wrapper_Event_decl.h"
#include "backend_types/Pythia_8_186/wrapper_Event_def.h"
#include "abstracttypedefs.h"
#include "wrappertypedefs.h"

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

