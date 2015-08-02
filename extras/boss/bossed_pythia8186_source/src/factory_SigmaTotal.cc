#include "Pythia8/SigmaTotal.h"
#include "backend_types/Pythia_8_186/wrapper_SigmaTotal.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_SigmaTotal* Factory_SigmaTotal_0()
    {
        return new SigmaTotal();
    }
    
}

