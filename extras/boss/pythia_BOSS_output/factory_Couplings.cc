#include "Pythia8/StandardModel.h"
#include "backend_types/Pythia_8_186/wrapper_Couplings.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_Couplings* Factory_Couplings_0()
    {
        return new Couplings();
    }
    
}

