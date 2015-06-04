#include "Pythia8/SusyCouplings.h"
#include "backend_types/Pythia_8_186/wrapper_CoupSUSY.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_CoupSUSY* Factory_CoupSUSY_0()
    {
        return new CoupSUSY();
    }
    
}

