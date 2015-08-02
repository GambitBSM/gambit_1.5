#include "Pythia8/StandardModel.h"
#include "backend_types/Pythia_8_209/wrapper_CoupSM.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_CoupSM* Factory_CoupSM_0()
    {
        return new CoupSM();
    }
    
}

