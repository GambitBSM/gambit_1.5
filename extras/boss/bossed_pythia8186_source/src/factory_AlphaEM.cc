#include "Pythia8/StandardModel.h"
#include "backend_types/Pythia_8_186/wrapper_AlphaEM.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_AlphaEM* Factory_AlphaEM_0()
    {
        return new AlphaEM();
    }
    
}

