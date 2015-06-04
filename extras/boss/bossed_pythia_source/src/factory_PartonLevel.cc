#include "Pythia8/PartonLevel.h"
#include "backend_types/Pythia_8_186/wrapper_PartonLevel.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_PartonLevel* Factory_PartonLevel_0()
    {
        return new PartonLevel();
    }
    
}

