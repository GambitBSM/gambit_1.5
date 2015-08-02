#include "Pythia8/ResonanceWidths.h"
#include "backend_types/Pythia_8_209/wrapper_ResonanceGmZ.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_ResonanceGmZ* Factory_ResonanceGmZ_0(int idResIn)
    {
        return new ResonanceGmZ(idResIn);
    }
    
}

