#include "Pythia8/StandardModel.h"
#include "backend_types/Pythia_8_209/wrapper_AlphaStrong.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_AlphaStrong* Factory_AlphaStrong_0()
    {
        return new AlphaStrong();
    }
    
}

