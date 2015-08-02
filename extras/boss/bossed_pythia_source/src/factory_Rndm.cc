#include "Pythia8/Basics.h"
#include "backend_types/Pythia_8_209/wrapper_Rndm.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_Rndm* Factory_Rndm_0()
    {
        return new Rndm();
    }
    
    Abstract_Rndm* Factory_Rndm_1(int seedIn)
    {
        return new Rndm(seedIn);
    }
    
}

