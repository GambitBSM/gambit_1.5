#include "Pythia8/ParticleDecays.h"
#include "backend_types/Pythia_8_186/wrapper_ParticleDecays.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_ParticleDecays* Factory_ParticleDecays_0()
    {
        return new ParticleDecays();
    }
    
}

