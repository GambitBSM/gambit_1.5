#include "Pythia8/BeamParticle.h"
#include "backend_types/Pythia_8_209/wrapper_BeamParticle.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_BeamParticle* Factory_BeamParticle_0()
    {
        return new BeamParticle();
    }
    
}

