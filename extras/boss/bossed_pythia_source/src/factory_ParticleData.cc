#include "Pythia8/ParticleData.h"
#include "backend_types/Pythia_8_209/wrapper_ParticleData.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_ParticleData* Factory_ParticleData_0()
    {
        return new ParticleData();
    }
    
}

