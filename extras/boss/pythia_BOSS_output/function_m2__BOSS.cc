

#include "Pythia8/Event.h"
#include "backend_types/Pythia_8_209/wrapper_Particle_decl.h"
#include "backend_types/Pythia_8_209/wrapper_Particle_def.h"
#include "gambit/Backends/function_return_utils.hpp"

namespace Pythia8
{
    double m2__BOSS(const Pythia8::Particle__BOSS& arg_1, const Pythia8::Particle__BOSS& arg_2)
    {
        return m2(dynamic_cast< const Pythia8::Particle& >(*arg_1.BEptr), dynamic_cast< const Pythia8::Particle& >(*arg_2.BEptr));
    }
    

}

