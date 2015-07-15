#include "backend_types/Pythia_8_209/wrapper_Info.h"
#include "backend_types/Pythia_8_209/wrapper_ParticleData.h"
#include "backend_types/Pythia_8_209/wrapper_Settings.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/SigmaTotal.h"

void Pythia8::SigmaTotal::init__BOSS(Pythia8::Abstract_Info* infoPtrIn, Pythia8::Abstract_Settings& settings, Pythia8::Abstract_ParticleData* particleDataPtrIn)
{
    init(dynamic_cast< Pythia8::Info* >(infoPtrIn), dynamic_cast< Pythia8::Settings& >(settings), dynamic_cast< Pythia8::ParticleData* >(particleDataPtrIn));
}




#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_SigmaTotal* Pythia8::SigmaTotal::pointerCopy__BOSS()
{
    Pythia8::Abstract_SigmaTotal* new_ptr = new Pythia8::SigmaTotal(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::SigmaTotal::pointerAssign__BOSS(Pythia8::Abstract_SigmaTotal* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::SigmaTotal* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<SigmaTotal*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
