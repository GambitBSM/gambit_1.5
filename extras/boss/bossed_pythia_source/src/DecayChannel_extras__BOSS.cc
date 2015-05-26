#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/ParticleData.h"

void Pythia8::DecayChannel::bRatio__BOSS(double bRatioIn)
{
    bRatio(bRatioIn);
}




#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_DecayChannel* Pythia8::DecayChannel::pointerCopy__BOSS()
{
    Pythia8::Abstract_DecayChannel* new_ptr = new Pythia8::DecayChannel(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::DecayChannel::pointerAssign__BOSS(Pythia8::Abstract_DecayChannel* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::DecayChannel* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<DecayChannel*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
