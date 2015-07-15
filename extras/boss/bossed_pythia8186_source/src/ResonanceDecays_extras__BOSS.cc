#include <vector>
#include "backend_types/Pythia_8_186/wrapper_Info.h"
#include "backend_types/Pythia_8_186/wrapper_ParticleData.h"
#include "backend_types/Pythia_8_186/wrapper_Rndm.h"
#include "backend_types/Pythia_8_186/wrapper_Event.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/ResonanceDecays.h"

void Pythia8::ResonanceDecays::init__BOSS(Pythia8::Abstract_Info* infoPtrIn, Pythia8::Abstract_ParticleData* particleDataPtrIn, Pythia8::Abstract_Rndm* rndmPtrIn)
{
    init(dynamic_cast< Pythia8::Info* >(infoPtrIn), dynamic_cast< Pythia8::ParticleData* >(particleDataPtrIn), dynamic_cast< Pythia8::Rndm* >(rndmPtrIn));
}


bool Pythia8::ResonanceDecays::next__BOSS(Pythia8::Abstract_Event& process, int iDecNow)
{
    return next(dynamic_cast< Pythia8::Event& >(process), iDecNow);
}


bool Pythia8::ResonanceDecays::next__BOSS(Pythia8::Abstract_Event& process)
{
    return next(dynamic_cast< Pythia8::Event& >(process));
}


bool Pythia8::ResonanceDecays::pickColours__BOSS(int iDec, Pythia8::Abstract_Event& process)
{
    return pickColours(iDec, dynamic_cast< Pythia8::Event& >(process));
}




#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_ResonanceDecays* Pythia8::ResonanceDecays::pointerCopy__BOSS()
{
    Pythia8::Abstract_ResonanceDecays* new_ptr = new Pythia8::ResonanceDecays(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::ResonanceDecays::pointerAssign__BOSS(Pythia8::Abstract_ResonanceDecays* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::ResonanceDecays* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<ResonanceDecays*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
