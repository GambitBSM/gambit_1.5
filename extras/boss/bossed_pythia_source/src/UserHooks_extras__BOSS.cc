#include "backend_types/Pythia_8_209/wrapper_Info.h"
#include "backend_types/Pythia_8_209/wrapper_Settings.h"
#include "backend_types/Pythia_8_209/wrapper_ParticleData.h"
#include "backend_types/Pythia_8_209/wrapper_Rndm.h"
#include "backend_types/Pythia_8_209/wrapper_BeamParticle.h"
#include "backend_types/Pythia_8_209/wrapper_CoupSM.h"
#include "backend_types/Pythia_8_209/wrapper_SigmaTotal.h"
#include "backend_types/Pythia_8_209/wrapper_Event.h"
#include "backend_types/Pythia_8_209/wrapper_SigmaProcess.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/UserHooks.h"

bool Pythia8::UserHooks::doVetoProcessLevel__BOSS(Pythia8::Abstract_Event& arg_1)
{
    return doVetoProcessLevel(dynamic_cast< Pythia8::Event& >(arg_1));
}


bool Pythia8::UserHooks::doVetoResonanceDecays__BOSS(Pythia8::Abstract_Event& arg_1)
{
    return doVetoResonanceDecays(dynamic_cast< Pythia8::Event& >(arg_1));
}


bool Pythia8::UserHooks::doVetoPT__BOSS(int arg_1, const Pythia8::Abstract_Event& arg_2)
{
    return doVetoPT(arg_1, dynamic_cast< const Pythia8::Event& >(arg_2));
}


bool Pythia8::UserHooks::doVetoStep__BOSS(int arg_1, int arg_2, int arg_3, const Pythia8::Abstract_Event& arg_4)
{
    return doVetoStep(arg_1, arg_2, arg_3, dynamic_cast< const Pythia8::Event& >(arg_4));
}


bool Pythia8::UserHooks::doVetoMPIStep__BOSS(int arg_1, const Pythia8::Abstract_Event& arg_2)
{
    return doVetoMPIStep(arg_1, dynamic_cast< const Pythia8::Event& >(arg_2));
}


bool Pythia8::UserHooks::doVetoPartonLevelEarly__BOSS(const Pythia8::Abstract_Event& arg_1)
{
    return doVetoPartonLevelEarly(dynamic_cast< const Pythia8::Event& >(arg_1));
}


bool Pythia8::UserHooks::doVetoPartonLevel__BOSS(const Pythia8::Abstract_Event& arg_1)
{
    return doVetoPartonLevel(dynamic_cast< const Pythia8::Event& >(arg_1));
}


double Pythia8::UserHooks::scaleResonance__BOSS(int arg_1, const Pythia8::Abstract_Event& arg_2)
{
    return scaleResonance(arg_1, dynamic_cast< const Pythia8::Event& >(arg_2));
}


bool Pythia8::UserHooks::doVetoISREmission__BOSS(int arg_1, const Pythia8::Abstract_Event& arg_2, int arg_3)
{
    return doVetoISREmission(arg_1, dynamic_cast< const Pythia8::Event& >(arg_2), arg_3);
}


bool Pythia8::UserHooks::doVetoFSREmission__BOSS(int arg_1, const Pythia8::Abstract_Event& arg_2, int arg_3, bool arg_4)
{
    return doVetoFSREmission(arg_1, dynamic_cast< const Pythia8::Event& >(arg_2), arg_3, arg_4);
}


bool Pythia8::UserHooks::doVetoFSREmission__BOSS(int arg_1, const Pythia8::Abstract_Event& arg_2, int arg_3)
{
    return doVetoFSREmission(arg_1, dynamic_cast< const Pythia8::Event& >(arg_2), arg_3);
}


bool Pythia8::UserHooks::doVetoMPIEmission__BOSS(int arg_1, const Pythia8::Abstract_Event& arg_2)
{
    return doVetoMPIEmission(arg_1, dynamic_cast< const Pythia8::Event& >(arg_2));
}


bool Pythia8::UserHooks::doReconnectResonanceSystems__BOSS(int arg_1, Pythia8::Abstract_Event& arg_2)
{
    return doReconnectResonanceSystems(arg_1, dynamic_cast< Pythia8::Event& >(arg_2));
}


void Pythia8::UserHooks::omitResonanceDecays__BOSS(const Pythia8::Abstract_Event& process, bool finalOnly)
{
    omitResonanceDecays(dynamic_cast< const Pythia8::Event& >(process), finalOnly);
}


void Pythia8::UserHooks::omitResonanceDecays__BOSS(const Pythia8::Abstract_Event& process)
{
    omitResonanceDecays(dynamic_cast< const Pythia8::Event& >(process));
}


void Pythia8::UserHooks::subEvent__BOSS(const Pythia8::Abstract_Event& event, bool isHardest)
{
    subEvent(dynamic_cast< const Pythia8::Event& >(event), isHardest);
}


void Pythia8::UserHooks::subEvent__BOSS(const Pythia8::Abstract_Event& event)
{
    subEvent(dynamic_cast< const Pythia8::Event& >(event));
}




#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_UserHooks* Pythia8::UserHooks::pointerCopy__BOSS()
{
    Pythia8::Abstract_UserHooks* new_ptr = new Pythia8::UserHooks(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::UserHooks::pointerAssign__BOSS(Pythia8::Abstract_UserHooks* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::UserHooks* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<UserHooks*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
