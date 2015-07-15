#include <vector>
#include "backend_types/Pythia_8_186/wrapper_Info.h"
#include "backend_types/Pythia_8_186/wrapper_ParticleData.h"
#include "backend_types/Pythia_8_186/wrapper_Rndm.h"
#include "backend_types/Pythia_8_186/wrapper_BeamParticle.h"
#include "backend_types/Pythia_8_186/wrapper_Couplings.h"
#include "backend_types/Pythia_8_186/wrapper_UserHooks.h"
#include "backend_types/Pythia_8_186/wrapper_ResonanceDecays.h"
#include "backend_types/Pythia_8_186/wrapper_Settings.h"
#include "backend_types/Pythia_8_186/wrapper_SigmaTotal.h"
#include "backend_types/Pythia_8_186/wrapper_Event.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/PartonLevel.h"

bool Pythia8::PartonLevel::next__BOSS(Pythia8::Abstract_Event& process, Pythia8::Abstract_Event& event)
{
    return next(dynamic_cast< Pythia8::Event& >(process), dynamic_cast< Pythia8::Event& >(event));
}


void Pythia8::PartonLevel::setupShowerSys__BOSS(Pythia8::Abstract_Event& process, Pythia8::Abstract_Event& event)
{
    setupShowerSys(dynamic_cast< Pythia8::Event& >(process), dynamic_cast< Pythia8::Event& >(event));
}


bool Pythia8::PartonLevel::resonanceShowers__BOSS(Pythia8::Abstract_Event& process, Pythia8::Abstract_Event& event, bool skipForR)
{
    return resonanceShowers(dynamic_cast< Pythia8::Event& >(process), dynamic_cast< Pythia8::Event& >(event), skipForR);
}


bool Pythia8::PartonLevel::wzDecayShowers__BOSS(Pythia8::Abstract_Event& event)
{
    return wzDecayShowers(dynamic_cast< Pythia8::Event& >(event));
}


void Pythia8::PartonLevel::statistics__BOSS()
{
    statistics();
}


int Pythia8::PartonLevel::decideResolvedDiff__BOSS(Pythia8::Abstract_Event& process)
{
    return decideResolvedDiff(dynamic_cast< Pythia8::Event& >(process));
}


bool Pythia8::PartonLevel::setupUnresolvedSys__BOSS(Pythia8::Abstract_Event& process, Pythia8::Abstract_Event& event)
{
    return setupUnresolvedSys(dynamic_cast< Pythia8::Event& >(process), dynamic_cast< Pythia8::Event& >(event));
}


void Pythia8::PartonLevel::setupHardSys__BOSS(Pythia8::Abstract_Event& process, Pythia8::Abstract_Event& event)
{
    setupHardSys(dynamic_cast< Pythia8::Event& >(process), dynamic_cast< Pythia8::Event& >(event));
}


void Pythia8::PartonLevel::setupResolvedDiff__BOSS(Pythia8::Abstract_Event& process)
{
    setupResolvedDiff(dynamic_cast< Pythia8::Event& >(process));
}


void Pythia8::PartonLevel::leaveResolvedDiff__BOSS(int iHardLoop, Pythia8::Abstract_Event& process, Pythia8::Abstract_Event& event)
{
    leaveResolvedDiff(iHardLoop, dynamic_cast< Pythia8::Event& >(process), dynamic_cast< Pythia8::Event& >(event));
}




#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_PartonLevel* Pythia8::PartonLevel::pointerCopy__BOSS()
{
    Pythia8::Abstract_PartonLevel* new_ptr = new Pythia8::PartonLevel(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::PartonLevel::pointerAssign__BOSS(Pythia8::Abstract_PartonLevel* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::PartonLevel* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<PartonLevel*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
