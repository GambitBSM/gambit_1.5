#include <vector>
#include "backend_types/Pythia_8_186/wrapper_Info.h"
#include "backend_types/Pythia_8_186/wrapper_ParticleData.h"
#include "backend_types/Pythia_8_186/wrapper_Rndm.h"
#include "backend_types/Pythia_8_186/wrapper_Couplings.h"
#include "backend_types/Pythia_8_186/wrapper_ParticleDataEntry.h"
#include "backend_types/Pythia_8_186/wrapper_Settings.h"
#include "backend_types/Pythia_8_186/wrapper_Event.h"
#include "backend_types/Pythia_8_186/wrapper_Particle.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/ParticleDecays.h"

bool Pythia8::ParticleDecays::decay__BOSS(int iDec, Pythia8::Abstract_Event& event)
{
    return decay(iDec, dynamic_cast< Pythia8::Event& >(event));
}


bool Pythia8::ParticleDecays::checkVertex__BOSS(Pythia8::Abstract_Particle& decayer)
{
    return checkVertex(dynamic_cast< Pythia8::Particle& >(decayer));
}


bool Pythia8::ParticleDecays::oscillateB__BOSS(Pythia8::Abstract_Particle& decayer)
{
    return oscillateB(dynamic_cast< Pythia8::Particle& >(decayer));
}


bool Pythia8::ParticleDecays::oneBody__BOSS(Pythia8::Abstract_Event& event)
{
    return oneBody(dynamic_cast< Pythia8::Event& >(event));
}


bool Pythia8::ParticleDecays::twoBody__BOSS(Pythia8::Abstract_Event& event)
{
    return twoBody(dynamic_cast< Pythia8::Event& >(event));
}


bool Pythia8::ParticleDecays::threeBody__BOSS(Pythia8::Abstract_Event& event)
{
    return threeBody(dynamic_cast< Pythia8::Event& >(event));
}


bool Pythia8::ParticleDecays::mGenerator__BOSS(Pythia8::Abstract_Event& event)
{
    return mGenerator(dynamic_cast< Pythia8::Event& >(event));
}


bool Pythia8::ParticleDecays::dalitzKinematics__BOSS(Pythia8::Abstract_Event& event)
{
    return dalitzKinematics(dynamic_cast< Pythia8::Event& >(event));
}


bool Pythia8::ParticleDecays::setColours__BOSS(Pythia8::Abstract_Event& event)
{
    return setColours(dynamic_cast< Pythia8::Event& >(event));
}




#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_ParticleDecays* Pythia8::ParticleDecays::pointerCopy__BOSS()
{
    Pythia8::Abstract_ParticleDecays* new_ptr = new Pythia8::ParticleDecays(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::ParticleDecays::pointerAssign__BOSS(Pythia8::Abstract_ParticleDecays* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::ParticleDecays* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<ParticleDecays*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
