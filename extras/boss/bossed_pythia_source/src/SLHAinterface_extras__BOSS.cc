#include "backend_types/Pythia_8_209/wrapper_CoupSUSY.h"
#include "backend_types/Pythia_8_209/wrapper_Couplings.h"
#include "backend_types/Pythia_8_209/wrapper_Info.h"
#include "backend_types/Pythia_8_209/wrapper_Settings.h"
#include "backend_types/Pythia_8_209/wrapper_Rndm.h"
#include "backend_types/Pythia_8_209/wrapper_ParticleData.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/SLHAinterface.h"

void Pythia8::SLHAinterface::setPtr__BOSS(Pythia8::Abstract_Info* infoPtrIn)
{
    setPtr(dynamic_cast< Pythia8::Info* >(infoPtrIn));
}


void Pythia8::SLHAinterface::init__BOSS(Pythia8::Abstract_Settings& settings, Pythia8::Abstract_Rndm* rndmPtr, Pythia8::Abstract_Couplings* couplingsPtrIn, Pythia8::Abstract_ParticleData* particleDataPtr, bool& useSHLAcouplings, std::basic_stringstream<char,std::char_traits<char>,std::allocator<char> >& ParticleDataBuffer)
{
    init(dynamic_cast< Pythia8::Settings& >(settings), dynamic_cast< Pythia8::Rndm* >(rndmPtr), dynamic_cast< Pythia8::Couplings* >(couplingsPtrIn), dynamic_cast< Pythia8::ParticleData* >(particleDataPtr), useSHLAcouplings, ParticleDataBuffer);
}


bool Pythia8::SLHAinterface::initSLHA__BOSS(Pythia8::Abstract_Settings& settings, Pythia8::Abstract_ParticleData* particleDataPtr)
{
    return initSLHA(dynamic_cast< Pythia8::Settings& >(settings), dynamic_cast< Pythia8::ParticleData* >(particleDataPtr));
}


void Pythia8::SLHAinterface::pythia2slha__BOSS(Pythia8::Abstract_ParticleData* particleDataPtr)
{
    pythia2slha(dynamic_cast< Pythia8::ParticleData* >(particleDataPtr));
}



Pythia8::Abstract_CoupSUSY& Pythia8::SLHAinterface::coupSUSY_ref__BOSS() { return coupSUSY; }

int& Pythia8::SLHAinterface::meMode_ref__BOSS() { return meMode; }


#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_SLHAinterface* Pythia8::SLHAinterface::pointerCopy__BOSS()
{
    Pythia8::Abstract_SLHAinterface* new_ptr = new Pythia8::SLHAinterface(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::SLHAinterface::pointerAssign__BOSS(Pythia8::Abstract_SLHAinterface* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::SLHAinterface* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<SLHAinterface*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
