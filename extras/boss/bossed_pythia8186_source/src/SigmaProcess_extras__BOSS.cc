#include <string>
#include "backend_types/Pythia_8_186/wrapper_Info.h"
#include "backend_types/Pythia_8_186/wrapper_Settings.h"
#include "backend_types/Pythia_8_186/wrapper_ParticleData.h"
#include "backend_types/Pythia_8_186/wrapper_Rndm.h"
#include "backend_types/Pythia_8_186/wrapper_BeamParticle.h"
#include "backend_types/Pythia_8_186/wrapper_Couplings.h"
#include "backend_types/Pythia_8_186/wrapper_SigmaTotal.h"
#include "backend_types/Pythia_8_186/wrapper_Particle.h"
#include "backend_types/Pythia_8_186/wrapper_Vec4.h"
#include "backend_types/Pythia_8_186/wrapper_SLHAinterface.h"
#include "backend_types/Pythia_8_186/wrapper_Event.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/SigmaProcess.h"

void Pythia8::SigmaProcess::init__BOSS(Pythia8::Abstract_Info* infoPtrIn, Pythia8::Abstract_Settings* settingsPtrIn, Pythia8::Abstract_ParticleData* particleDataPtrIn, Pythia8::Abstract_Rndm* rndmPtrIn, Pythia8::Abstract_BeamParticle* beamAPtrIn, Pythia8::Abstract_BeamParticle* beamBPtrIn, Pythia8::Abstract_Couplings* couplings, Pythia8::Abstract_SigmaTotal* sigmaTotPtrIn, Pythia8::Abstract_SLHAinterface* slhaInterfacePtrIn)
{
    init(dynamic_cast< Pythia8::Info* >(infoPtrIn), dynamic_cast< Pythia8::Settings* >(settingsPtrIn), dynamic_cast< Pythia8::ParticleData* >(particleDataPtrIn), dynamic_cast< Pythia8::Rndm* >(rndmPtrIn), dynamic_cast< Pythia8::BeamParticle* >(beamAPtrIn), dynamic_cast< Pythia8::BeamParticle* >(beamBPtrIn), dynamic_cast< Pythia8::Couplings* >(couplings), dynamic_cast< Pythia8::SigmaTotal* >(sigmaTotPtrIn), dynamic_cast< Pythia8::SLHAinterface* >(slhaInterfacePtrIn));
}


void Pythia8::SigmaProcess::init__BOSS(Pythia8::Abstract_Info* infoPtrIn, Pythia8::Abstract_Settings* settingsPtrIn, Pythia8::Abstract_ParticleData* particleDataPtrIn, Pythia8::Abstract_Rndm* rndmPtrIn, Pythia8::Abstract_BeamParticle* beamAPtrIn, Pythia8::Abstract_BeamParticle* beamBPtrIn, Pythia8::Abstract_Couplings* couplings, Pythia8::Abstract_SigmaTotal* sigmaTotPtrIn)
{
    init(dynamic_cast< Pythia8::Info* >(infoPtrIn), dynamic_cast< Pythia8::Settings* >(settingsPtrIn), dynamic_cast< Pythia8::ParticleData* >(particleDataPtrIn), dynamic_cast< Pythia8::Rndm* >(rndmPtrIn), dynamic_cast< Pythia8::BeamParticle* >(beamAPtrIn), dynamic_cast< Pythia8::BeamParticle* >(beamBPtrIn), dynamic_cast< Pythia8::Couplings* >(couplings), dynamic_cast< Pythia8::SigmaTotal* >(sigmaTotPtrIn));
}


void Pythia8::SigmaProcess::init__BOSS(Pythia8::Abstract_Info* infoPtrIn, Pythia8::Abstract_Settings* settingsPtrIn, Pythia8::Abstract_ParticleData* particleDataPtrIn, Pythia8::Abstract_Rndm* rndmPtrIn, Pythia8::Abstract_BeamParticle* beamAPtrIn, Pythia8::Abstract_BeamParticle* beamBPtrIn, Pythia8::Abstract_Couplings* couplings)
{
    init(dynamic_cast< Pythia8::Info* >(infoPtrIn), dynamic_cast< Pythia8::Settings* >(settingsPtrIn), dynamic_cast< Pythia8::ParticleData* >(particleDataPtrIn), dynamic_cast< Pythia8::Rndm* >(rndmPtrIn), dynamic_cast< Pythia8::BeamParticle* >(beamAPtrIn), dynamic_cast< Pythia8::BeamParticle* >(beamBPtrIn), dynamic_cast< Pythia8::Couplings* >(couplings));
}


void Pythia8::SigmaProcess::set3Kin__BOSS(double arg_1, double arg_2, double arg_3, Pythia8::Abstract_Vec4& arg_4, Pythia8::Abstract_Vec4& arg_5, Pythia8::Abstract_Vec4& arg_6, double arg_7, double arg_8, double arg_9, double arg_10, double arg_11, double arg_12)
{
    set3Kin(arg_1, arg_2, arg_3, dynamic_cast< Pythia8::Vec4& >(arg_4), dynamic_cast< Pythia8::Vec4& >(arg_5), dynamic_cast< Pythia8::Vec4& >(arg_6), arg_7, arg_8, arg_9, arg_10, arg_11, arg_12);
}


double Pythia8::SigmaProcess::sigmaHatWrap__BOSS(int id1in)
{
    return sigmaHatWrap(id1in);
}


double Pythia8::SigmaProcess::sigmaHatWrap__BOSS()
{
    return sigmaHatWrap();
}


void Pythia8::SigmaProcess::pickInState__BOSS(int id1in)
{
    pickInState(id1in);
}


void Pythia8::SigmaProcess::pickInState__BOSS()
{
    pickInState();
}


bool Pythia8::SigmaProcess::final2KinMPI__BOSS(int arg_1, int arg_2, Pythia8::Abstract_Vec4& arg_3, Pythia8::Abstract_Vec4& arg_4, double arg_5, double arg_6)
{
    return final2KinMPI(arg_1, arg_2, dynamic_cast< Pythia8::Vec4& >(arg_3), dynamic_cast< Pythia8::Vec4& >(arg_4), arg_5, arg_6);
}


bool Pythia8::SigmaProcess::final2KinMPI__BOSS(int arg_1, int arg_2, Pythia8::Abstract_Vec4& arg_3, Pythia8::Abstract_Vec4& arg_4, double arg_5)
{
    return final2KinMPI(arg_1, arg_2, dynamic_cast< Pythia8::Vec4& >(arg_3), dynamic_cast< Pythia8::Vec4& >(arg_4), arg_5);
}


bool Pythia8::SigmaProcess::final2KinMPI__BOSS(int arg_1, int arg_2, Pythia8::Abstract_Vec4& arg_3, Pythia8::Abstract_Vec4& arg_4)
{
    return final2KinMPI(arg_1, arg_2, dynamic_cast< Pythia8::Vec4& >(arg_3), dynamic_cast< Pythia8::Vec4& >(arg_4));
}


bool Pythia8::SigmaProcess::final2KinMPI__BOSS(int arg_1, int arg_2, Pythia8::Abstract_Vec4& arg_3)
{
    return final2KinMPI(arg_1, arg_2, dynamic_cast< Pythia8::Vec4& >(arg_3));
}


bool Pythia8::SigmaProcess::final2KinMPI__BOSS(int arg_1, int arg_2)
{
    return final2KinMPI(arg_1, arg_2);
}


bool Pythia8::SigmaProcess::final2KinMPI__BOSS(int arg_1)
{
    return final2KinMPI(arg_1);
}


bool Pythia8::SigmaProcess::final2KinMPI__BOSS()
{
    return final2KinMPI();
}


double Pythia8::SigmaProcess::weightDecayFlav__BOSS(Pythia8::Abstract_Event& arg_1)
{
    return weightDecayFlav(dynamic_cast< Pythia8::Event& >(arg_1));
}


double Pythia8::SigmaProcess::weightDecay__BOSS(Pythia8::Abstract_Event& arg_1, int arg_2, int arg_3)
{
    return weightDecay(dynamic_cast< Pythia8::Event& >(arg_1), arg_2, arg_3);
}


Pythia8::Abstract_Particle* Pythia8::SigmaProcess::getParton__BOSS(int i) const
{
    return new Pythia8::Particle(getParton(i));
}


void Pythia8::SigmaProcess::setId__BOSS(int id1in, int id2in, int id3in, int id4in)
{
    setId(id1in, id2in, id3in, id4in);
}


void Pythia8::SigmaProcess::setId__BOSS(int id1in, int id2in, int id3in)
{
    setId(id1in, id2in, id3in);
}


void Pythia8::SigmaProcess::setId__BOSS(int id1in, int id2in)
{
    setId(id1in, id2in);
}


void Pythia8::SigmaProcess::setId__BOSS(int id1in)
{
    setId(id1in);
}


void Pythia8::SigmaProcess::setId__BOSS()
{
    setId();
}


void Pythia8::SigmaProcess::setColAcol__BOSS(int col1, int acol1, int col2, int acol2, int col3, int acol3, int col4, int acol4, int col5)
{
    setColAcol(col1, acol1, col2, acol2, col3, acol3, col4, acol4, col5);
}


void Pythia8::SigmaProcess::setColAcol__BOSS(int col1, int acol1, int col2, int acol2, int col3, int acol3, int col4, int acol4)
{
    setColAcol(col1, acol1, col2, acol2, col3, acol3, col4, acol4);
}


void Pythia8::SigmaProcess::setColAcol__BOSS(int col1, int acol1, int col2, int acol2, int col3, int acol3, int col4)
{
    setColAcol(col1, acol1, col2, acol2, col3, acol3, col4);
}


void Pythia8::SigmaProcess::setColAcol__BOSS(int col1, int acol1, int col2, int acol2, int col3, int acol3)
{
    setColAcol(col1, acol1, col2, acol2, col3, acol3);
}


void Pythia8::SigmaProcess::setColAcol__BOSS(int col1, int acol1, int col2, int acol2, int col3)
{
    setColAcol(col1, acol1, col2, acol2, col3);
}


void Pythia8::SigmaProcess::setColAcol__BOSS(int col1, int acol1, int col2, int acol2)
{
    setColAcol(col1, acol1, col2, acol2);
}


void Pythia8::SigmaProcess::setColAcol__BOSS(int col1, int acol1, int col2)
{
    setColAcol(col1, acol1, col2);
}


void Pythia8::SigmaProcess::setColAcol__BOSS(int col1, int acol1)
{
    setColAcol(col1, acol1);
}


void Pythia8::SigmaProcess::setColAcol__BOSS(int col1)
{
    setColAcol(col1);
}


void Pythia8::SigmaProcess::setColAcol__BOSS()
{
    setColAcol();
}


double Pythia8::SigmaProcess::weightTopDecay__BOSS(Pythia8::Abstract_Event& process, int iResBeg, int iResEnd)
{
    return weightTopDecay(dynamic_cast< Pythia8::Event& >(process), iResBeg, iResEnd);
}


double Pythia8::SigmaProcess::weightHiggsDecay__BOSS(Pythia8::Abstract_Event& process, int iResBeg, int iResEnd)
{
    return weightHiggsDecay(dynamic_cast< Pythia8::Event& >(process), iResBeg, iResEnd);
}




#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_SigmaProcess* Pythia8::SigmaProcess::pointerCopy__BOSS()
{
    Pythia8::Abstract_SigmaProcess* new_ptr = new Pythia8::SigmaProcess(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::SigmaProcess::pointerAssign__BOSS(Pythia8::Abstract_SigmaProcess* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::SigmaProcess* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<SigmaProcess*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
