#ifndef __wrapperdeleter_Pythia_8_186_h__
#define __wrapperdeleter_Pythia_8_186_h__

#include "backend_types/Pythia_8_186/wrapper_CoupSM.h"
#include "backend_types/Pythia_8_186/wrapper_Pythia.h"
#include "backend_types/Pythia_8_186/wrapper_Event.h"
#include "backend_types/Pythia_8_186/wrapper_Couplings.h"
#include "backend_types/Pythia_8_186/wrapper_SLHAinterface.h"
#include "backend_types/Pythia_8_186/wrapper_Hist.h"
#include "backend_types/Pythia_8_186/wrapper_DecayChannel.h"
#include "backend_types/Pythia_8_186/wrapper_Vec4.h"
#include "backend_types/Pythia_8_186/wrapper_Settings.h"
#include "backend_types/Pythia_8_186/wrapper_ParticleDataEntry.h"
#include "backend_types/Pythia_8_186/wrapper_Info.h"
#include "backend_types/Pythia_8_186/wrapper_SlowJet.h"
#include "backend_types/Pythia_8_186/wrapper_CoupSUSY.h"
#include "backend_types/Pythia_8_186/wrapper_ParticleData.h"
#include "backend_types/Pythia_8_186/wrapper_Particle.h"
#include "backend_types/Pythia_8_186/wrapper_PartonLevel.h"
#include "backend_types/Pythia_8_186/wrapper_Rndm.h"
#include "gambit/Backends/wrappertypedefs.h"

void wrapper_deleter(Pythia8::Rndm__BOSS*);

void wrapper_deleter(Pythia8::PartonLevel__BOSS*);

void wrapper_deleter(Pythia8::Particle__BOSS*);

void wrapper_deleter(Pythia8::ParticleData__BOSS*);

void wrapper_deleter(Pythia8::CoupSUSY__BOSS*);

void wrapper_deleter(Pythia8::SlowJet__BOSS*);

void wrapper_deleter(Pythia8::Info__BOSS*);

void wrapper_deleter(Pythia8::ParticleDataEntry__BOSS*);

void wrapper_deleter(Pythia8::Settings__BOSS*);

void wrapper_deleter(Pythia8::Vec4__BOSS*);

void wrapper_deleter(Pythia8::DecayChannel__BOSS*);

void wrapper_deleter(Pythia8::Hist__BOSS*);

void wrapper_deleter(Pythia8::SLHAinterface__BOSS*);

void wrapper_deleter(Pythia8::Couplings__BOSS*);

void wrapper_deleter(Pythia8::Event__BOSS*);

void wrapper_deleter(Pythia8::Pythia__BOSS*);

void wrapper_deleter(Pythia8::CoupSM__BOSS*);

#endif /* __wrapperdeleter_Pythia_8_186_h__ */
