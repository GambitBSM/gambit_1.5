#ifndef __wrapperdeleter_Pythia_8_209_h__
#define __wrapperdeleter_Pythia_8_209_h__

#include "backend_types/Pythia_8_209/wrapper_ParticleDataEntry.h"
#include "backend_types/Pythia_8_209/wrapper_AlphaStrong.h"
#include "backend_types/Pythia_8_209/wrapper_Pythia.h"
#include "backend_types/Pythia_8_209/wrapper_Event.h"
#include "backend_types/Pythia_8_209/wrapper_Couplings.h"
#include "backend_types/Pythia_8_209/wrapper_AlphaEM.h"
#include "backend_types/Pythia_8_209/wrapper_LHdecayChannel.h"
#include "backend_types/Pythia_8_209/wrapper_ResonanceGmZ.h"
#include "backend_types/Pythia_8_209/wrapper_BeamParticle.h"
#include "backend_types/Pythia_8_209/wrapper_CoupSM.h"
#include "backend_types/Pythia_8_209/wrapper_SigmaTotal.h"
#include "backend_types/Pythia_8_209/wrapper_ParticleDecays.h"
#include "backend_types/Pythia_8_209/wrapper_Particle.h"
#include "backend_types/Pythia_8_209/wrapper_ResonanceWidths.h"
#include "backend_types/Pythia_8_209/wrapper_ResonanceDecays.h"
#include "backend_types/Pythia_8_209/wrapper_PartonLevel.h"
#include "backend_types/Pythia_8_209/wrapper_Rndm.h"
#include "backend_types/Pythia_8_209/wrapper_UserHooks.h"
#include "backend_types/Pythia_8_209/wrapper_Parm.h"
#include "backend_types/Pythia_8_209/wrapper_SigmaProcess.h"
#include "backend_types/Pythia_8_209/wrapper_LHdecayTable.h"
#include "backend_types/Pythia_8_209/wrapper_SusyLesHouches.h"
#include "backend_types/Pythia_8_209/wrapper_SLHAinterface.h"
#include "backend_types/Pythia_8_209/wrapper_SlowJet.h"
#include "backend_types/Pythia_8_209/wrapper_Hist.h"
#include "backend_types/Pythia_8_209/wrapper_Vec4.h"
#include "backend_types/Pythia_8_209/wrapper_Settings.h"
#include "backend_types/Pythia_8_209/wrapper_CoupSUSY.h"
#include "backend_types/Pythia_8_209/wrapper_DecayChannel.h"
#include "backend_types/Pythia_8_209/wrapper_Info.h"
#include "backend_types/Pythia_8_209/wrapper_ParticleData.h"
#include "gambit/Backends/wrappertypedefs.h"

void wrapper_deleter(Pythia8::ParticleData__BOSS*);

void wrapper_deleter(Pythia8::Info__BOSS*);

void wrapper_deleter(Pythia8::DecayChannel__BOSS*);

void wrapper_deleter(Pythia8::CoupSUSY__BOSS*);

void wrapper_deleter(Pythia8::Settings__BOSS*);

void wrapper_deleter(Pythia8::Vec4__BOSS*);

void wrapper_deleter(Pythia8::Hist__BOSS*);

void wrapper_deleter(Pythia8::SlowJet__BOSS*);

void wrapper_deleter(Pythia8::SLHAinterface__BOSS*);

void wrapper_deleter(Pythia8::SusyLesHouches__BOSS*);

void wrapper_deleter(Pythia8::LHdecayTable__BOSS*);

void wrapper_deleter(Pythia8::SigmaProcess__BOSS*);

void wrapper_deleter(Pythia8::Parm__BOSS*);

void wrapper_deleter(Pythia8::UserHooks__BOSS*);

void wrapper_deleter(Pythia8::Rndm__BOSS*);

void wrapper_deleter(Pythia8::PartonLevel__BOSS*);

void wrapper_deleter(Pythia8::ResonanceDecays__BOSS*);

void wrapper_deleter(Pythia8::ResonanceWidths__BOSS*);

void wrapper_deleter(Pythia8::Particle__BOSS*);

void wrapper_deleter(Pythia8::ParticleDecays__BOSS*);

void wrapper_deleter(Pythia8::SigmaTotal__BOSS*);

void wrapper_deleter(Pythia8::CoupSM__BOSS*);

void wrapper_deleter(Pythia8::BeamParticle__BOSS*);

void wrapper_deleter(Pythia8::ResonanceGmZ__BOSS*);

void wrapper_deleter(Pythia8::LHdecayChannel__BOSS*);

void wrapper_deleter(Pythia8::AlphaEM__BOSS*);

void wrapper_deleter(Pythia8::Couplings__BOSS*);

void wrapper_deleter(Pythia8::Event__BOSS*);

void wrapper_deleter(Pythia8::Pythia__BOSS*);

void wrapper_deleter(Pythia8::AlphaStrong__BOSS*);

void wrapper_deleter(Pythia8::ParticleDataEntry__BOSS*);

#endif /* __wrapperdeleter_Pythia_8_209_h__ */
