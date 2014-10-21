#ifndef __WRAPPERDELETER_Pythia_8_186_H__
#define __WRAPPERDELETER_Pythia_8_186_H__

#include "backend_types/Pythia_8_186/wrapper_Pythia_decl.h"
#include "backend_types/Pythia_8_186/wrapper_Event_decl.h"
#include "backend_types/Pythia_8_186/wrapper_Hist_decl.h"
#include "backend_types/Pythia_8_186/wrapper_Vec4_decl.h"
#include "backend_types/Pythia_8_186/wrapper_Info_decl.h"
#include "backend_types/Pythia_8_186/wrapper_Particle_decl.h"
#include "wrappertypedefs.h"

void wrapper_deleter(Pythia8::Particle__BOSS*);

void wrapper_deleter(Pythia8::Info__BOSS*);

void wrapper_deleter(Pythia8::Vec4__BOSS*);

void wrapper_deleter(Pythia8::Hist__BOSS*);

void wrapper_deleter(Pythia8::Event__BOSS*);

void wrapper_deleter(Pythia8::Pythia__BOSS*);

#endif /* __WRAPPERDELETER_Pythia_8_186_H__ */
