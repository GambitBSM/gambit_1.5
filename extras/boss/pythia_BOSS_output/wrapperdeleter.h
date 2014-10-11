#ifndef __WRAPPERDELETER_BOSSedPythia_1_0_H__
#define __WRAPPERDELETER_BOSSedPythia_1_0_H__

#include "backend_types/BOSSedPythia_1_0/wrapper_Pythia_decl.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Event_decl.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Hist_decl.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Vec4_decl.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Info_decl.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Particle_decl.h"
#include "wrappers_typedefs.hpp"

void wrapper_deleter(Pythia8::Particle_GAMBIT*);

void wrapper_deleter(Pythia8::Info_GAMBIT*);

void wrapper_deleter(Pythia8::Vec4_GAMBIT*);

void wrapper_deleter(Pythia8::Hist_GAMBIT*);

void wrapper_deleter(Pythia8::Event_GAMBIT*);

void wrapper_deleter(Pythia8::Pythia_GAMBIT*);

#endif /* __WRAPPERDELETER_BOSSedPythia_1_0_H__ */
