#ifndef __boss__SLHAinterface_Pythia_8_209_h__
#define __boss__SLHAinterface_Pythia_8_209_h__

// SLHAinterface.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// Main authors of this file: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for SUSY Les Houches Accord Interface.
// Handles the communication between PYTHIA and the SusyLesHouches classes.

#ifndef Pythia8_SLHAinterface_H
#define Pythia8_SLHAinterface_H

#include "Pythia8/Basics.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/SusyCouplings.h"
#include "Pythia8/SusyLesHouches.h"

namespace Pythia8 {

//==========================================================================

// The SLHAinterface class handles communication between Pythia and
// SusyLesHouches.

} 
#define ENUMS_DECLARED
#include "backend_types/Pythia_8_209/abstract_SLHAinterface.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
namespace Pythia8 { 
class SLHAinterface : public virtual Abstract_SLHAinterface {

public:

  // Constructor.
  SLHAinterface() {} ;

  // Set pointers
  void setPtr( Info* infoPtrIn ) {infoPtr     = infoPtrIn;}

  // Initialize and switch to SUSY couplings if reading SLHA spectrum
  void init( Settings& settings, Rndm* rndmPtr, Couplings* couplingsPtrIn,
    ParticleData* particleDataPtr, bool& useSHLAcouplings,
    stringstream& ParticleDataBuffer );

  // Initialize SUSY Les Houches Accord data.
  bool initSLHA(Settings& settings, ParticleData* particleDataPtr);

  // Initialize SLHA blocks SMINPUTS and MASS from PYTHIA SM parameter values.
  // E.g., to make sure that there are no important unfilled entries
  void pythia2slha(ParticleData* particleDataPtr);

  // SusyLesHouches - SLHA object for interface to SUSY spectra.
  SusyLesHouches slha;

  // SLHA derived couplings class and pointer to Couplings object
  CoupSUSY       coupSUSY;
  Couplings*     couplingsPtr;

  // Pointers to PYTHIA objects
  Info*          infoPtr;
  Settings*      settingsPtr;

  // Internal data members
  int            meMode;


        public:
            Abstract_SLHAinterface* pointerCopy__BOSS();

            void pointerAssign__BOSS(Abstract_SLHAinterface* in);

        public:
            Pythia8::Abstract_SusyLesHouches& slha_ref__BOSS();

            Pythia8::Abstract_CoupSUSY& coupSUSY_ref__BOSS();

            int& meMode_ref__BOSS();



        public:
            void setPtr__BOSS(Pythia8::Abstract_Info*);

            void init__BOSS(Pythia8::Abstract_Settings&, Pythia8::Abstract_Rndm*, Pythia8::Abstract_Couplings*, Pythia8::Abstract_ParticleData*, bool&, std::basic_stringstream<char,std::char_traits<char>,std::allocator<char> >&);

            bool initSLHA__BOSS(Pythia8::Abstract_Settings&, Pythia8::Abstract_ParticleData*);

            void pythia2slha__BOSS(Pythia8::Abstract_ParticleData*);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SLHAinterface_H

#endif /* __boss__SLHAinterface_Pythia_8_209_h__ */
