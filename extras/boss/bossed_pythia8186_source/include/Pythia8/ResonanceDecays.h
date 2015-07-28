#ifndef __boss__ResonanceDecays_Pythia_8_186_h__
#define __boss__ResonanceDecays_Pythia_8_186_h__

// ResonanceDecays.h is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main class for performing resonance decays.
// ResonanceDecays: handles the sequential decay of resonances in process.

#ifndef Pythia8_ResonanceDecays_H
#define Pythia8_ResonanceDecays_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/ResonanceWidths.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {
  
//==========================================================================

// The ResonanceDecays class handles the sequential decay of resonances
// that are part of the hard process (t, W, Z, H, SUSY,...).

} 
#define ENUMS_DECLARED
#include "backend_types/Pythia_8_186/abstract_ResonanceDecays.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
namespace Pythia8 { 
class ResonanceDecays : public virtual Abstract_ResonanceDecays {

public:

  // Constructor.
  ResonanceDecays() {}

  // Store pointers to Info and Rndm for error messages and random numbers.
  void init(Info* infoPtrIn,  ParticleData* particleDataPtrIn,
    Rndm* rndmPtrIn) {infoPtr = infoPtrIn;
    particleDataPtr = particleDataPtrIn; rndmPtr = rndmPtrIn;}
 
  // Generate the next decay sequence.
  bool next( Event& process, int iDecNow = 0);

private:

  // Constants: could only be changed in the code itself.
  static const int    NTRYCHANNEL, NTRYMASSES;
  static const double MSAFETY, WIDTHCUT, TINY, TINYBWRANGE,
                      WTCORRECTION[11];

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Rndm*         rndmPtr;

  // Select masses of decay products.
  bool pickMasses();

  // Select colours of decay products.
  bool pickColours(int iDec, Event& process);

  // Select kinematics isotropic in phase space.
  bool pickKinematics();

  // Flavour, colour and momentum information.
  int            id0, mult;
  double         m0;
  vector<int>    idProd, cols, acols;
  vector<double> mProd;
  vector<Vec4>   pProd;


        public:
            Abstract_ResonanceDecays* pointerCopy__BOSS();

            void pointerAssign__BOSS(Abstract_ResonanceDecays* in);


        public:
            void init__BOSS(Pythia8::Abstract_Info*, Pythia8::Abstract_ParticleData*, Pythia8::Abstract_Rndm*);

            bool next__BOSS(Pythia8::Abstract_Event&, int);

            bool next__BOSS(Pythia8::Abstract_Event&);

        private:
            bool pickColours__BOSS(int, Pythia8::Abstract_Event&);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_ResonanceDecays_H

#endif /* __boss__ResonanceDecays_Pythia_8_186_h__ */
