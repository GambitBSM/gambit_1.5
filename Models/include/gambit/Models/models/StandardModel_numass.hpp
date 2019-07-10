//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Standard Model parameters, defined in SLHA2
///  conventions. Corresponds directly to the
///  SLHA2 "SMINPUTS" block.
///
///  * Common models for massive neutrinos in Cosmology *
///
///  SMINPUTS is also a CAPABILITY that module
///  writers may be interested in using.
///
///  SLHA2 conventions:
///
///  SMINPUTS block description:
///
///    // SLHA1
///    alphainv;  // 1: Inverse electromagnetic coupling at the Z pole in the MSbar scheme (with 5 active flavours)
///    GF;        // 2: Fermi constant (in units of GeV^-2)
///    alphaS;    // 3: Strong coupling at the Z pole in the MSbar scheme (with 5 active flavours).
///    mZ;        // 4: Z pole mass
///    mBmB;      // 5: b quark running mass in the MSbar scheme (at mB)
///    mT;        // 6: Top quark pole mass
///    mTau;      // 7: Tau pole mass
///
///    // SLHA2
///    mNu3;      // 8: Heaviest neutrino pole mass
///
///    mE;        // 11: Electron pole mass
///    mNu1;      // 12: Lightest neutrino pole mass
///    mMu;       // 13: Muon pole mass
///    mNu2;      // 14: Second lightest neutrino pole mass
///
///    mD;        // 21: d quark running mass in the MSbar scheme at 2 GeV
///    mU;        // 22: u quark running mass in the MSbar scheme at 2 GeV
///    mS;        // 23: s quark running mass in the MSbar scheme at 2 GeV
///    mCmC;      // 24: c quark running mass in the MSbar scheme at mC
///
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (benjamin.farmer@fysik.su.se)
///    \date 2014 Sep - Dec, 2015 Jan - Mar
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2019 Feb
///
///  *********************************************

#include "gambit/Models/models/StandardModel_SLHA2.hpp"

#ifndef __StandardModel_SLHA2_numass_hpp__
#define __StandardModel_SLHA2_numass_hpp__

// Standard Model parameterisation in SLHA2 conventions (1 massive neutrino)
  #define MODEL StandardModel_numass_single
  #define PARENT StandardModel_SLHA2
  START_MODEL

  DEFINEPARS(alphainv, GF, alphaS) // 1,2,3
  DEFINEPARS(mZ)                   // 4

  DEFINEPARS(mE, mMu, mTau)        // 11,13,7
  DEFINEPARS(mNu)                  // 12 (14 and 8 are zero)

  DEFINEPARS(mD, mU)               // 21,22
  DEFINEPARS(mS, mCmC)             // 23,24
  DEFINEPARS(mBmB, mT)             // 5,6

  // CKM parameters
  DEFINEPARS(CKM_lambda, CKM_A, CKM_rhobar, CKM_etabar)

  // PMNS parameters
  DEFINEPARS(theta12, theta23, theta13)
  DEFINEPARS(delta13, alpha1, alpha2)

  // Translation function
  INTERPRET_AS_PARENT_FUNCTION(StandardModel_numass_single_to_StandardModel_SLHA2)
  #undef PARENT
  #undef MODEL

/////////////////

  // Standard Model parameterisation in SLHA2 conventions (3 massive degenerate neutrinos)
  #define MODEL StandardModel_numass_degenerate
  #define PARENT StandardModel_SLHA2
  START_MODEL

  DEFINEPARS(alphainv, GF, alphaS) // 1,2,3
  DEFINEPARS(mZ)                   // 4

  DEFINEPARS(mE, mMu, mTau)        // 11,13,7
  DEFINEPARS(Smu)                  // sum of the mases of 12,14,8

  DEFINEPARS(mD, mU)               // 21,22
  DEFINEPARS(mS, mCmC)             // 23,24
  DEFINEPARS(mBmB, mT)             // 5,6

  // CKM parameters
  DEFINEPARS(CKM_lambda, CKM_A, CKM_rhobar, CKM_etabar)

  // PMNS parameters
  DEFINEPARS(theta12, theta23, theta13)
  DEFINEPARS(delta13, alpha1, alpha2)

  // Translation function
  INTERPRET_AS_PARENT_FUNCTION(StandardModel_numass_degenerate_to_StandardModel_SLHA2)
  #undef PARENT
  #undef MODEL

/////////////////

  // Standard Model parameterisation in SLHA2 conventions (lightest neutrino + splitting)
  #define MODEL StandardModel_numass_split
  #define PARENT StandardModel_SLHA2
  START_MODEL

  DEFINEPARS(alphainv, GF, alphaS) // 1,2,3
  DEFINEPARS(mZ)                   // 4

  DEFINEPARS(mE, mMu, mTau)        // 11,13,7
  DEFINEPARS(mNu_light, dmNu21, dmNu3l)  // Differential neutrino masses

  DEFINEPARS(mD, mU)               // 21,22
  DEFINEPARS(mS, mCmC)             // 23,24
  DEFINEPARS(mBmB, mT)             // 5,6

  // CKM parameters
  DEFINEPARS(CKM_lambda, CKM_A, CKM_rhobar, CKM_etabar)

  // PMNS parameters
  DEFINEPARS(theta12, theta23, theta13)
  DEFINEPARS(delta13, alpha1, alpha2)

  INTERPRET_AS_PARENT_FUNCTION(StandardModel_numass_split_to_StandardModel_SLHA2)
  #undef PARENT
  #undef MODEL

#endif
