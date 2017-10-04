//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  MSSM20atMGUT_mA model definition.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2017 Aug
///
///  *********************************************

#ifndef __MSSM20atMGUT_mA_hpp__
#define __MSSM20atMGUT_mA_hpp__

// Parent model must be declared first! Include it here to ensure that this happens.
#include "gambit/Models/models/MSSM25atMGUT_mA.hpp"
#include "gambit/Models/models/MSSM20atMGUT.hpp"

#define MODEL MSSM20atMGUT_mA
#define PARENT MSSM25atMGUT_mA
  START_MODEL

  DEFINEPARS(TanBeta,
             mA,mu,M1,M2,M3)

  DEFINEPARS(mq2_12, mq2_3)

  DEFINEPARS(ml2_12, ml2_3)

  DEFINEPARS(md2_12, md2_3)

  DEFINEPARS(mu2_12, mu2_3)

  DEFINEPARS(me2_12, me2_3)

  DEFINEPARS(Ae_12, Ae_3)

  DEFINEPARS(Ad_3)

  DEFINEPARS(Au_3)

  INTERPRET_AS_PARENT_FUNCTION(MSSM20atMGUT_mA_to_MSSM25atMGUT_mA)
  INTERPRET_AS_X_FUNCTION(MSSM20atMGUT, MSSM20atMGUT_mA_to_MSSM20atMGUT)
  INTERPRET_AS_X_DEPENDENCY(MSSM20atMGUT, unimproved_MSSM_spectrum, Spectrum)

#undef PARENT
#undef MODEL

#endif
