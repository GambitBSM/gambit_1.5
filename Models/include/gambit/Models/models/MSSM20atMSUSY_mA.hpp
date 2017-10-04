//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  MSSM20atMSUSY_mA model definition.
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

#ifndef __MSSM20atMSUSY_mA_hpp__
#define __MSSM20atMSUSY_mA_hpp__

// Parent model must be declared first! Include it here to ensure that this happens.
#include "gambit/Models/models/MSSM25atMSUSY_mA.hpp"
#include "gambit/Models/models/MSSM20atMSUSY.hpp"

#define MODEL MSSM20atMSUSY_mA
#define PARENT MSSM25atMSUSY_mA
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

  INTERPRET_AS_PARENT_FUNCTION(MSSM20atMSUSY_mA_to_MSSM25atMSUSY_mA)
  INTERPRET_AS_X_FUNCTION(MSSM20atMSUSY, MSSM20atMSUSY_mA_to_MSSM20atMSUSY)
  INTERPRET_AS_X_DEPENDENCY(MSSM20atMSUSY, unimproved_MSSM_spectrum, Spectrum)

#undef PARENT
#undef MODEL

#endif
