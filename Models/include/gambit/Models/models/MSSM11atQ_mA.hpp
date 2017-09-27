//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  MSSM11atQ_mA model definition.
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

#ifndef __MSSM11atQ_mA_hpp__
#define __MSSM11atQ_mA_hpp__

// Parent model must be declared first! Include it here to ensure that this happens.
#include "gambit/Models/models/MSSM16atQ_mA.hpp"

#define MODEL MSSM11atQ_mA
#define PARENT MSSM16atQ_mA
  START_MODEL

  DEFINEPARS(Qin,TanBeta,
             mA,mu,M1,M2,M3)

  DEFINEPARS(mq2)

  DEFINEPARS(ml2)

  DEFINEPARS(Ae_3)

  DEFINEPARS(Ad_3)

  DEFINEPARS(Au_3)

  INTERPRET_AS_PARENT_FUNCTION(MSSM11atQ_mA_to_MSSM16atQ_mA)

#undef PARENT
#undef MODEL

#endif
