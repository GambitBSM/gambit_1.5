//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  MSSM9atQ_mA model definition.
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

#ifndef __MSSM9atQ_mA_hpp__
#define __MSSM9atQ_mA_hpp__

// Parent and friend models must be declared first! Include them here to ensure that this happens.
#include "gambit/Models/models/MSSM10atQ_mA.hpp"
#include "gambit/Models/models/MSSM10batQ_mA.hpp"

#define MODEL MSSM9atQ_mA
#define PARENT MSSM10atQ_mA
  START_MODEL

  DEFINEPARS(Qin,TanBeta,
             mA,mu,M1,M2,M3)

  DEFINEPARS(mf2)

  DEFINEPARS(Ad_3)

  DEFINEPARS(Au_3)

  INTERPRET_AS_PARENT_FUNCTION(MSSM9atQ_mA_to_MSSM10atQ_mA)
  INTERPRET_AS_X_FUNCTION(MSSM10batQ_mA, MSSM9atQ_mA_to_MSSM10batQ_mA)

#undef PARENT
#undef MODEL

#endif
