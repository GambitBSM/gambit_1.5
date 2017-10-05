//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  MSSM20atMSUSY model definition. 
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///   
///  \author Pat Scott  
///          (p.scott@imperial.ac.uk)
///  \date 2015 Sep
///
///  *********************************************

#ifndef __MSSM20atMSUSY_hpp__
#define __MSSM20atMSUSY_hpp__

// Parent model must be declared first! Include it here to ensure that this happens.
#include "gambit/Models/models/MSSM25atMSUSY.hpp"
#include "gambit/Models/models/MSSM20atQ.hpp"

#define MODEL MSSM20atMSUSY
#define PARENT MSSM25atMSUSY
  START_MODEL

  DEFINEPARS(TanBeta,SignMu,
             mHu2,mHd2,M1,M2,M3)

  DEFINEPARS(mq2_12, mq2_3)
 
  DEFINEPARS(ml2_12, ml2_3)

  DEFINEPARS(md2_12, md2_3)

  DEFINEPARS(mu2_12, mu2_3)

  DEFINEPARS(me2_12, me2_3)

  DEFINEPARS(Ae_12, Ae_3)
  
  DEFINEPARS(Ad_3)

  DEFINEPARS(Au_3)

  INTERPRET_AS_PARENT_FUNCTION(MSSM20atMSUSY_to_MSSM25atMSUSY)
  INTERPRET_AS_X_FUNCTION(MSSM20atQ, MSSM20atMSUSY_to_MSSM20atQ)
  INTERPRET_AS_X_DEPENDENCY(MSSM20atQ, unimproved_MSSM_spectrum, Spectrum)

#undef PARENT
#undef MODEL

#endif
