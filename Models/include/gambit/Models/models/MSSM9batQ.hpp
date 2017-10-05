//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  MSSM9batQ model definition. 
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///   
///  \author Peter Athron
///          (peter.athron@coepp.org.au)
///  \date 2015 Sep
///
///  *********************************************

#ifndef __MSSM9batQ_hpp__
#define __MSSM9batQ_hpp__

// Parent model must be declared first! Include it here to ensure that this happens.
#include "gambit/Models/models/MSSM15atQ.hpp" 

#define MODEL MSSM9batQ
#define PARENT MSSM15atQ
  START_MODEL

  DEFINEPARS(Qin,TanBeta,SignMu,
             mHu2,mHd2,M1,M2,M3)

  DEFINEPARS(mq2_3)
 
  DEFINEPARS(msf2)

  DEFINEPARS(At)

  INTERPRET_AS_PARENT_FUNCTION(MSSM9batQ_to_MSSM15atQ)

#undef PARENT
#undef MODEL

#endif
