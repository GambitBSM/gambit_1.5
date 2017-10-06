//  GAMBIT: Global and Modular BSM Inference Tool
//  *********************************************
//
//  MSSM25 translation function definitions
//  
//  Contains translation functions for
//  MSSM25atQ        --> MSSM30atQ
//  MSSM25atQ_mA     --> MSSM30atQ_mA
//  MSSM25atMGUT     --> MSSM30atMGUT
//  MSSM25atMGUT_mA  --> MSSM30atMGUT_mA
//  MSSM25atMSUSY    --> MSSM30atMSUSY
//  MSSM25atMSUSY_mA --> MSSM30atMSUSY_mA
//
//  *********************************************
//
//  Authors
//  =======
//
//  (add name and date if you modify)
//
//  Ben Farmer
//  2017 Oct
//
//  *********************************************

#include <string>
#include <vector>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"

#include "gambit/Models/models/MSSM25atQ.hpp"
#include "gambit/Models/models/MSSM25atQ_mA.hpp"
#include "gambit/Models/models/MSSM25atMGUT.hpp"
#include "gambit/Models/models/MSSM25atMGUT_mA.hpp"
#include "gambit/Models/models/MSSM25atMSUSY.hpp"
#include "gambit/Models/models/MSSM25atMSUSY_mA.hpp"


using namespace Gambit::Utils;

// General helper translation function definition
namespace Gambit { 
  void MSSM25atX_to_MSSM30atX(const ModelParameters &myP, ModelParameters &targetP)
  {
     // Copy all the common parameters of MSSM25atQ into MSSM30atQ
     targetP.setValues(myP,false);

     // Manually set the parameters which differ
     // slepton trilinear couplings
     // Off-diagonal elements set to zero by parent model
     // First and second generation elements set equal
     targetP.setValue("Ae_1",  myP["Ae_12"] ); // Ae2_11 in MSSM63
     targetP.setValue("Ae_2",  myP["Ae_12"] ); // Ae2_22   " "
     //targetP.setValue("Ae_3",  myP["Ae_3"]  ); // Ae2_33 // Taken care of by common parameter copy

     // down-type trilinear couplings
     // Off-diagonal elements set to zero by parent model
     // First and second generation to zero
     targetP.setValue("Ad_1",  0. );          // Ad2_11 in MSSM63
     targetP.setValue("Ad_2",  0. );          // Ad2_22   " "
     //targetP.setValue("Ad_3",  myP["Ad_3"] ); // Ad2_33 // Taken care of by common parameter copy

     // up-type trilinear couplings
     // Off-diagonal elements set to zero by parent model
     // First and second generation set to zero
     targetP.setValue("Au_1",  0. );          // Au2_11 in MSSM63
     targetP.setValue("Au_2",  0. );          // Au2_22   " "
     // targetP.setValue("Au_3",  myP["Au_3"] ); // Au2_33 // Taken care of by common parameter copy
     
     // Done  
  }
}

/// @{ Interpret-as-parent function definitions
/// These are particularly repetitive so let's define them with the help of a macro
#define DEFINE_IAPFUNC(PARENT) \
void MODEL_NAMESPACE::CAT_3(MODEL,_to_,PARENT) (const ModelParameters &myP, ModelParameters &targetP) \
{ \
   logger()<<"Running interpret_as_parent calculations for " STRINGIFY(MODEL) " --> " STRINGIFY(PARENT) "..."<<LogTags::info<<EOM; \
   MSSM25atX_to_MSSM30atX(myP, targetP); \
} \

#define MODEL MSSM25atQ
DEFINE_IAPFUNC(MSSM30atQ)
#undef MODEL
#define MODEL MSSM25atQ_mA
DEFINE_IAPFUNC(MSSM30atQ_mA)
#undef MODEL
#define MODEL MSSM25atMGUT
DEFINE_IAPFUNC(MSSM30atMGUT)
#undef MODEL
#define MODEL MSSM25atMGUT_mA
DEFINE_IAPFUNC(MSSM30atMGUT_mA)
#undef MODEL
#define MODEL MSSM25atMSUSY
DEFINE_IAPFUNC(MSSM30atMSUSY)
#undef MODEL
#define MODEL MSSM25atMSUSY_mA
DEFINE_IAPFUNC(MSSM30atMSUSY_mA)
#undef MODEL
/// @}

