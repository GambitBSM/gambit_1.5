//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  MSSM20atQ translation function definitions. 
///
///  Contains the interpret-as-parent translation 
///  functions for:
///
///  MSSM20atQ        --> MSSM25atQ
///  MSSM20atQ_mA     --> MSSM25atQ_mA
///  MSSM20atMGUT     --> MSSM25atMGUT
///  MSSM20atMGUT_mA  --> MSSM25atMGUT_mA
///  MSSM20atMSUSY    --> MSSM25atMSUSY
///  MSSM20atMSUSY_mA --> MSSM25atMSUSY_mA
///
///  As well as the interpret-as-friend translation
///  functions for:
///
///  MSSM20atQ_mA     --> MSSM20atQ
///  MSSM20atMGUT_mA  --> MSSM20atMGUT
///  MSSM20atMSUSY_mA --> MSSM20atMSUSY
///
///  and
///
///  MSSM20atMGUT  --> MSSM20atQ
///  MSSM20atMSUSY --> MSSM20atQ
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///   
///  \author Pat Scott  
///          (p.scott@imperial.ac.uk)
///  \date 2015 Sep
///
///  \author Ben Farmer
///          (benjamin.farmer@fysik.su.se)
///  \date 2017 Oct
///
///  *********************************************

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"

#include "gambit/Models/models/MSSM63atQ.hpp" // Contains declaration of MSSM_mA_to_MSSM_mhud and MSSMatX_to_MSSMatQ functions
#include "gambit/Models/models/MSSM20atQ.hpp"
#include "gambit/Models/models/MSSM20atQ_mA.hpp"
#include "gambit/Models/models/MSSM20atMGUT.hpp"
#include "gambit/Models/models/MSSM20atMGUT_mA.hpp"
#include "gambit/Models/models/MSSM20atMSUSY.hpp"
#include "gambit/Models/models/MSSM20atMSUSY_mA.hpp"

#include "gambit/Elements/spectrum.hpp"


// General helper translation function definition
namespace Gambit { 
  void MSSM20atX_to_MSSM25atX(const ModelParameters &myP, ModelParameters &targetP)
  {
     // Send all parameter values upstream to matching parameters in parent.
     // Ignore that some parameters don't exist in the parent, these are set below.
     targetP.setValues(myP,false);

     // RH squark soft masses, gen 1 and 2
     targetP.setValue("mq2_1",  myP["mq2_12"] ); // mq2_11 in MSSM63
     targetP.setValue("mq2_2",  myP["mq2_12"] ); // mq2_22   " "
     // RH slepton soft masses, gen 1 and 2
     targetP.setValue("ml2_1",  myP["ml2_12"] ); // ml2_11 in MSSM63
     targetP.setValue("ml2_2",  myP["ml2_12"] ); // ml2_22   " "
     // LH down-type squark soft masses
     targetP.setValue("md2_1",  myP["md2_12"] ); // ml2_11 in MSSM63
     targetP.setValue("md2_2",  myP["md2_12"] ); // ml2_22   " "
     // LH up-type squark soft masses
     targetP.setValue("mu2_1",  myP["mu2_12"] ); // mu2_11 in MSSM63
     targetP.setValue("mu2_2",  myP["mu2_12"] ); // mu2_22   " "
     // LH charged slepton soft masses
     targetP.setValue("me2_1",  myP["me2_12"] ); // me2_11 in MSSM63
     targetP.setValue("me2_2",  myP["me2_12"] ); // me2_22   " "
     // Done
  }
}
#undef MODEL

/// @{ Interpret-as-parent function definitions
/// These are particularly repetitive so let's define them with the help of a macro
#define DEFINE_IAPFUNC(PARENT) \
void MODEL_NAMESPACE::CAT_3(MODEL,_to_,PARENT) (const ModelParameters &myP, ModelParameters &targetP) \
{ \
   logger()<<"Running interpret_as_parent calculations for " STRINGIFY(MODEL) " --> " STRINGIFY(PARENT) "..."<<LogTags::info<<EOM; \
   MSSM20atX_to_MSSM25atX(myP, targetP); \
} \

#define MODEL MSSM20atQ
DEFINE_IAPFUNC(MSSM25atQ)
#undef MODEL
#define MODEL MSSM20atQ_mA
DEFINE_IAPFUNC(MSSM25atQ_mA)
#undef MODEL
#define MODEL MSSM20atMGUT
DEFINE_IAPFUNC(MSSM25atMGUT)
#undef MODEL
#define MODEL MSSM20atMGUT_mA
DEFINE_IAPFUNC(MSSM25atMGUT_mA)
#undef MODEL
#define MODEL MSSM20atMSUSY
DEFINE_IAPFUNC(MSSM25atMSUSY)
#undef MODEL
#define MODEL MSSM20atMSUSY_mA
DEFINE_IAPFUNC(MSSM25atMSUSY_mA)
#undef MODEL
/// @}

/// @{ Interpret-as-friend (mA parameterisations to primary parameterisations)
#define MODEL MSSM20atQ_mA
void MODEL_NAMESPACE::MSSM20atQ_mA_to_MSSM20atQ(const ModelParameters &myP, ModelParameters &targetP)
{
   logger()<<"Running interpret_as_X calculations for MSSM20atQ_mA --> MSSM20atQ."<<LogTags::info<<EOM;
   USE_MODEL_PIPE(MSSM20atQ) // Need the pipe for the TARGET model
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSM_mA_to_MSSM_mhud(myP, targetP, HE);
}
#undef MODEL

#define MODEL MSSM20atMGUT_mA
void MODEL_NAMESPACE::MSSM20atMGUT_mA_to_MSSM20atMGUT(const ModelParameters &myP, ModelParameters &targetP)
{
   logger()<<"Running interpret_as_X calculations for MSSM20atMGUT_mA --> MSSM20atMGUT."<<LogTags::info<<EOM;
   USE_MODEL_PIPE(MSSM20atMGUT) // Need the pipe for the TARGET model
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSM_mA_to_MSSM_mhud(myP, targetP, HE);
}
#undef MODEL

#define MODEL MSSM20atMSUSY_mA
void MODEL_NAMESPACE::MSSM20atMSUSY_mA_to_MSSM20atMSUSY(const ModelParameters &myP, ModelParameters &targetP)
{
   logger()<<"Running interpret_as_X calculations for MSSM20atMSUSY_mA --> MSSM20atMSUSY."<<LogTags::info<<EOM;
   USE_MODEL_PIPE(MSSM20atMSUSY) // Need the pipe for the TARGET model
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSM_mA_to_MSSM_mhud(myP, targetP, HE);
}
#undef MODEL
/// @}

/// @{ Interpret-as-friend (MGUT and MSUSY to Q, in primary parameterisation only)
#define MODEL MSSM20atMGUT
void MODEL_NAMESPACE::MSSM20atMGUT_to_MSSM20atQ (const ModelParameters &myP, ModelParameters &targetP)
{
   USE_MODEL_PIPE(MSSM20atQ) // Need pipe for TARGET model
   logger()<<"Running interpret_as_X calculations for MSSM20atMGUT --> MSSM20atQ..."<<LogTags::info<<EOM;
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSMatX_to_MSSMatQ(myP, targetP, HE);
}
#undef MODEL

#define MODEL MSSM20atMSUSY
void MODEL_NAMESPACE::MSSM20atMSUSY_to_MSSM20atQ (const ModelParameters &myP, ModelParameters &targetP)
{
   USE_MODEL_PIPE(MSSM20atQ) // Need pipe for TARGET model
   logger()<<"Running interpret_as_X calculations for MSSM20atMSUSY --> MSSM20atQ..."<<LogTags::info<<EOM;
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSMatX_to_MSSMatQ(myP, targetP, HE);
}
#undef MODEL
