///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
///
///  MSSM25 translation function definitions
///
///  *********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2017 Aug
///
///  *********************************************

#include <string>
#include <vector>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"

#include "gambit/Models/models/MSSM25atQ_mA.hpp"


// Activate debug output
//#define MSSM25atQ_mA_DBUG

using namespace Gambit::Utils;

#define MODEL MSSM25atQ_mA
  void MODEL_NAMESPACE::MSSM25atQ_mA_to_MSSM30atQ_mA (const ModelParameters &myP, ModelParameters &targetP)
  {
     logger()<<"Running interpret_as_parent calculations for MSSM25atQ_mA --> MSSM30atQ_mA..."<<LogTags::info<<EOM;

     //std::cout<<"Running interpret_as_parent calculations for MSSM25atQ_mA --> MSSM30atQ_mA..."<<std::endl;
     //std::cout<<"myP     = "<<myP.getModelName()<<std::endl;
     //std::cout<<"targetP = "<<targetP.getModelName()<<std::endl;

     // Copy all the common parameters of MSSM25atQ into MSSM30atQ
     targetP.setValues(myP,false); // Set "missing_is_error" flag to false since some MSSM25atQ_mA parameters are not in MSSM30atQ_mA (or rather they are named differently) 
     
     // Manually set the parameters which differ
     // slepton trilinear couplings
     // Off-diagonal elements set to zero by parent model
     // First and second generation elements set equal, third generation left free
     targetP.setValue("Ae_1",  myP["Ae_12"] ); // Ae2_11 in MSSM63
     targetP.setValue("Ae_2",  myP["Ae_12"] ); // Ae2_22   " "

     // down-type trilinear couplings
     // Off-diagonal elements set to zero by parent model
     // First and second generation to zero, third generation left free
     targetP.setValue("Ad_1",  0. );          // Ad2_11 in MSSM63
     targetP.setValue("Ad_2",  0. );          // Ad2_22   " "

     // up-type trilinear couplings
     // Off-diagonal elements set to zero by parent model
     // First and second generation set to zero, third generation left free
     targetP.setValue("Au_1",  0. );          // Au2_11 in MSSM63
     targetP.setValue("Au_2",  0. );          // Au2_22   " "

     #ifdef MSSM25atQ_mA_DBUG
       std::cout << "MSSM25atQ_mA parameters:" << myP << std::endl;
       std::cout << "MSSM30atQ_mA parameters:" << targetP << std::endl;
     #endif
  }

#undef MODEL
