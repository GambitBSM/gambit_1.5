///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
///
///  DiracDM model source file. 
///  The Dirac fermion mass mV is the tree-level mass, 
///  lambda_hV is the tree-level DM-Higgs coupling and
///  cosXI = 0 (1) corresponds to scalar (pseudoscalar) 
///  DM-Higgs interaction.
///  (file format is based on SingletDM.cpp)
///
///  *********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Ankit Beniwal
///          (ankit.beniwal@adelaide.edu.com)
///  \date 2016 Aug
///
///  *********************************************

#include <string>
#include <vector>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"

#include "gambit/Models/models/DiracDM.hpp"
#include "gambit/Elements/spectrum.hpp"

using namespace Gambit::Utils;

// Need to define MODEL and PARENT in order for helper macros to work correctly
#define MODEL DiracDM

  // No translation function exist for the model

#undef MODEL
