// GAMBIT: Global and Modular BSM Inference Tool
//
// *********************************************
//
// Sterile RH Neutrino Model with differential masses
//
// *********************************************
//
// Authors
// =======
//
// (add name and date if you modify)
//
// \author Tomas Gonzalo
//         (t.e.gonzalo@fys.uio.no)
// \date 2017 December
//
// *********************************************

#ifndef __SN_differential_hpp__
#define __SN_differential_hpp__

#include "gambit/Models/models/SN_dev.hpp"

#define MODEL SN_differential
#define PARENT SN_dev
  START_MODEL
  DEFINEPARS(M_1, delta_M_2, delta_M_3, ReOm23, ImOm23, ReOm13, ImOm13, ReOm12, ImOm12, L_Ge, L_Xe)
  INTERPRET_AS_PARENT_FUNCTION(SN_differential_to_SN_dev)
#undef PARENT
#undef MODEL

#endif
