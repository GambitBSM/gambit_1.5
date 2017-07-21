//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Standard Model Neutrino parameters.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2017 July
///
///  *********************************************


#ifndef __StandardModel_Neutrinos_hpp__
#define __StandardModel_Neutrinos_hpp__

#include "gambit/Models/models/StandardModel_SLHA2.hpp"

#define MODEL StandardModel_Neutrinos

  START_MODEL

  // Mass parameters
  DEFINEPARS(min_mass, ordering, md21, md31, md23)

  // PMNS parameters
  DEFINEPARS(theta12, theta23, theta13)
  DEFINEPARS(delta, alpha1, alpha2)

#undef MODEL

#endif /* __StandardModel_Neutrinos_hpp__ */
