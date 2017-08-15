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

#define MODEL StandardModel_Neutrinos

  START_MODEL

  // Mass parameters
  DEFINEPARS(min_mass, ordering, md21, md31, md23)

  // PMNS parameters
  // TODO: Change this names when the conflict with StandardModel_SLHA2 is fixed
  DEFINEPARS(t12, t23, t13)
  DEFINEPARS(delta, a1, a2)

#undef MODEL

#endif /* __StandardModel_Neutrinos_hpp__ */
