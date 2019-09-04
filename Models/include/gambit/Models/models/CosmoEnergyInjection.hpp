//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Model definitions for Cosmological Energy injection
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2019 Sep
///
///  *********************************************

#ifndef __CosmoEnergyInjection_hpp__
#define __CosmoEnergyInjection_hpp__

// A dummy model for decaying dark matter
#define MODEL TestDecayingDM
  START_MODEL
  DEFINEPARS(mass,lifetime,fraction,BR)
  MAP_TO_CAPABILITY(mass, DM_mass)
  MAP_TO_CAPABILITY(lifetime, lifetime)
  MAP_TO_CAPABILITY(fraction, DM_fraction)
#undef MODEL

#endif
