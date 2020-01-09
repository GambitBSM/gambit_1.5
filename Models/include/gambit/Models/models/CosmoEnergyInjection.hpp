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
///  \date 2020 Jan
///
///  *********************************************

#ifndef __CosmoEnergyInjection_hpp__
#define __CosmoEnergyInjection_hpp__

// General model of annihilating dark matter with monochromatic injection of e+/- and photons
#define MODEL AnnihilatingDM_general
  START_MODEL
  DEFINEPARS(mass)      // mass of dcdm candidate [GeV]
  DEFINEPARS(sigmav)    // thermally averaged cross section [cm^3 s^-1]
  DEFINEPARS(BR)        // Branching ratio into electrons (1-BR into photons)

  MAP_TO_CAPABILITY(mass, DM_mass)
  MAP_TO_CAPABILITY(sigmav, sigmav)
#undef MODEL

// Special case of 'AnnihilatingDM_general': Only annihilate into photons
#define MODEL AnnihilatingDM_photon
#define PARENT AnnihilatingDM_general
  START_MODEL
  DEFINEPARS(mass)      // mass of dcdm candidate [GeV]
  DEFINEPARS(sigmav)    // thermally averaged cross section [cm^3 s^-1]

  INTERPRET_AS_PARENT_FUNCTION(AnnihilatingDM_photon_to_AnnihilatingDM_general)
#undef PARENT
#undef MODEL

// Special case of 'AnnihilatingDM_general': Only annihilate into electrons
#define MODEL AnnihilatingDM_electron
#define PARENT AnnihilatingDM_general
  START_MODEL
  DEFINEPARS(mass)      // mass of dcdm candidate [GeV]
  DEFINEPARS(sigmav)    // thermally averaged cross section [cm^3 s^-1]

  INTERPRET_AS_PARENT_FUNCTION(AnnihilatingDM_electron_to_AnnihilatingDM_general)
#undef PARENT
#undef MODEL

// General model of decaying dark matter with monochromatic injection of e+/- and photons
#define MODEL DecayingDM_general
  START_MODEL
  DEFINEPARS(mass)      // mass of dcdm candidate [GeV]
  DEFINEPARS(lifetime)  // lifetime of dcdm candiate [s]
  DEFINEPARS(fraction)  // rho_dcdm / rho_cdm in the infinite past
  DEFINEPARS(BR)        // Branching ratio into electrons (1-BR into photons)

  MAP_TO_CAPABILITY(mass, DM_mass)
  MAP_TO_CAPABILITY(lifetime, lifetime)
  MAP_TO_CAPABILITY(fraction, DM_fraction)
#undef MODEL

// Special case of 'DecayingDM_general': Only decay into photons
#define MODEL DecayingDM_photon
#define PARENT DecayingDM_general
  START_MODEL
  DEFINEPARS(mass)      // mass of dcdm candidate [GeV]
  DEFINEPARS(lifetime)  // lifetime of dcdm candiate [s]
  DEFINEPARS(fraction)  // rho_dcdm / rho_cdm in the infinite past

  INTERPRET_AS_PARENT_FUNCTION(DecayingDM_photon_to_DecayingDM_general)
#undef PARENT
#undef MODEL

// Special case of 'DecayingDM_general': Only decay into electrons
#define MODEL DecayingDM_electron
#define PARENT DecayingDM_general
  START_MODEL
  DEFINEPARS(mass)      // mass of dcdm candidate [GeV]
  DEFINEPARS(lifetime)  // lifetime of dcdm candiate [s]
  DEFINEPARS(fraction)  // rho_dcdm / rho_cdm in the infinite past

  INTERPRET_AS_PARENT_FUNCTION(DecayingDM_electron_to_DecayingDM_general)
#undef PARENT
#undef MODEL

#endif
