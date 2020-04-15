///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
///
///  Model translation functions for Cosmological Energy injection
///
///  *********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Patrick Stöcker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2019 Sep
///  \date 2020 Jan
///
///  *********************************************

#include <string>
#include <vector>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "gambit/Utils/numerical_constants.hpp"

#include "gambit/Models/models/CosmoEnergyInjection.hpp"

#define MODEL AnnihilatingDM_mixture
void MODEL_NAMESPACE::AnnihilatingDM_mixture_to_AnnihilatingDM_general (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(AnnihilatingDM_general)
  logger()<<"Running interpret_as_parent calculations for AnnihilatingDM_mixture --> AnnihilatingDM_general ..."<<LogTags::info<<EOM;

  targetP.setValue("mass", myP.getValue("mass"));
  targetP.setValue("sigmav", myP.getValue("sigmav"));
}

void MODEL_NAMESPACE::energy_injection_spectrum_AnnihilatingDM_mixture(DarkAges::Energy_injection_spectrum& spectrum)
{
  using namespace Pipes::energy_injection_spectrum_AnnihilatingDM_mixture;

  double m = *Param["mass"];
  double BR_el = *Param["BR"];
  double BR_ph = 1.0 - BR_el;

  if (m <= m_electron && BR_el >= std::numeric_limits<double>::epsilon())
  {
    std::ostringstream err;
    err << "The mass of the annihilating dark matter candiate is below the electron mass.";
    err << " No production of e+/e- is possible.";
    model_error().raise(LOCAL_INFO,err.str());
  }

  spectrum.E_el.clear();
  spectrum.E_ph.clear();
  spectrum.spec_el.clear();
  spectrum.spec_ph.clear();

  spectrum.E_el.resize(1,std::max(m-m_electron, std::numeric_limits<double>::min()));
  spectrum.E_ph.resize(1,m);
  spectrum.spec_el.resize(1,BR_el*2e9);
  spectrum.spec_ph.resize(1,BR_ph*2e9);
}
#undef MODEL

#define MODEL AnnihilatingDM_photon
void MODEL_NAMESPACE::AnnihilatingDM_photon_to_AnnihilatingDM_mixture (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(AnnihilatingDM_mixture)
  logger()<<"Running interpret_as_parent calculations for AnnihilatingDM_photon --> AnnihilatingDM_mixture ..."<<LogTags::info<<EOM;

  targetP.setValue("mass", myP.getValue("mass"));
  targetP.setValue("sigmav", myP.getValue("sigmav"));

  // All decays into photons (BR into electrons is 0.0)
  targetP.setValue("BR", 0.0);
}
#undef MODEL

#define MODEL AnnihilatingDM_electron
void MODEL_NAMESPACE::AnnihilatingDM_electron_to_AnnihilatingDM_mixture (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(AnnihilatingDM_mixture)
  logger()<<"Running interpret_as_parent calculations for AnnihilatingDM_electron --> AnnihilatingDM_mixture ..."<<LogTags::info<<EOM;

  targetP.setValue("mass", myP.getValue("mass"));
  targetP.setValue("sigmav", myP.getValue("sigmav"));

  // All decays into electrons
  targetP.setValue("BR", 1.0);
}
#undef MODEL

#define MODEL DecayingDM_mixture
void MODEL_NAMESPACE::DecayingDM_mixture_to_DecayingDM_general (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(DecayingDM_general)
  logger()<<"Running interpret_as_parent calculations for DecayingDM_mixture --> DecayingDM_general ..."<<LogTags::info<<EOM;

  targetP.setValue("mass", myP.getValue("mass"));
  targetP.setValue("lifetime", myP.getValue("lifetime"));
  targetP.setValue("fraction", myP.getValue("fraction"));
}

void MODEL_NAMESPACE::energy_injection_spectrum_DecayingDM_mixture(DarkAges::Energy_injection_spectrum& spectrum)
{
  using namespace Pipes::energy_injection_spectrum_DecayingDM_mixture;

  double m = *Param["mass"];
  double BR_el = *Param["BR"];
  double BR_ph = 1.0 - BR_el;

  if (m <= 2*m_electron && BR_el >= std::numeric_limits<double>::epsilon())
  {
    std::ostringstream err;
    err << "The mass of the decaying dark matter candiate is below twice the electron mass.";
    err << " No production of e+/e- is possible.";
    model_error().raise(LOCAL_INFO,err.str());
  }

  spectrum.E_el.clear();
  spectrum.E_ph.clear();
  spectrum.spec_el.clear();
  spectrum.spec_ph.clear();

  spectrum.E_el.resize(1,std::max(m*0.5-m_electron, std::numeric_limits<double>::min()));
  spectrum.E_ph.resize(1,m*0.5);
  spectrum.spec_el.resize(1,BR_el*2e9);
  spectrum.spec_ph.resize(1,BR_ph*2e9);
}
#undef MODEL

#define MODEL DecayingDM_photon
void MODEL_NAMESPACE::DecayingDM_photon_to_DecayingDM_mixture (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(DecayingDM_mixture)
  logger()<<"Running interpret_as_parent calculations for DecayingDM_photon --> DecayingDM_mixture ..."<<LogTags::info<<EOM;

  targetP.setValue("mass", myP.getValue("mass"));
  targetP.setValue("lifetime", myP.getValue("lifetime"));
  targetP.setValue("fraction", myP.getValue("fraction"));

  // All decays into photons (BR into electrons is 0.0)
  targetP.setValue("BR", 0.0);
}
#undef MODEL

#define MODEL DecayingDM_electron
void MODEL_NAMESPACE::DecayingDM_electron_to_DecayingDM_mixture (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(DecayingDM_mixture)
  logger()<<"Running interpret_as_parent calculations for DecayingDM_electron --> DecayingDM_mixture ..."<<LogTags::info<<EOM;

  targetP.setValue("mass", myP.getValue("mass"));
  targetP.setValue("lifetime", myP.getValue("lifetime"));
  targetP.setValue("fraction", myP.getValue("fraction"));

  // All decays into electrons
  targetP.setValue("BR", 1.0);
}
#undef MODEL
