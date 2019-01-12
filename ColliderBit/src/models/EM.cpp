//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  External Model-specific sources for ColliderBit.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2019 Jan
///
///  *********************************************

#include "gambit/ColliderBit/getColliderPythia.hpp"
#include "gambit/ColliderBit/generateEventColliderPythia.hpp"
#include "gambit/ColliderBit/getBuckFast.hpp"
#include "gambit/ColliderBit/smearEvent.hpp"


namespace Gambit
{
  namespace ColliderBit
  {

    // Get Monte Carlo event generator
    GET_SPECIFIC_PYTHIA(getPythia_EM, Pythia_EM_default, EM)
    GET_SPECIFIC_PYTHIA_FROM_SLHA(getPythia_EMFileReader, Pythia_EM_default, EM)
    GET_PYTHIA_AS_BASE_COLLIDER(getPythia_EMAsBase)

    // Run event generator
    GET_PYTHIA_EVENT(generateEventPythia_EM, Pythia_EM_default::Pythia8::Event)

    // Get detector simulations
    GET_BUCKFAST_AS_BASE_DETECTOR(getBuckFastATLASPythia_EM, Pythia_EM_default::Pythia8::Event, ATLAS, )
    GET_BUCKFAST_AS_BASE_DETECTOR(getBuckFastATLASmultieffPythia_EM, Pythia_EM_default::Pythia8::Event, ATLAS, multieff)
    GET_BUCKFAST_AS_BASE_DETECTOR(getBuckFastCMSPythia_EM, Pythia_EM_default::Pythia8::Event, CMS, )
    GET_BUCKFAST_AS_BASE_DETECTOR(getBuckFastCMSmultieffPythia_EM, Pythia_EM_default::Pythia8::Event, CMS, multieff)
    GET_BUCKFAST_AS_BASE_DETECTOR(getBuckFastIdentityPythia_EM, Pythia_EM_default::Pythia8::Event, Identity, )

    // Run detector simulations
    SMEAR_EVENT(smearEventATLAS_EM, ATLAS)
    SMEAR_EVENT(smearEventATLASmultieff_EM, ATLASmultieff)
    SMEAR_EVENT(smearEventCMS_EM, CMS)
    SMEAR_EVENT(smearEventCMSmultieff_EM, CMSmultieff)
    SMEAR_EVENT(copyEvent_EM, Identity)

  }
}
