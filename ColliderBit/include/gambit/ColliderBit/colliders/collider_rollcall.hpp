//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  A list of all colliders recognised by ColliderBit,
///  and specialisations of event generators for each one.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///  (p.scott@imperial.ac.uk)
///  \date Jan 2019
///
///  *********************************************


#pragma once

// Forward declare the SpecializablePythia template; we don't care about its implementations at this stage.
template<T1, T2>
class SpecializablePythia;

#define DECLARE_COLLIDER(NS, TYPE) namespace NS { void init(TYPE*); }

namespace Gambit
{

  namespace ColliderBit
  {

    DECLARE_COLLIDER(Pythia_external,       SpecializablePythia<Pythia_EM_default::Pythia8::Pythia, Pythia_EM_default::Pythia8::Event>)
    DECLARE_COLLIDER(Pythia_SUSY_LHC_8TeV,  SpecializablePythia<Pythia_default::Pythia8::Pythia, Pythia_default::Pythia8::Event>)
    DECLARE_COLLIDER(Pythia_glusq_LHC_8TeV, SpecializablePythia<Pythia_default::Pythia8::Pythia, Pythia_default::Pythia8::Event>)
    DECLARE_COLLIDER(Pythia_SUSY_LHC_13TeV, SpecializablePythia<Pythia_default::Pythia8::Pythia, Pythia_default::Pythia8::Event>)

  }
}
