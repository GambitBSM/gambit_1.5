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

#define DECLARE_COLLIDER(NS, TYPE) namespace NS { void init(TYPE*); }


namespace Gambit
{

  namespace ColliderBit
  {

    template<typename T1, typename T2>
    class ColliderPythia;

    /// Typedefs for each Pythia collider
    /// @{
    typedef ColliderPythia<Pythia_default::Pythia8::Pythia, Pythia_default::Pythia8::Event>       ColliderPythia_defaultversion;
    typedef ColliderPythia<Pythia_EM_default::Pythia8::Pythia, Pythia_EM_default::Pythia8::Event> ColliderPythia_EM_defaultversion;
    /// @{

    DECLARE_COLLIDER(Pythia_external,       ColliderPythia_EM_defaultversion)
    DECLARE_COLLIDER(Pythia_SUSY_LHC_8TeV,  ColliderPythia_defaultversion)
    DECLARE_COLLIDER(Pythia_glusq_LHC_8TeV, ColliderPythia_defaultversion)
    DECLARE_COLLIDER(Pythia_SUSY_LHC_13TeV, ColliderPythia_defaultversion)

  }
}
