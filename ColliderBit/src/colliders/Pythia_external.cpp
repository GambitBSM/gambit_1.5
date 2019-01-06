//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Initialisation function for generic external
///  model collider.
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


#include "gambit/ColliderBit/colliders/ColliderPythia.hpp"

namespace Gambit
{

  namespace ColliderBit
  {

    namespace Pythia_external
    {
      void init(ColliderPythia<Pythia_EM_default::Pythia8::Pythia, Pythia_EM_default::Pythia8::Event>*) {}
    }

  }

}
