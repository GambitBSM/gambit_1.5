//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  A ColliderPythia init function for a basic
///  SUSY @ 13TeV LHC scenario. Note that additional
///  Pythia settings may still be applied externally
///  via yaml file input.
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

    namespace Pythia_SUSY_LHC_13TeV
    {

      void init(ColliderPythia<Pythia_default::Pythia8::Pythia, Pythia_default::Pythia8::Event>* specializeMe)
      {
        specializeMe->addToSettings("Beams:eCM = 13000");
        specializeMe->addToSettings("Main:timesAllowErrors = 1000");
        specializeMe->addToSettings("SUSY:all = on");
        specializeMe->addToSettings("Random:setSeed = on");
      }

    }

  }

}
