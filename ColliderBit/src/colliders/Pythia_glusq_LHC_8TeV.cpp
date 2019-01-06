//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  A ColliderPythia init function for
///  gluino-squark production @ 8TeV LHC scenario.
///  This "inherits" Pythia_SUSY_LHC_8TeV by
///  explicitly calling its init before changing
///  additional settings. Note that additional
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

    namespace Pythia_glusq_LHC_8TeV
    {

      void init(ColliderPythia<Pythia_default::Pythia8::Pythia, Pythia_default::Pythia8::Event>* specializeMe)
      {
        Pythia_SUSY_LHC_8TeV::init(specializeMe);
        specializeMe->addToSettings("SUSY:idA = 1000021");
        specializeMe->addToSettings("SUSY:idVecB = 1000001, 1000002, 1000003, 1000004, 2000001, 2000002, 2000003, 2000004");
      }

    }

  }

}
