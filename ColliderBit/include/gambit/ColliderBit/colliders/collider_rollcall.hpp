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

namespace Gambit
{

  namespace ColliderBit
  {

    template<typename ColliderPythiaT>
    void Pythia_external_init(ColliderPythiaT*) {}

    template<typename ColliderPythiaT>
    void Pythia_SUSY_LHC_8TeV_init(ColliderPythiaT* specializeMe)
    {
      specializeMe->addToSettings("Beams:eCM = 8000");
      specializeMe->addToSettings("Main:timesAllowErrors = 1000");
      specializeMe->addToSettings("SUSY:all = on");
      specializeMe->addToSettings("Random:setSeed = on");
    }

    template<typename ColliderPythiaT>
    void Pythia_glusq_LHC_8TeV_init(ColliderPythiaT* specializeMe)
    {
      Pythia_SUSY_LHC_8TeV_init(specializeMe);
      specializeMe->addToSettings("SUSY:idA = 1000021");
      specializeMe->addToSettings("SUSY:idVecB = 1000001, 1000002, 1000003, 1000004, 2000001, 2000002, 2000003, 2000004");
    }

    template<typename ColliderPythiaT>
    void Pythia_SUSY_LHC_13TeV_init(ColliderPythiaT* specializeMe)
    {
      specializeMe->addToSettings("Beams:eCM = 13000");
      specializeMe->addToSettings("Main:timesAllowErrors = 1000");
      specializeMe->addToSettings("SUSY:all = on");
      specializeMe->addToSettings("Random:setSeed = on");
    }


  }
}
