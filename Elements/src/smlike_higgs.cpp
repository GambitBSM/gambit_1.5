//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Helper function to determine which Higgs is
///  most SM-like.
///
///  *********************************************
///
///  Authors:
///  <!-- add name and date if you modify -->
///
///  \author Peter Athron
///          (peter.athron@coepp.org.au)
///  \date 2017
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2017
///
///  *********************************************

#include "gambit/Elements/smlike_higgs.hpp"

namespace Gambit
{

  /// Determine which MSSM higgs is most SM-like.
  /// Needs expansion to work with non-MSSM (e.g. *HDM) models
  int SMlike_higgs_PDG_code(const SubSpectrum& mssm_spec)
  {
    double sa =  - mssm_spec.get(Par::Pole_Mixing,"h0",1,1);
    double ca = mssm_spec.get(Par::Pole_Mixing,"h0",1,2);
    double tb = mssm_spec.get(Par::dimensionless, "tanbeta" );
    double sb = sin(atan(tb));
    double cb = cos(atan(tb));
    //cos (beta - alpha) and sin(beta-alpha)
    double cbma = cb * ca + sb * sa;
    double sbma = sb * ca - cb * ca;
    if(sbma > cbma) return 25;
    return 35;
  }

}
