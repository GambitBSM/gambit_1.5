//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Definitions of container classes
///  for the DarkAges backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stpecker@physik.rwth-aachen.de)
///  \date 2018 Mar
///
///  *********************************************

#ifndef __DarkAges_types_hpp__
#define __DarkAges_types_hpp__

namespace Gambit
{
  namespace DarkAges
  {
    struct injectionSpectrum
    {
      std::vector<double> E;
      std::vector<double> spec_el;
      std::vector<double> spec_ph;
    };

    struct fz_table
    {
      std::vector<double> redshift;
      std::vector<double> f_heat;
      std::vector<double> f_lya;
      std::vector<double> f_hion;
      std::vector<double> f_heion;
      std::vector<double> f_lowe;
    };
  }
}

#endif // defined __DarkAges_types_hpp__
