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
///          (stoecker@physik.rwth-aachen.de)
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
      // Flag whether DarkAges calculates f_c(z) per injection channel or an
      // effective f(z) which is multiplied with repartition fucntions
      // chi_c(x_e,z) later on (in classy)
      bool f_eff_mode = false;

      // Redshift vector (filled in both modes)
      std::vector<double> redshift;

      // f_c(z) seperated by injection channels
      // (filled only if f_eff_mode is false)
      std::vector<double> f_heat;
      std::vector<double> f_lya;
      std::vector<double> f_hion;
      std::vector<double> f_heion;
      std::vector<double> f_lowe;

      // f_eff(z)
      // (filled only if f_eff_mode is true)
      std::vector<double> f_eff;
    };
  }
}

#endif // defined __DarkAges_types_hpp__
