//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Source code for types for module CosmoBit.
///  For instructions on adding new types, see
///  the corresponding header.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2017 Nov
///
///  \author Selim Hotinli
///  \date 2018 Jan
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 Oct
///
///  \author Sebastian Hoof
///          (hoof@uni-goettingen.de)
///  \date 2020 Mar
///
///  \author Pat Scott
///          (pat.scott@uq.edu.au)
///  \date 2020 Apr
///
///  *********************************************

#include <string>
#include <iostream>
#include <valarray>

#include "gambit/CosmoBit/CosmoBit_types.hpp"
#include "gambit/CosmoBit/CosmoBit_utils.hpp"
#include "gambit/Utils/numerical_constants.hpp"

namespace Gambit
{
  namespace CosmoBit
  {

    // Default constructor for multimode inputs
    // N.B. The constructor is set up to support single field inflation models with instant reheating but allow the inclusion
    //      of more complex models. Some default MultiModeCode (MMC) inputs below might have to be adjusted to facilitate this.
    Multimode_inputs::Multimode_inputs()
    {
      num_inflaton = 1; // Assume single field inflation
      instreheat = 1; // Use the instant reheating approximation
      // Using the instant reheating approximation makes N_pivot a derived parameter and the input N_pivot has no effect
      // We fix it to the dummy value below
      N_pivot = 50;
      // CAVE Changing the default value of 'slowroll_infl_end' requires defining a custom condition for the end of inflation in MMC!
      slowroll_infl_end = 1; // = true, i.e. stop inflation when slow roll parameters = 1
      // Control the output of analytic approximations for comparison. We do not use these.
      use_deltaN_SR = 0; // = false, i.e. MMC will not calculate deltaN observables (assumes slow roll & sum-separable potentials) at the pivot scale
      use_horiz_cross_approx = 0; // = false, i.e. do not ignore the horizon-crossing-approximation for the above
      evaluate_modes = 1; // = true, i.e. evalute modes and do not just rely on background evolution
      get_runningofrunning = 0; // = false, i.e. do not compute the dervative of the spectral index w.r.t. ln(k)
      // Set the initial conditions for the inflation field(s).
      // N.B. For single field inflation, MMC determines the parameters below self-consistenly; choose sensible entries as starting point
      phi_init0 = {10.0};
      dphi_init0 = {1.0};
    };

    // return Parametrised_ps members A_s, n_s, r, and N_pivot as str to double map for printing
    map_str_dbl Parametrised_ps::get_parametrised_ps_map()
    {
      map_str_dbl result;
      result["ln10A_s"] = ln10A_s;
      result["n_s"] = n_s;
      result["r"] = r;
      result["N_pivot"] = N_pivot;

      return result;

    };

    void Primordial_ps::fill_k(double *k_array, int len)
    {
      std::vector<double> K(k_array, k_array+len);
      k = std::move(K);
      vec_size = len;
      //let's leave the debug print a while longer though, you never know ..
      //for( int i =0; i<len;i++) {std::cout << k[i] << std::endl;};
    }

    void Primordial_ps::fill_P_s(double *P_s_array, int len)
    {
      std::vector<double> ps(P_s_array, P_s_array+len);
      P_s = std::move(ps);
    }
    void Primordial_ps::fill_P_s_iso(double *P_s_iso_array, int len)
    {
      std::vector<double> psi(P_s_iso_array, P_s_iso_array+len);
      P_s_iso = std::move(psi);
    }
    void Primordial_ps::fill_P_t(double *P_t_array, int len)
    {
      std::vector<double> pt(P_t_array, P_t_array+len);
      P_t = std::move(pt);
    }

  }
}
