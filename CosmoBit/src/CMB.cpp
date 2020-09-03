//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  CosmoBit routines relating to the CMB.
///
///  Routines include extracting CMB spectra and 
///  computing effects of energy injections.
///
///  Does *not* contain the Planck likelihoods; 
///  these live in CosmoBit/src/Planck.cpp.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Selim C. Hotinli
///          (selim.hotinli14@pimperial.ac.uk)
///  \date 2017 Jul
///  \date 2018 May
///  \date 2018 Aug - Sep
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2017 Nov
///  \date 2018 Jan - May
///  \date 2019 Jan, Feb, June, Nov
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 June
///  \date 2019 Mar,June
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date 2019 June, Nov
///
///  \author Sebastian Hoof
///          (hoof@uni-goettingen.de)
///  \date 2020 Mar
///
///  \author Pat Scott
///          (pat.scott@uq.edu.au)
///  \date 2018 Mar
///  \date 2019 Jul
///  \date 2020 Apr
///
///  *********************************************

#include <gsl/gsl_spline.h>

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/CosmoBit/CosmoBit_rollcall.hpp"
#include "gambit/CosmoBit/CosmoBit_types.hpp"
#include "gambit/CosmoBit/CosmoBit_utils.hpp"

namespace Gambit
{

  namespace CosmoBit
  {
    using namespace LogTags;

    /**********/
    /* Classy */
    /**********/

    /// Getter functions for CL spectra from classy.

    /* UNLENSED SPECTRA */

    /// Temperature autocorrelation 
    void class_get_unlensed_Cl_TT(std::vector<double>& result)
    {
      using namespace Pipes::class_get_unlensed_Cl_TT;
      result = BEreq::class_get_unlensed_cl("tt");

      // Loop through the TT spectrum and check if there is a negative value. If so, invalidate.
      for (auto it=result.begin(); it != result.end(); it++)
      {
        if (*it < 0.0)
        {
          invalid_point().raise("Found a negative value in the TT spectrum.");
        }
      }
    }

    /// Temperature & E-mode cross-correlation
    void class_get_unlensed_Cl_TE(std::vector<double>& result)
    {
      using namespace Pipes::class_get_unlensed_Cl_TE;
      result = BEreq::class_get_unlensed_cl("te");
    }

    /// E-mode autocorrelation 
    void class_get_unlensed_Cl_EE(std::vector<double>& result)
    {
      using namespace Pipes::class_get_unlensed_Cl_EE;
      result = BEreq::class_get_unlensed_cl("ee");

      // Loop through the EE spectrum and check if there is a negative value. If so, invalidate.
      for (auto it=result.begin(); it != result.end(); it++)
      {
        if (*it < 0.0)
        {
          invalid_point().raise("Found a negative value in the EE spectrum.");
        }
      }
    }

    /// B-mode autocorrelation 
    void class_get_unlensed_Cl_BB(std::vector<double>& result)
    {
      using namespace Pipes::class_get_unlensed_Cl_BB;
      result = BEreq::class_get_unlensed_cl("bb");

      // Loop through the BB spectrum and check if there is a negative value. If so, invalidate.
      for (auto it=result.begin(); it != result.end(); it++)
      {
        if (*it < 0.0)
        {
          invalid_point().raise("Found a negative value in the BB spectrum.");
        }
      }
    }

    /// Lensing (Phi) autocorrelation
    void class_get_unlensed_Cl_PhiPhi(std::vector<double>& result)
    {
      using namespace Pipes::class_get_unlensed_Cl_PhiPhi;
      result = BEreq::class_get_unlensed_cl("pp");

      // Loop through the PhiPhi spectrum and check if there is a negative value. If so, invalidate.
      for (auto it=result.begin(); it != result.end(); it++)
      {
        if (*it < 0.0)
        {
          invalid_point().raise("Found a negative value in the PhiPhi spectrum.");
        }
      }
    }

    /* LENSED SPECTRA */

    /// Temperature autocorrelation
    void class_get_lensed_Cl_TT(std::vector<double>& result)
    {
      using namespace Pipes::class_get_lensed_Cl_TT;
      result = BEreq::class_get_lensed_cl("tt");

      // Loop through the TT spectrum and check if there is a negative value. If so, invalidate.
      for (auto it=result.begin(); it != result.end(); it++)
      {
        if (*it < 0.0)
        {
          invalid_point().raise("Found a negative value in the TT spectrum.");
        }
      }
    }

    /// Temperature & E-mode cross-correlation
    void class_get_lensed_Cl_TE(std::vector<double>& result)
    {
      using namespace Pipes::class_get_lensed_Cl_TE;
      result = BEreq::class_get_lensed_cl("te");
    }

    /// E-mode autocorrelation
    void class_get_lensed_Cl_EE(std::vector<double>& result)
    {
      using namespace Pipes::class_get_lensed_Cl_EE;
      result = BEreq::class_get_lensed_cl("ee");

      // Loop through the EE spectrum and check if there is a negative value. If so, invalidate.
      for (auto it=result.begin(); it != result.end(); it++)
      {
        if (*it < 0.0)
        {
          invalid_point().raise("Found a negative value in the EE spectrum.");
        }
      }
    }

    /// B-mode autocorrelation
    void class_get_lensed_Cl_BB(std::vector<double>& result)
    {
      using namespace Pipes::class_get_lensed_Cl_BB;
      result = BEreq::class_get_lensed_cl("bb");

      // Loop through the BB spectrum and check if there is a negative value. If so, invalidate.
      for (auto it=result.begin(); it != result.end(); it++)
      {
        if (*it < 0.0)
        {
          invalid_point().raise("Found a negative value in the BB spectrum.");
        }
      }
    }

    /// Lensing (Phi) autocorrelation
    void class_get_lensed_Cl_PhiPhi(std::vector<double>& result)
    {
      using namespace Pipes::class_get_lensed_Cl_PhiPhi;
      result = BEreq::class_get_lensed_cl("pp");

      // Loop through the PhiPhi spectrum and check if there is a negative value. If so, invalidate.
      for (auto it=result.begin(); it != result.end(); it++)
      {
        if (*it < 0.0)
        {
          invalid_point().raise("Found a negative value in the PhiPhi spectrum.");
        }
      }
    }

    ////////////////
    /// DarkAges ///
    ////////////////

    /// Get the energy injection efficiency table from DarkAges.
    void energy_injection_efficiency_func(DarkAges::Energy_injection_efficiency_table& result)
    {
      using namespace Pipes::energy_injection_efficiency_func;
      result = BEreq::get_energy_injection_efficiency_table();
    }

  } // namespace CosmoBit

} // namespace Gambit
