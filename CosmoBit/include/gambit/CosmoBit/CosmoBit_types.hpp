//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Type definition header for module CosmoBit.
///
///  Compile-time registration of type definitions
///  required for the rest of the code to
///  communicate with CosmoBit.
///
///  Add to this if you want to define a new type
///  for the functions in CosmoBit to return, but
///  you don't expect that type to be needed by
///  any other modules.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2017 Nov
///
///  *********************************************


#ifndef __CosmoBit_types_hpp__
#define __CosmoBit_types_hpp__

#include "gambit/Backends/backend_types/class.hpp"

namespace Gambit
{

  namespace CosmoBit
  {

    // Forward declaration of warnings and errors
    error& CosmoBit_error();
    warning& CosmoBit_warning();

    // Container for the structs of Class
    struct Class_container
    {
      Class_container();
      ~Class_container();

      Class::file_content fc;     /* for input parameters */
      Class::precision pr;        /* for precision parameters */
      Class::background ba;       /* for cosmological background */
      Class::thermo th;           /* for thermodynamics */
      Class::perturbs pt;         /* for source functions */
      Class::transfers tr;        /* for transfer functions */
      Class::primordial pm;       /* for primordial spectra */
      Class::spectra sp;          /* for output spectra */
      Class::nonlinear nl;        /* for non-linear spectra */
      Class::lensing le;          /* for lensed spectra */
      Class::output op;           /* for output files */
      Class::ErrorMsg class_errmsg;      /* for error messages */

      bool non_free_pointer;

      int lmax;
      std::vector<double> Cl_TT;
      std::vector<double> Cl_TE;
      std::vector<double> Cl_EE;
      std::vector<double> Cl_BB;
      std::vector<double> Cl_PhiPhi;
		
	  std::vector<double> Pk_S; // Primordial Scalar Power Spectrum
	  std::vector<double> Pk_T; // Primordial Tensor Power Spectrum
	  std::vector<double> k_ar; // Corresponding wavenumbers.

    };

    // Generic class for cosmological likelihoods
/*    class CosmoLike
    {
      public:
        double get_l_max() const;
        bool needs_TT() const;
        bool needs_Polarization() const;
        bool needs_lensing() const;
    }; */
  }
}

#endif // defined __CosmoBit_types_hpp__
