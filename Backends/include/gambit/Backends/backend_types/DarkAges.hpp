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
  // Shorthand for map from string to pointer to double (here atm, should probably move into classy or DarkAges types) TODO
  typedef std::map<std::string,double*> map_str_dblptr;
  
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
      // Added for classy interface to exoclass
      int num_lines;
      map_str_dblptr ptrs_to_member_vecs; 
      
    };
  }
}

#endif // defined __DarkAges_types_hpp__
