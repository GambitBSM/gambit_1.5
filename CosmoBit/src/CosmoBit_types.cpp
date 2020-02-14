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
///	 \date 2018 Jan
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 Oct
///  *********************************************

#include <string>
#include <iostream>
#include <stdlib.h>     /* malloc, free, rand */
#include <valarray>

#include "gambit/CosmoBit/CosmoBit_types.hpp"
#include "gambit/CosmoBit/CosmoBit_utils.hpp"
#include "gambit/Utils/numerical_constants.hpp"

namespace Gambit
{
  namespace CosmoBit
  {

    /// helper to convert the memory address a double pointer points to
    /// to an integer (-> uintptr_t, size of type depends on system & ensures 
    /// it is big enough to store memory addresses of the underlying setup)
    /// (implemented to pass contents of arrays to CLASS) 
    uintptr_t memaddress_to_uint(double* ptr)
    {
      uintptr_t addr;
      addr = reinterpret_cast<uintptr_t>(ptr);
      return addr; 
    }
    

    BBN_container::BBN_container()
    { 
      // maps elements to their position in 'ratioH' array in AlterBBN holding 
      // primordial element abundances relative to H abundance
      abund_map["H2"] = 3;
      abund_map["D"] = 3;
      abund_map["H3"] = 4;
      abund_map["He3"] = 5;
      abund_map["He4"] = 6;   
      abund_map["Yp"] = 6;
      abund_map["Li6"] = 7;
      abund_map["Li7"] = 8;
      abund_map["Be7"] = 9;
      abund_map["Li8"] = 10;
    }

    void BBN_container::init_arr_size(int nnuc)
    {
      NNUC = nnuc;
      BBN_abund.resize(NNUC+1, 0.);
      BBN_covmat.resize(NNUC+1, std::vector<double>(NNUC+1,0.));
    }
    
    MPLike_data_container::MPLike_data_container(){}
    MPLike_data_container::MPLike_data_container(pybind11::object &data_in, map_str_pyobj likelihoods_in): data(data_in), likelihoods(likelihoods_in){}

    SM_time_evo::SM_time_evo(double t0, double tf, double N_t) : grid_size(N_t), t_grid(N_t), T_evo(N_t), Tnu_evo(N_t), H_evo(N_t), H_int(N_t)
    {
      
      // check if implemented routines are valid for given initial time
      if(t0 < 1e3) 
      {
        std::ostringstream err;
        err << "Requested initial time for evolution of Temperature & Hubble rate for SM for Temperatures t_initial was "<< t0<<". Implemented routines are only valid after e+/- annihilation (t > 10^3). Aborting now.";
        CosmoBit_error().raise(LOCAL_INFO, err.str());
      }

      // initialize time grid in log space
      double Delta_logt = (log(tf) - log(t0))/(grid_size-1);      
      for (int jj = 0; jj<grid_size; ++jj) 
      {
           t_grid[jj] = exp(log(t0) + jj*Delta_logt);
      }
      double Neff_SM = 3.046;
      double g_star_SM = 2.+2.*7./8.*Neff_SM*pow(4./11.,4./3.); // contribution from photons & neutrinos with Neff = 3.046
      
      // factor needed to calculate temperature evolution. For details see definition of functions set_T_evo(),.. in CosmoBit_types.hpp header
      factor_T_evo = 1./sqrt(2.*sqrt(8.*pi*pi*pi*_GN_SI_ *pow(_kB_SI_,4.)*g_star_SM/90./_c_SI_/_c_SI_/pow(_hP_SI_/2./pi*_c_SI_,3.)))*_kB_eV_over_K_/1e3;
      // TODO: different from T_nu/T_gamma from nu oscillation results?
      factor_Tnu_evo = pow(Neff_SM/3.,1./4.)* pow(4./11.,1./3.)*factor_T_evo; // = T_nu/T_gamma * factor_T_evo
      factor_HT_evo = sqrt(8.*pi*_GN_SI_/3.*pi*pi/30.*g_star_SM/pow(_hP_SI_/2./pi*_c_SI_,3.))*(1e6*_eV_to_J_*_eV_to_J_)/_c_SI_;

      set_T_evo();
      set_Tnu_evo();
      set_Ht_evo();
    }


    void SM_time_evo::calc_H_int()
    {
      /*
      This calculates int(y(x') dx', {x', x0, x}) for each x in x_grid (with x0 = x_grid[0]), using a simple trapezoidal integration.
      This function explicitly assumes that the logarithms of the values in x_grid are equally spaced.
      This is very simplified and designed to work fast in case of the Hubble rate. Can go wrong with other functions, so it should
      really only be used in this context (or if you exactly know what you are doing..).
      */
      std::valarray<double> g_grid(grid_size);
      double Delta_z = (log(t_grid[grid_size-1]) - log(t_grid[0]))/(grid_size - 1);
      
      g_grid = t_grid*H_evo;
      
      H_int[0] = 0.0;
      H_int[1] = 0.5*Delta_z*(g_grid[0] + g_grid[1]);
      
      double cumsum = g_grid[1];
      for (int i = 2; i<grid_size; ++i) 
      {
          H_int[i] = 0.5*Delta_z*(g_grid[0] + g_grid[i] + 2.0*cumsum);
          cumsum = cumsum + g_grid[i];
      }
    }

/*
    ----------  ClassyInput Methods ---------
    In the following the methods of the Class 'ClassyInput' are implemented. 
    ClassyInput has the attribute 'input_dict' which is a python dictionary 
    containing the input parameters for CLASS

*/


    /// add all entries from extra_entries to input_dict, concatenates and returns all 
    /// keys that are contained in both dictionaries:
    /// -> no keys in common: returns empty string ("")
    /// -> else: returns sting containing all duplicated keys
    /// need to check after use of this function if returned string was empty to avoid overwriting of 
    /// input values & inconsistencies. 
    std::string ClassyInput::addDict(pybind11::dict extra_dict)
    {
      // string to be returned -- stays empty if no duplicated dict entries are found
      std::string common_keys ("");

      for (auto item : extra_dict)
      {
        pybind11::str key = pybind11::str(item.first);
        pybind11::str arg = pybind11::str(item.second);

        if(!input_dict.attr("has_key")(key).cast<bool>())
        {
          input_dict[key] = arg;
          //std::cout << "Adding key = " << std::string(pybind11::str(item.first)) << ", "<< "value=" << std::string(pybind11::str(item.second)) << std::endl;
        }
        else
        {
          // add duplicated strings 
          common_keys = common_keys + " " + key.cast<std::string>();
          //std::cout << "Found duplicated key = " << std::string(pybind11::str(item.first)) << ", "<< "value=" << std::string(pybind11::str(item.second)) << std::endl;
        }
      }
      return common_keys;
    }

    // function to merge python dictionary extra_dict into input_dict (member of ClassyInput). 
    void ClassyInput::merge_input_dicts(pybind11::dict extra_dict) 
    {
      static bool first_run = true;
      // call function to merge classy pybind dicts defined in CosmoBit_utils.cpp
      // this function takes care of merging the input. If a key is contained in both
      // dictionaries it has to be decided on a case-by-case basis how to deal with that
      // some examples are already implemented, might have to be extended in the future
      // if more CLASS features are used. 
      merge_pybind_dicts(input_dict, extra_dict, first_run);
      first_run = false;
    }

    // return stringstream to print current entries of 
    // input_dict (can be send to logger)
    std::string ClassyInput::print_entries_to_logger()
    {
      using namespace LogTags;
      std::ostringstream log_msg;
      log_msg << "Parameters passed to class: \n \n";
      for (auto item : input_dict)
      {  
        pybind11::str key = pybind11::str(item.first);
        pybind11::str value = pybind11::str(item.second);

        log_msg << std::string(key) << " = " << std::string(value) << " \n";
      }
      return log_msg.str();
    }

    // Default constructor for multimode inputs...
    multimode_inputs::multimode_inputs() {};

    // // Default constructor for primordial and parametrised power spectra
    // primordial_ps::primordial_ps() 
    // {
    //   std::cout << "Calling primordial_ps constructor..." << std::endl;
    // };

    void primordial_ps::fill_k(double *k_array, int len)
    {
      std::vector<double> K(k_array, k_array+len);
      k = std::move(K);
      vec_size = len;
      // (JR) the vector gets just filled with copies of kmin # todo
      // for testing -> atm 
      // issue fixed -- thanks to whoever did it :) let's leave the debug print
      // a while longer though, you never know .. 
      //for( int i =0; i<len;i++) {std::cout << k[i] << std::endl;};
    }

    void primordial_ps::fill_P_s(double *P_s_array, int len)
    {
      std::vector<double> ps(P_s_array, P_s_array+len);
      P_s = std::move(ps);
    }    
    void primordial_ps::fill_P_s_iso(double *P_s_iso_array, int len)
    {
      std::vector<double> psi(P_s_iso_array, P_s_iso_array+len);
      P_s_iso = std::move(psi);
    }    
    void primordial_ps::fill_P_t(double *P_t_array, int len)
    {
      std::vector<double> pt(P_t_array, P_t_array+len);
      P_t = std::move(pt);
    }


    parametrised_ps::parametrised_ps() {};

  }
}
