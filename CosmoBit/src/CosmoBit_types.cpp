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
#include "gambit/Utils/numerical_constants.hpp"

namespace Gambit
{
  namespace CosmoBit
  {

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

    void BBN_container::init_arr(int nnuc)
    {
      NNUC = nnuc;
      BBN_abund.resize(NNUC+1, 0.);
      BBN_covmat.resize(NNUC+1, std::vector<double>(NNUC+1,0.));
    }
    
    MPLike_data_container::MPLike_data_container(){}
    MPLike_data_container::MPLike_data_container(pybind11::object &data_in, map_str_pyobj likelihoods_in): data(data_in), likelihoods(likelihoods_in){}

    void Classy_cosmo_container::set_input_dict(pybind11::dict input_dict) {cosmo_input_dict = input_dict;}

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


    Class_container::Class_container() : lmax(2508)
    {
      input.clear();
      Cl_TT.resize(lmax+1, 0.);
      Cl_TE.resize(lmax+1, 0.);
      Cl_EE.resize(lmax+1, 0.);
      Cl_BB.resize(lmax+1, 0.);
      Cl_PhiPhi.resize(lmax+1, 0.);
      //std::cout << "Hello it's me. I am a Class_container" << std::endl;

      Pk_S.resize(lmax+1, 0.); // Primordial Scalar Power Spectrum
      Pk_T.resize(lmax+1, 0.); // Primordial Tensor Power Spectrum
      k_ar.resize(lmax+1, 0.); // Corresponding wavenumbers.
    }
/*
    Class_container::~Class_container()
    {
      //std::cout << "And I am out" << std::endl;
    }
*/
    void ClassInput::addEntry(std::string key, std::string val)
    {
      input_list[key] = val;
    }

    void ClassInput::addEntry(std::string key, double val)
    {
      std::ostringstream stringified_val;
      stringified_val << val;
      addEntry(key,stringified_val.str());
    }

    void ClassInput::addEntry(std::string key, int val)
    {
      std::ostringstream stringified_val;
      stringified_val << val;
      addEntry(key,stringified_val.str());
    }

    void ClassInput::addEntry(std::string key, std::vector<double> val)
    {
      std::ostringstream stringified_val;
      bool first = true;
      for (auto it = val.begin(); it != val.end(); it++)
      {
	if (not first)
	  stringified_val << ", ";
	stringified_val << *it;
	first = false;
      }
      addEntry(key,stringified_val.str());
    }

    void ClassInput::addEntry(std::string key, std::vector<int> val)
    {
      std::ostringstream stringified_val;
      bool first = true;
      for (auto it = val.begin(); it != val.end(); it++)
      {
	if (not first)
	  stringified_val << ", ";
	stringified_val << *it;
	first = false;
      }
      addEntry(key,stringified_val.str());
    }


    std::string ClassInput::print_entries_to_logger()
    {
      using namespace LogTags;
      std::ostringstream log_msg;
      log_msg << "Parameters passed to class: \n \n";
      std::map<std::string,std::string> in_map = get_map();
      for(auto iter=in_map.begin(); iter != in_map.end(); iter++)
      {
        log_msg << iter->first.c_str() << " = " << iter->second.c_str() << " \n"; // uncomment if you want class input to be printed in terminal before class call 
      }
      return log_msg.str();
    }

    void ClassInput::clear()
    {
      input_list.clear();
    }

    std::map<std::string,std::string> ClassInput::get_map()
    {
      return input_list;
    }

/*
    ----------  ClassyInput Methods ---------
    In the following the methods of the Class 'ClassyInput' are implemented. 
    ClassyInput has the attribute 'input_dict' which is a python dictionary 
    containing the input parameters for CLASS

*/


    // add all entries from extra_dict to input_dict, will throw an error if 
    // one entry is contained in both dictionaries
    // returns 1 if adding was successful, -1 if an entries 
    // appeared twice -> need to throw an error since overwriting CLASS input
    // without realising it can be dangerous
    int ClassyInput::addDict(pybind11::dict extra_dict)
    {
      int success = 1;
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
          success = -1;
          break;
        }
      }
      return success;
    }

    // function to merge python dictionary extra_dict into input_dict. If both dictionaries have the same key
    // the values of the keys will be concatenated (without duplicating an entry)
    // -> very specific for merging the classy input dictionaries
    // typical example for this would be 
    // dict input_dict: 'output' : 'tCl nCl' and 
    // dict extra_dict: 'output' : 'tCl mPk' => mPk has to be added such that 'output' : 'tCl nCl mPk'
    void ClassyInput::merge_input_dicts(pybind11::dict extra_dict) 
    {
      // loop through extra_dict (should typically be the shorter one)
      static bool first_run = true;
      for (auto item : extra_dict)
      {  
        pybind11::str key = pybind11::str(item.first);
        pybind11::str arg = pybind11::str(item.second);
        
        // if item not contained in extra_dict but not in input_dict it will be added to input_dict
        if(!input_dict.attr("has_key")(key).cast<bool>())
        {
          input_dict[key] = arg;
          //std::cout << "Adding key = " << std::string(pybind11::str(item.first)) << ", "<< "value=" << std::string(pybind11::str(item.second)) << std::endl;
        }

        // if item contained in both: have to decide on a case-by-case basis what to do! 
        else
        { 
          // if 'output' is defined twice the entries have to be merged
          // e.g. input_dict['output'] = 'A B C'
          //      extra_dict['output'] = 'A B X Y'
          // should result in input_dict['output'] = 'A B C X Y'
          if(key.attr("find")(pybind11::str("output")).cast<int>()!=-1) // (python string.find("x") returns -1 if "x" not contained)
          {
            // split string of extra_dict['output'] by spaces into list 
            // (-> in the example above this would give list = ['A', 'B', 'X', 'Y'])
            pybind11::list list = extra_dict[key].attr("split")();
            // iterate through list and check if current item is also contained in input_dict['output']
            // if it is not it will be added
            for(auto it : list)
            { 
              std::string list_entry = std::string(pybind11::str(it));

              // add entry if it is not already in input_dict['output']
              if(input_dict[key].attr("find")(pybind11::str(list_entry)).cast<int>()==-1)
              { 
                // add part of extra_dict[key] string that is not contained in input_dict[key] string to input_dict[key]
                std::string new_arg=std::string(pybind11::str(input_dict[key]))+std::string(" ")+std::string(list_entry);
                input_dict[key]= new_arg;
              }
            }
          }

          // if 'l_max...' (scalars or tensors) is defined twice use the maximum requested value
          // e.g. input_dict['l_max_scalars'] = '500'
          //      extra_dict['l_max_scalars'] = '2500'
          // should result in input_dict['l_max_scalars'] = '2500'
          else if(key.attr("find")(pybind11::str("l_max")).cast<int>()!=-1)
          {
            // cast pybind11::detail::item_accessor to pybind11::str, to c++ string and then to int
            //  (I know... the problem is that the yaml file entries are all parsed as strings so 
            //  we don't have to distinguish between the different types of different CLASS input keys. 
            //  The python wrapper parses everything as string so it is fine. Only here, where we actually
            //  have to compare the passed numbers this causes troubles and we have to cast from a python 
            //  string to a c++ int to be able to tell which entry is the larger one. If you have an idea
            //  to make this nicer --> go for it!!)
            int lmax_input_dict = std::stoi(std::string(pybind11::str(input_dict[key])));
            int lmax_extra_dict = std::stoi(std::string(pybind11::str(extra_dict[key])));

            // if lmax_extra_dict is higher than the entry in the input_dict, replace is
            if (lmax_input_dict < lmax_extra_dict){input_dict[key] = lmax_extra_dict;}
            // if not 'input_dict'already contains the higher value and there is nothing to do here. 
          }
          else
          {
            if(first_run)
            {
              // (JR) TODO see what other combinations could go wrong here -> different likelihoods
              // asking for different redshift bins?
              // But need this check only on the first run (the inputs we are worried about are 
              // parameters specifying which output CLASS should produce so they won't change 
              // in-between the calculation of different points in one scan)
              std::cout <<"___________________________________________________________________"<<std::endl;
              std::cout <<"___________________________________________________________________"<<std::endl;
              std::cout << "Both dictionaries to merge contain key" << std::string(key) << "with entries "
                  << std::string(pybind11::str(input_dict[key]))<< " and " << std::string(pybind11::str(extra_dict[key])) << ". Don't know how to deal with that, yet. Will be taken care of soon." << std::endl;
              std::cout <<"___________________________________________________________________"<<std::endl;
              std::cout <<"___________________________________________________________________"<<std::endl;
            }
          }
        }
      }
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

        log_msg << std::string(key) << " = " << std::string(value) << " \n"; // uncomment if you want class input to be printed in terminal before class call 
      }
      return log_msg.str();
    }

    primordial_ps::primordial_ps(){};
    parametrised_ps::parametrised_ps(){};

  }
}
