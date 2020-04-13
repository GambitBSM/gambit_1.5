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

    /// Constructor
    BBN_container::BBN_container() : abund_map{{"H2",3}, {"D",3}, {"H3",4}, {"He3",5}, {"He4",6}, {"Yp",6}, {"Li6",7}, {"Li7",8}, {"Be7",9}, {"Li8",10}}
    {}

    /// Initialize sizes of vectors (get NNUC, number of computed element abundances, from AlterBBN)
    void BBN_container::init_arr_size(int nnuc)
    {
      NNUC = nnuc;
      BBN_abund.resize(NNUC+1, 0.);
      BBN_covmat.resize(NNUC+1, std::vector<double>(NNUC+1,0.));
    }

    /// Initialise the translation map from element name to position in abundance vector
    void BBN_container::set_abund_map(map_str_int map_in) {abund_map = map_in;}

    /// Setter functions for abundance vector
    void BBN_container::set_BBN_abund(int pos, double val) {BBN_abund[pos] = val;}

    /// Setter function for covariance matrix
    void BBN_container::set_BBN_covmat(int row, int col, double val) {BBN_covmat[row][col] = val;}

    /// Global parameter in AlterBBN; holds number of computed element abundances
    int BBN_container::get_NNUC() {return NNUC;};

    /// Getter for map from isotope names to position in BBN_abundance vector
    const std::map<std::string,int>& BBN_container::get_abund_map() {return abund_map;};

    /// Getter for abundance vector
    double BBN_container::get_BBN_abund(int pos) {return BBN_abund[pos];}

    /// Getter for covariance matrix
    double BBN_container::get_BBN_covmat(int row, int col) {return BBN_covmat[row][col];}

    /// Setter for active isotopes
    void BBN_container::set_active_isotopes(std::set<str> isos)
    {
      active_isotopes = isos;
      active_isotope_indices.clear();
      for (const str& s : active_isotopes) active_isotope_indices.insert(abund_map.at(s));
    }

    /// Getter for active isotopes
    const std::set<str>& BBN_container::get_active_isotopes() {return active_isotopes;}

    /// Getter for indices of active isotopes in BBN_abundance vector
    const std::set<int>& BBN_container::get_active_isotope_indices() {return active_isotope_indices;}


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
      // (JR) the vector gets just filled with copies of kmin # todo
      // for testing -> atm
      // issue fixed -- thanks to whoever did it :) let's leave the debug print
      // a while longer though, you never know ..
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
