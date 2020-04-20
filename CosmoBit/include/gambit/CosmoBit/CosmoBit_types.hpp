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
///  \date 2018 May
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 Oct
///  \date 2019 Mar
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


#ifndef __CosmoBit_types_hpp__
#define __CosmoBit_types_hpp__

#include "gambit/Utils/util_types.hpp"
#include "gambit/CosmoBit/CosmoBit_utils.hpp"
#include "gambit/Backends/backend_types/MontePythonLike.hpp"
#include <valarray>
#include <tuple>

#include <pybind11/stl.h>

namespace Gambit
{

  namespace CosmoBit
  {

    // Forward declaration of warnings and errors
    error& CosmoBit_error();
    warning& CosmoBit_warning();

    typedef std::map< str,std::valarray < double > > map_str_valarray_dbl;
    typedef std::tuple<pybind11::object, map_str_str, map_str_pyobj> MPLike_objects_container;

    class SM_time_evo
    {
     /* Time evolution of photon & neutrino temperature as well as Hubble rate in SM for time t > 1e3 s.
        Heads-up: explicit SM-specific assumptions enter. I.e.
            - Neff = 3.046 -> T_nu/T_gamma = (Neff/Nnu)^(1/4) (4/11)^(1/4) = 0.716486
            - Only contribution to Omega_rad in early Universe are photons & neutrinos
            - g_star (relativistic degrees of freedom) is constant in time
        => Use these routines only for times t >~ 10^3 s, where electrons and positrons are
        already completely annihilated
        For times after CMB release you can get these values from the class background structure.
     */
      public:
        SM_time_evo(double t0, double tf, double grid_size);

        // SM photon temperature (keV) as a function of time (seconds)
        // T_g [keV] = T_cmb [keV] / sqrt(2*H0*(Omega_r0)^0.5
        //           = { [2 sqrt(8pi^3 G[_SI_] (kB [J/K])^4] *g_star_SM / [90 c^2 (c hbar [SI])^3] }^-0.5
        // with G: Newton's constant in m^3/kg/s^2, kB: Boltzmann const in J/K, c: speed of light in m/s,
        // hbar x c: reduced Planck const in Jm
        // H0 and T_cmb cancel out in SM where Omega_r0 = [1 + Neff] Omega_g0, the factor that stays besides natural constants is
        // g_star = g_gamma + 7/8 * g_nu * Neff * (4/11)^(4/3) = 3.38354
        void set_T_evo() {T_evo = factor_T_evo/sqrt(t_grid);}

        // SM neutrino temperature (keV) as a function of time (seconds)
        // Defined in this way, one gets Neff = 3.046, using Nnu = 3 SM neutrinos (Tnu = 0.716486 T(t0))
        void set_Tnu_evo() {Tnu_evo = factor_Tnu_evo/sqrt(t_grid);}

        // SM Hubble rate (1/s) as a function of time (seconds)
        // Solve Friedmann eq. for rad. dominated Universe -> H(t) = (da/dt)/a = H0 sqrt(Omega_g0 a^-4)
        void set_Ht_evo() {H_evo= 1./(2.0*t_grid);}

        // SM Hubble rate (1/s) as a function of temperature (keV)
        // H^2 = 8 pi G/3 rho = 8 pi G/3  g_star pi^2/30 (k_b T)^4/(hbar c)^3
        void set_HT_evo(std::valarray<double> T) {H_evo = factor_HT_evo*T*T;}

        // fast integration of Hubble rate using trapezoidal integration
        void calc_H_int();

        std::valarray<double> get_t_grid() const {return t_grid;}
        std::valarray<double> get_T_evo() const {return T_evo;}
        std::valarray<double> get_Tnu_evo() const {return Tnu_evo;}
        std::valarray<double> get_H_evo() const {return H_evo;}
        std::valarray<double> get_H_int() const {return H_int;}

    private:
        int grid_size;
        std::valarray<double> t_grid;
        std::valarray<double> T_evo;
        std::valarray<double> Tnu_evo;
        std::valarray<double> H_evo;
        std::valarray<double> H_int;

        double factor_T_evo;
        double factor_Tnu_evo;
        double factor_HT_evo;
    };

    /// Class containing the inputs used for inputs to MultiModeCode
    class Multimode_inputs
    {
        public:
            // Constructor
            Multimode_inputs();
            // Debugging options
            int silence_output;
            // k values where to evaluate the power spectrum
            double k_min;
            double k_max;
            int numsteps;
            // Parameters realted to the pivot scale
            double k_pivot;
            double N_pivot;
            double dlnk;
            // Parameters related to the potential and initial condidtions
            int num_inflaton = -1;
            int potential_choice = -1;
            int vparam_rows = -1;
            std::vector<double> vparams;
            std::vector<double> phi_init0;
            std::vector<double> dphi_init0;
            // Parameters realted to the scenario for initial conditions
            int slowroll_infl_end;
            int instreheat;
            // Parameters related to approximations and observables
            int use_deltaN_SR;
            int evaluate_modes;
            int use_horiz_cross_approx;
            int get_runningofrunning;
    };


    /// Class containing the primordial power spectrum.
    /// Members:
    /// - vector of modes k (1/Mpc)
    /// - scalar power spectrum of these modes P_s(k) (dimensionless)
    /// - tensor power spectrum of these modes P_t(k) (dimensionless)
    class Primordial_ps
    {
        public:
            Primordial_ps() {};
            ~Primordial_ps() {};

            // Fill k from an array of doubles
            void set_N_pivot(double npiv) { N_pivot = npiv; }
            void fill_k(double*, int);
            void fill_P_s(double*, int);
            void fill_P_s_iso(double*, int);
            void fill_P_t(double*, int);

            double get_N_pivot() { return N_pivot; }
            std::vector<double>& get_k() { return k; }
            std::vector<double>& get_P_s() { return P_s; }
            std::vector<double>& get_P_t() { return P_t; }
            int get_vec_size() { return vec_size; }

        private:
            double N_pivot;
            std::vector<double> k;
            std::vector<double> P_s;
            std::vector<double> P_s_iso;
            std::vector<double> P_t;
            // needed to pass vector length to CLASS,
            // set in 'fill_k' method
            int vec_size;
    };

    /// Class containing the *parametrised* primordial power spectrum.
    /// Members:
    /// - spectral tilt n_s
    /// - amplitude of scalar perturbations A_s
    /// - scalar to tensor ratio r
    class Parametrised_ps
    {
        public:
            Parametrised_ps() {};
            ~Parametrised_ps() {};

            void set_N_pivot(double npiv) { N_pivot = npiv; }
            void set_n_s(double ns) { n_s = ns; }
            void set_ln10A_s(double ln10As) { ln10A_s = ln10As; }
            void set_r(double rr) { r = rr; }

            double get_N_pivot() { return N_pivot; }
            double get_n_s() { return n_s; }
            double get_ln10A_s() { return ln10A_s; }
            double get_r() { return r; }

            // return members as str to double map for printing
            map_str_dbl get_parametrised_ps_map();


        private:
            double N_pivot;
            double n_s;
            double ln10A_s;
            double r;
    };
  }
}

#endif // defined __CosmoBit_types_hpp__
