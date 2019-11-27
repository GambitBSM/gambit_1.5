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
///  *********************************************


#ifndef __CosmoBit_types_hpp__
#define __CosmoBit_types_hpp__

#include "gambit/Utils/util_types.hpp"
#include "gambit/Backends/backend_types/MontePythonLike.hpp"
#include <valarray>
#include <stdint.h>  // save memory address as int

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

namespace Gambit
{

  namespace CosmoBit
  {

    typedef std::map< std::string,std::valarray < double > > map_str_valarray_dbl;
    
    class MPLike_data_container
    {
    /* Class holing MPLIke data structure & map with initialised Likelihoods objects; this is
        separated form the Classy_cosmo_container since it needs to be initialised as 'static const'
        such that the initialisation and reading in of data only happens once. 
        This is essential since the parsing of the data at initialisation of a Likelihood object can take
        much longer than the actual Likelihood calculation. 
        --
        Memebers:
        --
        pybind11::object data : MPLike data structure
        map_str_pyobj likelihoods : map likelihood name to initialised MPLike likelihood object
    */
    public:

        MPLike_data_container();
        MPLike_data_container(pybind11::object &data_in, map_str_pyobj likelihoods_in);
        
        pybind11::object data;
        map_str_pyobj likelihoods;

    };

    class Classy_cosmo_container
    {
    /* Class holding python object cosmo which is an instance of the class Class() from classy (Python wrapper for CLASS)
        This needs to be passed around between the classy frontend, CosmoBit & the MPLike frontend. 
    */
    public:
        Classy_cosmo_container(){}
        pybind11::object cosmo;
        pybind11::dict cosmo_input_dict;
        pybind11::dict cosmo_prec_dict;

        void set_input_dict(pybind11::dict input_dict);
    };

    class BBN_container
    {
      public:
        BBN_container();

        std::vector<double> BBN_abund;
        std::vector< std::vector<double> > BBN_covmat;
        std::map<std::string, int> abund_map;

        void init_arr(int nnuc);
        int get_NNUC(){return NNUC;}; // global parameter in AlterBBN, holds number of computed element abundances
        std::map<std::string,int> get_map(){return abund_map;};

      private:
        int NNUC;
    };

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


    // Forward declaration of warnings and errors
    error& CosmoBit_error();
    warning& CosmoBit_warning();

/*
    ----------  ClassyInput Methods ---------
    In the following the methods of the Class 'ClassyInput' are implemented. 
    ClassyInput has the attribute 'input_dict' which is a python dictionary 
    containing the input parameters for CLASS

*/
    // Class that manages the input dictionary for classy
    class ClassyInput
    {
      public:
        // add all entries from extra_entries to input_dict, will throw an error if 
        // one entry is contained in both dictionaries
        // returns 1 if adding was successful, -1 if an entries 
        // appeared twice -> need to throw an error since overwriting CLASS input
        // without realising it can be dangerous // TODO: actually implement the check somewhere.. 
        int addDict(pybind11::dict extra_entries);

        void addEntry(str key, double value) {input_dict[key.c_str()]=value;};
        void addEntry(str key, int    value) {input_dict[key.c_str()]=value;};
        void addEntry(str key, str    value) {input_dict[key.c_str()]=value.c_str();};
        void addEntry(str key, std::vector<double>& values)
        {
            // get pointers to arrays holding the information that needs
            // to be passed on to class, convert to uintptr_t (type large enough
            // to store memory address of the used system) and pass to class
            uintptr_t addr;
            addr = reinterpret_cast<uintptr_t>(&values[0]);
            input_dict[key.c_str()] = addr; 
        };

        bool hasKey(str key){return input_dict.contains(key.c_str());};
        //int addEntry(str key,std::ostringstream value){input_dict[key.c_str()]=value.c_str()};
        
        // merge dictionaries with overwriting/combining rules that only
        // apply for CLASS input dictionaries
        void merge_input_dicts(pybind11::dict extra_dict);
        
        // routine to print CLASS input values to logger
        std::string print_entries_to_logger();

        // clears all entries from input_dict
        void clear(){input_dict.attr("clear")();};

        // return input_dict
        pybind11::dict get_input_dict(){return input_dict;};

      private:
        pybind11::dict input_dict;
    };


    /// Class containing the primordial power spectrum.
    /// Members:
    /// - vector of modes k (1/Mpc)
    /// - scalar power spectrum of these modes P_s(k) (dimensionless)
    /// - tensor power spectrum of these modes P_t(k) (dimensionless)
    class primordial_ps
    {
        public:
            primordial_ps();

            // Fill k from an array of doubles
            void fill_k(double *k_array, int len) 
            {
                std::vector<double> K(k_array, k_array+len);
                k = std::move(K);
            };
            void fill_P_s(double *P_s_array, int len)
            {
                std::vector<double> ps(P_s_array, P_s_array+len);
                P_s = std::move(ps);
            };
            void fill_P_s_iso(double *P_s_iso_array, int len)
            {
                std::vector<double> ps_iso(P_s_iso_array, P_s_iso_array+len);
                P_s_iso = std::move(ps_iso);
            };
            void fill_P_t(double *P_t_array, int len)
            {
                std::vector<double> pt(P_t_array, P_t_array+len);
                P_t = std::move(pt);
            };

            std::vector<double>& get_k() { return k; }
            std::vector<double>& get_P_s() { return P_s; }
            std::vector<double>& get_P_t() { return P_t; }

        private:
            std::vector<double> k;
            std::vector<double> P_s;
            std::vector<double> P_s_iso;
            std::vector<double> P_t;
    };

    /// Class containing the *parametrised* primordial power spectrum.
    /// Members:
    /// - spectral tilt n_s
    /// - amplitude of scalar perturbations A_s
    /// - scalar to tensor ratio r
    class parametrised_ps
    {
        public:
            parametrised_ps();

            void set_ns(double ns) {n_s = ns;};
            void set_As(double As) {A_s = As;};
            void set_r(double R) {r = R;};

            double get_ns() {return n_s;};
            double get_As() {return A_s;};
            double get_r() {return r;};


        private:
            double n_s;
            double A_s;
            double r;
    };
  }
}

#endif // defined __CosmoBit_types_hpp__
