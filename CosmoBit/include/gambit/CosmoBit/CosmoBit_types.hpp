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
#include "gambit/CosmoBit/CosmoBit_utils.hpp"
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

    /// helper to convert the memory address a double pointer points to
    /// to an integer (-> uintptr_t, size of type depends on system & ensures
    /// it is big enough to store memory addresses of the underlying setup)
    /// (implemented to pass contents of arrays to CLASS)
    uintptr_t memaddress_to_uint(double* ptr);

    typedef std::map< std::string,std::valarray < double > > map_str_valarray_dbl;

    /// type to store likelihood calculation results from MontePython
    /// has two members: logLike_results & obs_results, both string to double map
    /// mapping the MontePython experiment/likelihood name to the calculated LogLike
    ///     - logLike_results:  contains results from all experiment names tagged with 'Likelihood',
    ///         i.e. the likelihoods being added to total LogLike to drive scan
    ///     - obs_results:  contains results from all experiment names tagged with 'Observable',
    ///         i.e. the likelihoods NOT added to total LogLike to drive scan
    class MPLike_result_container
    {
    public:

        MPLike_result_container(){};

        // add map entry to likelihood/observables map
        void add_logLike(std::string experiment_name, double logLike) {logLike_results[experiment_name] = logLike;}
        void add_obs(std::string experiment_name, double obs) {obs_results[experiment_name] = obs;}

        // return likelihood/observable result mpas
        map_str_dbl get_logLike_results () {return logLike_results;}
        map_str_dbl get_obs_results     () {return obs_results;}

    private:
        map_str_dbl logLike_results;
        map_str_dbl obs_results;

    };

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


    /// Class to store all results from an AlterBBN run
    /// -> element abundances stored in BBN_nuc (length NNUC+1), 
    /// -> covariance matrix in BBN_covmat ( dim NNUC+1 x NNUC+1)
    /// -> abund_map maps name of element to position in BBN_abundance vector
    ///    see constructor of BBN_container
    class BBN_container
    {
      public:
        BBN_container();

        /// TODO make these private and check that it does not break anything.. (I think it does -> add setter and getter functions)
        std::vector<double> BBN_abund;
        std::vector< std::vector<double> > BBN_covmat;
        std::map<std::string, int> abund_map; // TODO make map return a BE convinience function (just in case something chances in AlterBBN)

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

        /// add all entries from extra_entries to input_dict, concatenates and returns all
        /// keys that are contained in both dictionaries:
        /// -> no keys in common: returns empty string ("")
        /// -> else: returns sting containing all duplicated keys
        /// need to check after use of this function if returned string was empty to avoid overwriting of
        /// input values & inconsistencies.
        std::string addDict(pybind11::dict extra_entries);

        void addEntry(str key, double value) {input_dict[key.c_str()]=std::to_string(value).c_str();};
        void addEntry(str key, int    value) {input_dict[key.c_str()]=std::to_string(value).c_str();};
        void addEntry(str key, str    value) {input_dict[key.c_str()]=value.c_str();};
        void addEntry(str key, std::vector<double>& values)
        {
            // get pointers to arrays holding the information that needs
            // to be passed on to class, convert to uintptr_t (type large enough
            // to store memory address of the used system) and pass to class
            //uintptr_t addr;
            //addr = reinterpret_cast<uintptr_t>(&values[0]);
            input_dict[key.c_str()] = memaddress_to_uint(&values[0]);

        };
        //void addEntry(str key, double* ptr)
        //{
        //    // get pointers to arrays holding the information that needs
        //    // to be passed on to class, convert to uintptr_t (type large enough
        //    // to store memory address of the used system) and pass to class
        //    uintptr_t addr;
        //    addr = reinterpret_cast<uintptr_t>(ptr);
        //    input_dict[key.c_str()] = addr;
        //};

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

    /// Class containing the inputs used for inputs to MultiModeCode
    class multimode_inputs
    {
        public:
            multimode_inputs();

            // K array.
            double k_min;
            double k_max;
            int numsteps;

            double N_pivot;

            // Potential parameters
            std::vector<double> vparams;
            std::vector<double> phi_init0;
            std::vector<double> dphi_init0; // TODO probably remove this
            int num_inflaton = -1;
            int potential_choice = -1;
            int vparam_rows = -1;

            // Initial conditions stuff
            int slowroll_infl_end;
            int instreheat;
            int ic_sampling;
            double energy_scale;
            int numb_samples;
            int save_iso_N;
            double N_iso_ref;

            // Analytic approximations
            int use_deltaN_SR;
            int evaluate_modes;
            int use_horiz_cross_approx;
            int get_runningofrunning;

            double k_pivot;
            double dlnk;

    };


    /// Class containing the primordial power spectrum.
    /// Members:
    /// - vector of modes k (1/Mpc)
    /// - scalar power spectrum of these modes P_s(k) (dimensionless)
    /// - tensor power spectrum of these modes P_t(k) (dimensionless)
    class primordial_ps
    {
        public:
            primordial_ps() {};
            ~primordial_ps() {};

            // Fill k from an array of doubles
            void fill_k(double*, int);
            void fill_P_s(double*, int);
            void fill_P_s_iso(double*, int);
            void fill_P_t(double*, int);
            void set_kpivot(double k) {k_pivot = k;};
            void set_Npivot(double N) {N_pivot = N;};


            std::vector<double>& get_k() { return k; }
            std::vector<double>& get_P_s() { return P_s; }
            std::vector<double>& get_P_t() { return P_t; }
            int get_vec_size() { return vec_size; }
            double get_kpivot() {return k_pivot;};
            double get_Npivot() {return N_pivot;};

        private:
            std::vector<double> k;
            std::vector<double> P_s;
            std::vector<double> P_s_iso;
            std::vector<double> P_t;
            // needed to pass vector length to CLASS,
            // set in 'fill_k' method
            int vec_size;
            // k_pivot used by CLASS in spectra module,
            // independent of primordial module, therefore
            // need to pass it consistently
            double k_pivot;
            double N_pivot; // called N_star in CLASS
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
            void set_kpivot(double k) {k_pivot = k;};
            void set_Npivot(double N) {N_pivot = N;};

            double get_ns() const {return n_s;};
            double get_As() const {return A_s;};
            double get_r() const {return r;};
            double get_kpivot() const {return k_pivot;};
            double get_Npivot() const {return N_pivot;};


        private:
            double n_s;
            double A_s;
            double r;
            double k_pivot;
            double N_pivot; // called N_star in CLASS
    };
  }
}

#endif // defined __CosmoBit_types_hpp__
