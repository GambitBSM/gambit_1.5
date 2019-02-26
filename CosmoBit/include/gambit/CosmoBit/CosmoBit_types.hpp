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
///  *********************************************


#ifndef __CosmoBit_types_hpp__
#define __CosmoBit_types_hpp__

#include "gambit/Backends/backend_types/class.hpp"
#include <valarray>

namespace Gambit
{

  namespace CosmoBit
  {
    typedef std::map< std::string,std::valarray < double > > map_str_valarray_dbl;
    
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
     /* Use these routines for times t before that but t >~ 10^3 s, where electrons and positrons are 
        already completely annihilated 
        For times after CMB release you can get these values from the class background structure.
     */
      public:
        SM_time_evo(double t0, double tf, double grid_size);
       
        // SM photon temperature (keV) as a function of time (seconds)  
        void set_T_evo() {T_evo = 1147.40/sqrt((get_t_grid()));}

        // SM neutrino temperature (keV) as a function of time (seconds)
        // Defined in this way, one gets Neff = 3.046, using Nnu = 3 SM neutrinos
        void set_Tnu_evo() {Tnu_evo = 822.0965/sqrt(t_grid);}
        
        // SM Hubble rate (1/s) as a function of time (seconds)
        void set_Ht_evo() {H_evo= 1/(2.0*t_grid);} 
        
        // SM Hubble rate (1/s) as a function of temperature (keV)
        void set_HT_evo(std::valarray<double> T) {H_evo = 3.7978719e-7*T*T;}

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
    };


    // Forward declaration of warnings and errors
    error& CosmoBit_error();
    warning& CosmoBit_warning();

    // map (dictionary) for name and value of the inputs for Class
    class ClassInput
    {
      public:
        //classInput();
        //~classInput();
        void addEntry(std::string key,std::string val);
        void addEntry(std::string key,double val);
        void addEntry(std::string key,int val);
        void addEntry(std::string key, std::vector<double> val);
        void addEntry(std::string key, std::vector<int> val);

        void clear();
        std::map<std::string,std::string> get_map();

      private:
	std::map<std::string,std::string> input_list;
    };

    // Container for the structs of Class
    class Class_container
    {
      public:
        Class_container();
        //~Class_container();

        ClassInput input;

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
  }
}

#endif // defined __CosmoBit_types_hpp__
