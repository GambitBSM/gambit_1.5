//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Function definitions of CosmoBit.
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
///  *********************************************
#include <cmath>
#include <functional>
#include <iostream>
#include <omp.h>
#include <stdlib.h>     /* malloc, free, rand */
#include <string>
#include <valarray>
#include <memory>  // make_unique pointers
#include <stdint.h> // save memory addresses as int
#include <boost/algorithm/string/trim.hpp>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_spline.h>

#include "gambit/Utils/yaml_options.hpp"
#include "gambit/Utils/statistics.hpp"
#include "gambit/Utils/ascii_table_reader.hpp"
#include "gambit/Utils/ascii_dict_reader.hpp"
#include "gambit/Utils/numerical_constants.hpp"
#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/CosmoBit/CosmoBit_rollcall.hpp"
#include "gambit/CosmoBit/CosmoBit_types.hpp"
#include "gambit/CosmoBit/CosmoBit_utils.hpp"

/// Helper function for diagnosing MultiModeCode errors
std::string multimode_error_handling(int& err)
{

  std::string message = "MultiModeCode error: ";
  switch(err)
  {

    /// > 0 = "failure; not fatal"
    case 1:
      message = "Inflation did not start.";
      break;
    case 2:
      message = "The pivot scale didn't leave the horizon.";
      break;
    case 3:
      message = "A modes' IC couldn't be consistently set.";
      break;
    case 4:
      message = "Too many e-folds; can't initialize the scale factor.";
      break;
    case 5:
      message = "Trying to save the field values at some reference N, but got less evolution than that.";
      break;
    case 6:
      message = "Didn't satisfy reheating bounds.";
      break;


    /// < 0 = "fatal"
    case -1: 
      message = "Numerical underflow error in odeint.";
      break;


    // Specifically added for SMASH
    case 999:
      message = "ICs for SMASH potential resulted in too long(or short) inflation.";
      break;


    // Otherwise -- who knows.
    default: 
      message = "GAMBIT caught an error in MultiMode. Check the MultiModeCode output for more info.";
  }
  return message;
}

namespace Gambit
{

  namespace CosmoBit
  {
    using namespace LogTags;

    struct my_f_params { double a; double b; double c; double d; };

    // Utility functions (not rolcalled)
    double fn1 (double x, void * p)
    {
      struct my_f_params * params
      = (struct my_f_params *)p;
      double a = (params->a);
      double b = (params->b);
      double c = (params->c);
      double d = (params->d);

      return exp(4.0/3.0*(d-0.125*(pow(10.0,b)+6.0*pow(10.0,a))*(-c*c+x*x)))
      -(1.0+pow(10.0,a)*c*c)/(1.0+pow(10.0,a)*x*x);
    }


    // Returns SMASH potential value given parameters and field amplitude
    double pot_SMASH(double lambda, double phi, double xi)
    {
      double Mpl = 1.0;
      double pot;

      pot = lambda/4.0*pow(phi,4.0)*pow(1.0+xi*(pow(phi,2.0))/Mpl,-2.0);

      return pot;
    }

    double SRparameters_epsilon_SMASH(double phi, double b, double xi)
    {
      double Mpl = 1.0;
      double eps = 0.0;

      eps = 8.0*pow(Mpl,4.0)/(b*phi*phi*Mpl*Mpl+xi*(b+6.0*xi)*pow(phi,4.0));

      std::cout << "eps = " << eps << std::endl;

      return eps;
    }

    double SRparameters_eta_SMASH(double phi, double b, double xi)
    {
      double Mpl = 1.0;
      double eta = 0.0;

      eta = ((12.0*b*pow(Mpl,6.0)+4.0*pow(Mpl,4.0)*xi*phi*phi*(b+12.0*xi)
      -8.0*Mpl*Mpl*xi*xi*pow(phi,4.0)*(b+6.0*xi))
      /pow(b*Mpl*Mpl*phi+xi*pow(phi,3.0)*(b+6.0*xi),2.0));

      std::cout << "eta = " << eta << std::endl;

      return eta;
    }

    double ns_SR(double eps,double eta)
    {
      double ns = 0.0;

      ns = 1.0 + 2.0*eta - 6.0*eps;

      return ns;
    }

    double r_SR(double eps)
    {
      double r = 0.0;

      r = 16.0*eps;

      return r;
    }

    double As_SR(double eps,double pot)
    {
      const double pi = std::acos(-1);
      double Mpl = 1.0;
      double As = 0.0;

      As = pot/24.0/pi/pi/eps/pow(Mpl,4.0);

      return As;
    }

    void injection_spectrum_annihilatingDM(DarkAges::injectionSpectrum& spectrum)
    {
      using namespace Pipes::injection_spectrum_annihilatingDM;

      double m = *Param["mass"];
      double BR_el = *Param["BR"];
      double BR_ph = 1.0 - BR_el;

      if (m <= m_electron && BR_el >= std::numeric_limits<double>::epsilon())
      {
        std::ostringstream err;
        err << "The mass of the annihilating dark matter candiate is below the electron mass.";
        err << " No production of e+/e- is possible.";
        CosmoBit_error().raise(LOCAL_INFO,err.str());
      }

      spectrum.E_el.clear();
      spectrum.E_ph.clear();
      spectrum.spec_el.clear();
      spectrum.spec_ph.clear();

      spectrum.E_el.resize(1,std::max(m-m_electron, std::numeric_limits<double>::min()));
      spectrum.E_ph.resize(1,m);
      spectrum.spec_el.resize(1,BR_el*2e9);
      spectrum.spec_ph.resize(1,BR_ph*2e9);
    }

    void injection_spectrum_decayingDM(DarkAges::injectionSpectrum& spectrum)
    {
      using namespace Pipes::injection_spectrum_decayingDM;

      double m = *Param["mass"];
      double BR_el = *Param["BR"];
      double BR_ph = 1.0 - BR_el;

      if (m <= 2*m_electron && BR_el >= std::numeric_limits<double>::epsilon())
      {
        std::ostringstream err;
        err << "The mass of the decaying dark matter candiate is below twice the electron mass.";
        err << " No production of e+/e- is possible.";
        CosmoBit_error().raise(LOCAL_INFO,err.str());
      }

      spectrum.E_el.clear();
      spectrum.E_ph.clear();
      spectrum.spec_el.clear();
      spectrum.spec_ph.clear();

      spectrum.E_el.resize(1,std::max(m*0.5-m_electron, std::numeric_limits<double>::min()));
      spectrum.E_ph.resize(1,m*0.5);
      spectrum.spec_el.resize(1,BR_el*2e9);
      spectrum.spec_ph.resize(1,BR_ph*2e9);
    }

    void lifetime_ALP_agg(double& result)
    {
      // lifetime in s if onlz the decay a -> g g is open.
      using namespace Pipes::lifetime_ALP_agg;

      double gagg = *Param["gagg"]; // in GeV^-1
      double ma = *Param["ma0"]; // in eV

      // Calculate the decay width (in GeV)
      // (It's maybe worth being a seperate capability)
      double Gamma = 1/64./pi * pow(gagg,2) * pow(ma,3) * 1e-27;

      // For the lifetime take 1/Gamma and translate GeV^-1 into s by multiplication with "hbar"
      result = 1./Gamma * hbar;

      // Reject points which have a lifetime bigger than 1e17s (or whatever the user chooses)
      // Gets only triggered if the user wishes to do so.
      // !! This is not a real physical bound but it is more to deal with the lack of likelihoods so far. !!
      static bool do_rejection = runOptions->getValueOrDef<bool>(false,"do_rejection");
      static double tdec_max = runOptions->getValueOrDef<double>(1.0e17,"reject_tau_bigger_than");
      if (do_rejection && result > tdec_max)
      {
        std::ostringstream err;
        err << "ALP lifetime (" << result << " [s]) exceeds the threshold of " << tdec_max <<" [s].";
        invalid_point().raise(err.str());
      }
    }

    // Compute the abundance of ALPs expected from thermal production via Primakoff processes
    void minimum_abundance_ALP(double& result)
    {
      using namespace Pipes::minimum_abundance_ALP;

      double gagg = *Param["gagg"];                 // Read axion-photon coupling in GeV^-1
      double T_R_in_GeV = 1e-3 * (*Param["T_R"]);   // Read reheating temperature in MeV and convert it to GeV

      // Check for stupid input (T_R < m_e) and throw an error if the user really pushed it that far.
      if (m_electron >= T_R_in_GeV)
        CosmoBit_error().raise(LOCAL_INFO,"The reheating temperature is below the electron mass.");

      result = 1.56e-5 * pow(gagg,2) * m_planck * (T_R_in_GeV - m_electron);

    }

    // Compute the minimal fraction of dark matter in ALPs expected from thermal production via Primakoff processes
    void minimum_fraction_ALP(double& result)
    {
      using namespace Pipes::minimum_fraction_ALP;

      const double rho0_crit_by_h2 = 3.*pow(m_planck_red*1e9,2) * pow((1e5*1e9*hbar/_Mpc_SI_),2); // rho0_crit/(h^2)

      double ma0 = *Param["ma0"];                    // non-thermal ALP mass in eV
      double T = *Param["T_cmb"];                        // CMB temperature in K
      double omega_cdm = *Param["omega_cdm"];        // omega_cdm = Omega_cdm * h^2

      // Consistency check: if the ALP abundance from all thermal processes is less than that expected just from Primakoff processes, invalidate this point.
      double Ya0_min = *Dep::minimum_abundance;
      double ssm0 = entropy_density_SM(T);           // SM entropy density today in eV^3 (cf. footnote 24 of PDG2018-Astrophysical parameters)
      double rho0_cdm = omega_cdm * rho0_crit_by_h2; // rho0_cdm = Omega_cdm * rho0_crit;
      double rho0_min = Ya0_min * ma0 * ssm0;        // energy density of axions today in eV^4
      result = rho0_min / rho0_cdm;
    }

    // The fraction of dark matter in decaying ALPs at the time of production
    void DM_fraction_ALP(double& result)
    {
      using namespace Pipes::DM_fraction_ALP;

      const double t_rec = 1e12;                     // Age of the Universe at recombination in s

      double f0_thermal = *Param["f0_thermal"];      // Fraction of DM in ALPs at production, due to thermal production
      double omega_ma = *Dep::RD_oh2;                // Cosmological density of ALPs from vacuum misalignment

      double omega_cdm = *Param["omega_cdm"];        // omega_cdm = Omega_cdm * h^2
      double tau_a = *Dep::lifetime;                 // ALP lifetime in s

      // Consistency check: if the ALP abundance from all thermal processes is less than that expected just from Primakoff processes, invalidate this point.
      double f0_min = *Dep::minimum_fraction;
      if (f0_thermal < f0_min)
      {
        std::ostringstream err;
        err << "The choice of f0_thermal (" << f0_thermal;
        err << ") is in contradiction with the minimum dark matter fraction f0_min (";
        err << f0_min << ") produced via Primakoff processes.";
        invalid_point().raise(err.str());
      }

      // Compute the total fraction of DM in ALPs at production
      double xi_ini = f0_thermal + omega_ma/omega_cdm;

      // Consistency check: invalidate if there are more ALPs than dark matter at the time of recombination (t ~ 1e12s)
      double xi_at_rec = xi_ini * exp(-t_rec/tau_a );
      if (xi_at_rec  > 1.)
      {
        std::ostringstream err;
        err << "ALPs are over-abundant (n_a > n_cdm) at t = 10^12 s. (n_a/n_cdm = "<< xi_at_rec <<")";
        invalid_point().raise(err.str());
      }

      result = xi_ini;
    }

    void total_DM_abundance_ALP(double& result)
    {
      using namespace Pipes::total_DM_abundance_ALP;

      const double rho0_crit_by_h2 = 3.*pow(m_planck_red*1e9,2) * pow((1e5*1e9*hbar/_Mpc_SI_),2); // rho0_crit/(h^2)
      double omega_cdm = *Param["omega_cdm"];  // omega_cdm = Omega_cdm * h^2
      double fraction = *Dep::DM_fraction;
      double rho0_ALP = omega_cdm * fraction* rho0_crit_by_h2; // rho0_cdm = Omega_cdm * rho0_crit

      double ma0 = *Param["ma0"];
      double TCMB = *Param["T_cmb"];
      double ssm0 = entropy_density_SM(TCMB);

      result = rho0_ALP / (ma0 * ssm0);
    }

    void energy_injection_efficiency_func(DarkAges::fz_table& result)
    {
      using namespace Pipes::energy_injection_efficiency_func;
      result = BEreq::DA_efficiency_function();
    }


    void f_effective_func(double& result)
    {
      using namespace Pipes::f_effective_func;

      bool silent = runOptions->getValueOrDef<bool>(false,"silent_mode");
      int last_steps = runOptions->getValueOrDef<int>(4,"show_last_steps");
      double z_eff = 0.01;
      if(ModelInUse("DecayingDM_general"))
      {
        z_eff = runOptions->getValueOrDef<double>(300.,"z_eff");
      }
      else if (ModelInUse("AnnihilatingDM_general"))
      {
        z_eff = runOptions->getValueOrDef<double>(600.,"z_eff");
      }

      DarkAges::fz_table fzt = *Dep::energy_injection_efficiency;

      bool f_eff_mode = fzt.f_eff_mode;
      std::vector<double> z = fzt.redshift;
      std::vector<double> fh = fzt.f_heat;
      std::vector<double> fly = fzt.f_lya;
      std::vector<double> fhi = fzt.f_hion;
      std::vector<double> fhei = fzt.f_heion;
      std::vector<double> flo = fzt.f_lowe;
      std::vector<double> feff = fzt.f_eff;

      int npts = z.size();
      double ftot[npts];
      double red[npts];
      for (int i = 0; i < npts; i++)
      {
        if (f_eff_mode)
          ftot[i] = feff.at(i);
        else
          ftot[i] = fh.at(i)+fly.at(i)+fhi.at(i)+fhei.at(i)+flo.at(i);
        red[i] = z.at(i);
      }

      gsl_interp_accel *gsl_accel_ptr = gsl_interp_accel_alloc();
      gsl_spline *spline_ptr = gsl_spline_alloc(gsl_interp_cspline, npts);

      gsl_spline_init(spline_ptr, red, ftot, npts);

      result = gsl_spline_eval(spline_ptr, z_eff, gsl_accel_ptr);

      gsl_spline_free(spline_ptr);
      gsl_interp_accel_free(gsl_accel_ptr);

      if (!silent)
      {
        std::cout << "################" << std::endl;
        if (ModelInUse("DecayingDM_general"))
        {
          std::cout << "Scenario: Dark Matter Decay" << std::endl;
          std::cout << "tau = " << *Param["lifetime"] << std::endl;
        }
        else if (ModelInUse("AnnihilatingDM_general"))
          std::cout << "Scenario: Dark matter (s-wave) annihilation" << std::endl;
        std::cout << "m = " << *Param["mass"] << std::endl;
        std::cout << "BR (electron) = " << *Param["BR"] << std::endl;
        std::cout << "---------------" << std::endl;
        if (f_eff_mode)
          std::cout << "z\tf_eff" << std::endl;
        else
          std::cout << "z\tf_heat\tf_lya\tf_hion\tf_heion\tf_lowe" << std::endl;
        for (unsigned int i = z.size() - last_steps; i < z.size(); i++)
        {
          if (f_eff_mode)
            std::cout << z.at(i) << "\t" << feff.at(i) << std::endl;
          else
            std::cout << z.at(i) << "\t" << fh.at(i) << "\t" << fly.at(i) << "\t" << fhi.at(i) << "\t" << fhei.at(i) << "\t" << flo.at(i)  << std::endl;
        }
        std::cout << "f_eff (sum of all channels at z = "<< z_eff << ") = " << result << std::endl;
        std::cout << "################\n" << std::endl;
      }
    }


    /// If no value for N_star (number of e folds before inflation ends) is specified in the yaml
    /// file AND no inflation model is scanned fall back to default value of CLASS; otherwise 
    /// it is a model parameter of the 
    /// (heads-up: name in CLASS is N_star, not N_pivot)
    void set_N_pivot(double &result)
    {

      using namespace Pipes::set_N_pivot;

      // if N_pivot is not in Parameter map either read in yaml file Rule if it exists or 
      // fall back to default value in CLASS
      if( Param.find("N_pivot") == Param.end())
      {
        result = runOptions->getValueOrDef<double>(60, "N_pivot");
      }
      // else "N_pivot" is contained in parameter map return the parameter value
      else
      {
        result = *Param["N_pivot"];
      }

    }
    
    /// capability to set k_pivot for consistent use within cosmobit 
    /// (to make sure it is consistent between CLASS and multimodecode)
    void set_k_pivot(double &result)
    {

      using namespace Pipes::set_k_pivot;

      result = runOptions->getValueOrDef<double>(0.05, "k_pivot");

    }


    void set_NuMasses_SM_baseline(map_str_dbl &result)
    {
      // Fallback for NuMasses_SM capability if Neutrino masses are not set, i.e. neither StandardModel_SLHA2 nor any child is in use
      // The Planck baseline analysis assumes a single massive neutrino with a fixed masss of 0.06 eV.
      using namespace Pipes::set_NuMasses_SM_baseline;

      result["mNu1"] = 0.06;
      result["mNu2"] = 0.0;
      result["mNu3"] = 0.0;

      result["N_ncdm"] = 1;

      result["mNu_tot_eV"] = 0.06;
    }

    void set_NuMasses_SM(map_str_dbl &result)
    {
      using namespace Pipes::set_NuMasses_SM;

      double mNu1, mNu2, mNu3;
      int N_ncdm = 0;

      // (PS) Heads up! The units in StandardModel_SLHA2 are GeV
      // Here we are using eV
      mNu1 = 1e9*(*Param["mNu1"]);
      mNu2 = 1e9*(*Param["mNu2"]);
      mNu3 = 1e9*(*Param["mNu3"]);

      if(mNu1 > 0.)
        N_ncdm++;
      if(mNu2 > 0.)
        N_ncdm++;
      if(mNu3 > 0.)
        N_ncdm++;

      result["mNu1"]=mNu1;
      result["mNu2"]=mNu2;
      result["mNu3"]=mNu3;

      result["N_ncdm"] = N_ncdm;

      result["mNu_tot_eV"] = mNu1 + mNu2 + mNu3;
    }

    void get_mNu_tot(double& result)
    {
      using namespace Pipes::get_mNu_tot;

      map_str_dbl NuMasses = *Dep::NuMasses_SM;
      result = NuMasses["mNu_tot_eV"];
    }

    void get_N_ur(double& result)
    {
      // Returns the effective number of ultrarelativistic species today
      using namespace Pipes::get_N_ur;

      map_str_dbl NuMassInfo = *Dep::NuMasses_SM;

      int N_ncdm = static_cast<int>(NuMassInfo["N_ncdm"]);
      double N_ur{};
      switch (N_ncdm)
      {
        case 1:
          N_ur = 2.0328;  // N_ur (today) = 2.0328 for 1 massive neutrino at CMB release
          break;
        case 2:
          N_ur= 1.0196;  // N_ur (today) = 1.0196 for 2 massive neutrino at CMB release
          break;
        case 3:
          N_ur = 0.00641;  // N_ur (today) = 0.00641 for 3 massive neutrinos at CMB release
          break;
        case 0:
          {
            std::ostringstream err;
            err << "Warning: All your neutrino masses are zero. The Planck baseline LCDM model assumes at least one massive neutrino.\n";
            err << "A cosmological model without massive neutrinos is not implemented in CosmoBit. If you want to consider this you can add it to the function 'get_N_ur' of the capability 'N_ur' ";
            CosmoBit_warning().raise(LOCAL_INFO, err.str());
          }
          break;
        default:
          {
            std::ostringstream err;
            err << "You are asking for more than three massive neutrino species.\n";
            err << "Such a case is not implemented in CosmoBit. If you want to consider this you can add it to the function 'get_N_ur' of the capability 'N_ur'.";
            CosmoBit_error().raise(LOCAL_INFO, err.str());
          }
      }

      result = N_ur;
      if (ModelInUse("etaBBN_rBBN_rCMB_dNurBBN_dNurCMB"))
      {
        // Check if the input for dNeff is negative (unphysical)
        static bool allow_negative_dNeff = runOptions->getValueOrDef<bool>(false,"allow_negative_delta_N_ur");
        const ModelParameters& NP_params = *Dep::etaBBN_rBBN_rCMB_dNurBBN_dNurCMB_parameters;
        double dNurCMB =  NP_params.at("dNur_CMB");
        double rCMB =  NP_params.at("r_CMB");
        if ( (!allow_negative_dNeff) && (dNurCMB < 0.0) )
        {
          std::string err = "A negative value for \"dNur_CMB\" is unphysical and is not allowed in CosmoBit by default!\n\n";
          err += "If you want to proceed with megative values, please set the \"allow_negative_delta_N_ur\"-rule to \"true\" within the yaml-file.";
          CosmoBit_error().raise(LOCAL_INFO,err.c_str());
        }

        // If the check is passed, set the result.
        result = pow(rCMB,4)*(N_ur) + dNurCMB;
      }
      logger() << "N_ur calculated to be " << result << EOM;
    }

    /// create a python dictionary with the inputs that have to be passed to class 
    /// setting parameters related to (massive) neutrinos & ncdm components
    void set_classy_NuMasses_Nur_input(pybind11::dict &result)
    {
      using namespace Pipes::set_classy_NuMasses_Nur_input;

      // make sure dict is empty
      result.clear();

      map_str_dbl NuMasses_SM = *Dep::NuMasses_SM;
      int N_ncdm = static_cast<int>(NuMasses_SM["N_ncdm"]);

      // set number of ultra relativistic species
      result["N_ur"] = *Dep::N_ur;

      // Number of non-cold DM species
      // if non-zero a mass & temperature for each species has to be 
      // passed to class
      if (N_ncdm > 0.)
      {
        result["N_ncdm"] = N_ncdm;   

        std::vector<double> m_ncdm = m_ncdm_classInput(NuMasses_SM);
        // head up: this explicitly assumed that all ncdm components have the same temperature!! todo ? 
        std::vector<double> T_ncdm(N_ncdm,*Dep::T_ncdm); 

        // create one string with m_ncdm masses and 
        // T_ncdm temperatures separated by commas (to match the input
        // format of class)
        std::ostringstream ss1, ss2;
        std::string separator;
        for (auto x : m_ncdm) 
        {
          ss1 << separator << x;
          separator = ",";
        }
        separator = "";
        for (auto x : T_ncdm) 
        {
          ss2 << separator << x;
          separator = ",";
        }
        
        result["m_ncdm"] = ss1.str();
        result["T_ncdm"] = ss2.str();
      }
      else
      {
        result["T_ncdm"] = *Dep::T_ncdm;
      }
    }

    /// create a python dictionary with the standard inputs that have to be passed
    /// to class: cosmological parameters (H0,omega_b,tau_reio,omega_cdm) & add 
    /// model dependent results for N_ur, neutrino masses & helium abundance
    /// here potential extra input options given in the yaml file are read in
    void set_classy_baseline_params(pybind11::dict &result)
    {
      using namespace Pipes::set_classy_baseline_params;

      //std::cout << " enter " << __PRETTY_FUNCTION__ << std::endl;

      // make sure dict is empty
      result.clear();

      // keep track if it is the first run -- if so 
      // some extra consistency checks to make sure no contradicting
      // values are in the classy python input dictionary
      static bool first_run = true;

      // Get the dictionary with inputs for the neutrino masses and merge it
      // into the empty results dictionary.
      pybind11::dict NuMasses_In = *Dep::classy_NuMasses_Nur_input;
      merge_pybind_dicts(result, NuMasses_In, first_run);

      // standard cosmological parameters (common to all CDM -like models)
      result["H0"] =            *Param["H0"];
      result["T_cmb"] =         *Param["T_cmb"];
      result["omega_b"] =       *Param["omega_b"];
      result["tau_reio"] =      *Param["tau_reio"];
      result["omega_cdm"] =     *Param["omega_cdm"];

      // set helium abundance (vector Helium_abundance.at(0): mean, .at(1): uncertainty)
      std::vector<double> Helium_abundance = *Dep::Helium_abundance;
      result["YHe"] = Helium_abundance.at(0); 

      // TODO: need to test if class or exo_class in use! does not work -> (JR) should be fixed with classy implementation
      // -> (JR again) not sure if that is actually true.. need to test.  
      if (ModelInUse("DecayingDM_general") || ModelInUse("AnnihilatingDM_general"))
      {
        // add Decaying/annihilating DM specific options to python dictionary passed to CLASS, consistency checks only executed in first run
        merge_pybind_dicts(result,*Dep::classy_parameters_EnergyInjection, first_run);
      }
      

      // Other Class input direct from the YAML file 
      // check if these are already contained in the input dictionary -- if so throw an error
      // only do it for the first run though
      static pybind11::dict yaml_input;
      if(first_run)
      {
        YAML::Node classy_dict;
        if (runOptions->hasKey("classy_dict"))
        {
          classy_dict = runOptions->getValue<YAML::Node>("classy_dict");
          for (auto it=classy_dict.begin(); it != classy_dict.end(); it++)
          {
            std::string name = it->first.as<std::string>();
            std::string value = it->second.as<std::string>();
            
            // Check if the key exists in the dictionary
            if (not result.contains(name.c_str()))
            {
              yaml_input[name.c_str()] = value;
            }
            // If it does, throw an error, there's some YAML conflict going on.
            else
            { 
              std::cout << pybind11::repr(result) << std::endl;
              CosmoBit_error().raise(LOCAL_INFO, "The key '" + name + "' already "
                "exists in the CLASSY dictionary. You are probably trying to override a model parameter. If you really"
                "want to do this you should define an extra function to set the class parameters for the model you "
                "are considering.");
            }
          }
        }
      }

      /// check that user did not try to pass k_pivot or N_star through class dictionary -- these are 
      /// fixed by capabilities to ensure consistent use throughout the code
      if (yaml_input.contains("k_pivot") || yaml_input.contains("N_star"))
      {

        CosmoBit_error().raise(LOCAL_INFO, "You tried to pass (one of) the parameters - k_pivot\n  - N_star\n"
              "to CLASS. If an inflationary model is in use these have to be consistent throughout the code."
              "Hence, the values have to be set via a capability. Add "
              "  - capability: k_pivot\n    function: set_k_pivot\n    options:\n"
              "      k_pivot: 0.02\n"
              "or\n"
              "  - capability: N_pivot\n    function: set_N_pivot\n    options:\n"
              "      N_pivot: 40\n"
              "to the Rules section of your yaml file. If you are not scanning an inflation model "
              "in combination with LCDM_no_primordial the rule for N_pivot will be ignored as it is an "
              "inflationary model parameter. It will therefore be set to the according value.");

      }
      // If the model LCDM_no_primordial is used the primordial power spectrum is computed independent of CLASS. Therefore
      // this option should not be passed by the user
      if (ModelInUse("LCDM_no_primordial"))
      {
        if (yaml_input.contains("P_k_ini type"))
        {
          CosmoBit_error().raise(LOCAL_INFO, "You are using the model 'LCDM_no_primordial' which means the (shape of the) primordial "
                "power spectrum is set by an additional inflation model. Therefore the computation of the pps happens "
                "independent of CLASS.\n"
                "GAMBIT will take care of setting all CLASS inputs regarding that consistently. Please remove the option "
                "'P_k_ini type' for the capability 'classy_baseline_params'");
        }
      }

      // add yaml options to python dictionary passed to CLASS, consistency checks only executed in first run
      merge_pybind_dicts(result,yaml_input, first_run);
      

      // At last: if the Planck likelihood is used add all relevant input parameters to the CLASS dictionary:
      // output: 'lCl, pCl, tCl', lensing: yes, non linear: halofit , l_max_scalar: 2508
      // Note: this has to be done *after* the yaml options are read it since it would trigger an error of duplicated 
      //   keys in case the user specifies one of the options themselves. The 'merge_pybind_dicts' routine will properly
      //   with concatennating the output values and choosing the maximum passed value for l_max_scalars. Contradictions
      //   in the lensing or non linear choice will rightfully trigger an error.
      if (ModelInUse("cosmo_nuisance_Planck_lite") || ModelInUse("cosmo_nuisance_Planck_TTTEEE")|| ModelInUse("cosmo_nuisance_Planck_TT"))
      {
        // add Planck like specific options to python dictionary passed to CLASS, consistency checks only executed in first run
        merge_pybind_dicts(result,*Dep::classy_PlanckLike_input, first_run);
      }

      first_run = false;
    }

    /// Set the LCDM_no_primordial_ps in classy. Depends on a parametrised primordial_ps.
    /// Looks at the parameters used in a run,
    /// and passes them to classy in the form of a Python dictionary.
    void set_classy_parameters_parametrised_ps(CosmoBit::Classy_input& result)
    {
      using namespace Pipes::set_classy_parameters_parametrised_ps;

      //std::cout << " enter " << __PRETTY_FUNCTION__ << std::endl;

      // Clean the input container
      result.clear();

      // Now need to pass the primordial power spectrum
      // TODO check whether A_s from MultiModeCode has the log taken or not.
      Parametrised_ps pps = *Dep::parametrised_power_spectrum;
      result.add_entry("n_s", pps.get_n_s());
      result.add_entry("ln10^{10}A_s", pps.get_A_s());

      // add k_pivot and N_pivot entries
      result.add_entry("P_k_ini type", "analytic_Pk");
      result.add_entry("k_pivot", *Dep::k_pivot);
      

      // if r = 0 only compute scalar modes, else tensor modes as well
      //
      // => don't explicitly set "modes" to 's' since it defaults to it. If you set it here anyways
      // you won't be able to run CLASS when only requesting background quantities (e.g. for BAO & SNe likelihoods)
      // as the perturbations module won't run and therefore the entry "modes" won't be read. 
      if(pps.get_r() == 0){} 
      else
      {
        // don't set to zero in CLASS dict as it won't be read if no tensor modes are requested
        result.add_entry("r", pps.get_r()); 
        result.add_entry("modes","t,s");
      }

      // Get standard cosmo parameters, nu masses, helium abundance &
      // extra run options for class passed in yaml file to capability
      // 'classy_baseline_params'
      // Note: this should only contain a value for 'k_pivot' if the model 'LCDM'
      //    is in use.
      pybind11::dict classy_base_dict = *Dep::classy_baseline_params;
      
      // Add classy_base_dict entries to the result dictionary of the type Classy_input
      std::string common_keys = result.add_dict(classy_base_dict);
      if(common_keys != "")
      {
        CosmoBit_error().raise(LOCAL_INFO, "The key(s) '" + common_keys + "' already "
                "exists in the CLASSY dictionary. You are probably trying to override a CLASS setting. Check that none "
                "of the parameters you pass through your yaml file through RunOptions for the capability 'classy_baseline_params' "
                "is in contradiction with any settings made via the dependency resolution by CosmoBit in the function '"+__func__+"'.");
      }
    }

    /// Set the LCDM_no_primordial_ps in classy. Depends on a primordial_ps.
    /// Looks at the parameters used in a run,
    /// and passes them to classy in the form of a Python dictionary.
    void set_classy_parameters_primordial_ps(CosmoBit::Classy_input& result)
    {
      using namespace Pipes::set_classy_parameters_primordial_ps;

      //std::cout << " enter " << __PRETTY_FUNCTION__ << std::endl;

      // Clean the input container
      result.clear();

      // Now need to pass the primordial power spectrum 
      static Primordial_ps pps{};
      pps = *Dep::primordial_power_spectrum;
      result.add_entry("modes", "t,s");

      result.add_entry("P_k_ini type", "pointer_to_Pk");
      result.add_entry("k_array", pps.get_k());
      result.add_entry("pks_array", pps.get_P_s());
      result.add_entry("pkt_array", pps.get_P_t());
      result.add_entry("lnk_size" , pps.get_vec_size()); // don't hard code but somehow make consistent with multimode @TODO -> test
      // pass pivot scale of external spectrum to CLASS
      result.add_entry("k_pivot", *Dep::k_pivot);
      result.add_entry("N_star",  *Dep::N_pivot);
      
      // Get standard cosmo parameters, nu masses, helium abundance &
      // extra run options for class passed in yaml file to capability
      // 'classy_baseline_params'
      // Note: this should only contain a value for 'k_pivot' if the model 'LCDM'
      //    is in use, which is not allowed in this function as the full pk
      //    is provided by an external model scanned in combination with 'LCDM_no_primordial'
      pybind11::dict classy_base_dict = *Dep::classy_baseline_params;

      // Add classy_base_dict entries to the result dictionary of the type Classy_input
      // The string common_keys will be empty if the two dictionaries 'result'
      // and 'classy_base_dict' have no keys in common. If so 'common_keys'
      // will be a concatenation of the duplicated entries.
      std::string common_keys = result.add_dict(classy_base_dict);

      // No check if something went wrong and some parameters were defined twice
      if(common_keys != "")
      {
        CosmoBit_error().raise(LOCAL_INFO, "The key(s) '" + common_keys + "' already "
                "exists in the CLASSY dictionary. You are probably trying to override a CLASS setting. Check that none "
                "of the parameters you pass through your yaml file through RunOptions for the capability 'classy_baseline_params' "
                "is in contradiction with any settings made via the dependency resolution by CosmoBit in the function '"+__func__+"'.");
      }
    }

    /// Function to set the generic inputs used for MultiModeCode.
    void set_multimode_inputs(Multimode_inputs &result)
    {
      using namespace Pipes::set_multimode_inputs;

      // Clear anything from previous run
      result = Multimode_inputs();

      // N_pivot is model parameter in all inflation models
      result.N_pivot = *Param["N_pivot"];
      
      // NOTE: if the flag 'ic_samp_flags' is set to 1 then phi0 is fixed to phi_ini0
      // and the fields velocities are assumed to be slow roll => no sampling over 
      // phi_init0 and dphi_init0. Set them as GAMBIT model parameters if you want to 
      // sample over these as well (and adopt the choice for the sampling flag s.t.
      // the entries for dphi_init0 are taken into account. You will also need to set
      // the prior ranges to be equal to the central value to prevent multimode from 
      // doing some sampling on its own.. )


      // Go through each allowed inflationary model and set the choice of potential parameters (vparams), 
      // the number of inflatons, and the initial conditions from the GAMBIT inputs.
      
      // N.B. Potential definitions are in the MultiModeCode file modpk_potential.f90, if you wish to 
      // study a new inflationary potential then implement it there.
      if (ModelInUse("Inflation_SR1quad"))
      {
        result.vparams.push_back(*Param["m2_inflaton"]);
        result.phi_init0.push_back(*Param["phi_init0"]);
        result.potential_choice = 1; // -> 0.5 m_i^2 phi_i^2 --- N-quadratic
        result.num_inflaton = 1;
        result.vparam_rows = 1;
      }
      else if (ModelInUse("Inflation_1natural"))
      {
        result.vparams.push_back(*Param["lambda"]); 
        result.vparams.push_back(*Param["faxion"]); // finv = 1.e0_dp/(10.e0_dp**faxion))
        result.phi_init0.push_back(*Param["phi_init0"]);
        result.potential_choice = 2; // -> lambda**4*(1.e0_dp+cos(finv*phi))
        result.num_inflaton = 1;
        result.vparam_rows = 2;
      }
      else if (ModelInUse("Inflation_1quar"))
      {
        result.vparams.push_back(*Param["lambda"]);
        result.phi_init0.push_back(*Param["phi_init0"]);
        result.potential_choice = 3; // -> 0.25 lambda_i phi_i^4 --- N-quartic
        result.num_inflaton = 1;
        result.vparam_rows = 1;
      }
      else if (ModelInUse("Inflation_1linear"))
      {
        result.vparams.push_back(*Param["lambda"]);
        result.phi_init0.push_back(*Param["phi_init0"]);
        result.potential_choice = 4; // -> lambda_i phi_i --- N-linear
        result.num_inflaton = 1;
        result.vparam_rows = 1;
      }
      else if (ModelInUse("Inflation_1mono23"))
      {
        result.vparams.push_back(*Param["lambda"]);
        result.phi_init0.push_back(*Param["phi_init0"]);
        result.potential_choice = 5; // -> 1.5 lambda_i phi_i^(2/3)
        result.num_inflaton = 1;
        result.vparam_rows = 1;
      }
      else if (ModelInUse("Inflation_1hilltopInf"))
      {
        // @Selim, I don't know the form of the potential so I can't figure 
        // out which one the corresponding one in multimode would be
        result.vparams.push_back(*Param["lambda"]); 
        result.vparams.push_back(*Param["mu"]);
        result.phi_init0.push_back(*Param["phi_init0"]);
        result.potential_choice = 5; 
        result.num_inflaton = 1;
        result.vparam_rows = 2;
      }

      // (JR) TODO: MultiMode segFaults if this is empty have to do this properly though 
      result.dphi_init0.push_back(1.);

      static bool first_run = true;
      // do some consistency check for inputs passed to multimode before the first run
      if(first_run)
      {
        int size = result.vparams.size();
        if (result.num_inflaton*result.vparam_rows != size)
        {
          std::ostringstream err;
          err << "Error in MultiModecode settings: the number of free parameters in a inflation model";
          err << " set through the vector 'vparams' has to match (num_inflaton) x (vparam_rows). In this case the ";
          err << "length is " << result.vparams.size() << " and the product is "  << result.num_inflaton*result.vparam_rows <<  " -- double check the model dependent settings in 'get_multimode_results'";
          CosmoBit_error().raise(LOCAL_INFO, err.str());
        }
        if (result.potential_choice == -1 || result.num_inflaton == -1 || result.vparam_rows == -1)
        {
          std::ostringstream err;
          err << "Error in MultiModecode settings: you did not set one (or more) of the parameters";
          err << "'potential_choice', 'num_inflaton' or 'vparam_rows' for your inflation model. ";
          err << "Double check the model dependent settings in 'get_multimode_results'";
          CosmoBit_error().raise(LOCAL_INFO, err.str());
        }
        first_run = false;
      }


      /// @TODO Separate into cases where we want the full ps and the parametrised ps

      result.k_min = runOptions->getValueOrDef<double>(1e-6,"k_min");
      result.k_max = runOptions->getValueOrDef<double>(1e+6,"k_max");
      result.numsteps = runOptions->getValueOrDef<int>(100, "numsteps");

      if (result.numsteps > 1000)
        CosmoBit_error().raise(LOCAL_INFO, "Currently MultiModeCode supports a maximum k-array size of 1000. Please change your yaml file settings.");


      //  Read in MultiMode settings from the yaml file
      // TODO use getValueOrDef as much as possible

      //-------------------------------------------------------------
      // Parameters to control how the ICs are sampled.
      //-------------------------------------------------------------
      result.ic_sampling = runOptions->getValue<int> ("ic_sampling");
      result.energy_scale = runOptions->getValue<double> ("energy_scale");
      result.numb_samples = runOptions->getValue<int> ("numb_samples");
      result.save_iso_N = runOptions->getValue<int> ("save_iso_N");
      result.N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this
      result.slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      result.instreheat = runOptions->getValue<int> ("instreheat");


      // (JR) @Selim: we should probably double check if we need
      // to set these parameters or if we can/should fix them as well. 
      //-------------------------------------------------------------
      // Control the output of analytic approximations for comparison.
      //-------------------------------------------------------------
      result.use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      result.evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      result.use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      result.get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");

      // has to be set to be the same as the scale entering class. todo
      // fix through some capability (probably the same setting the clac_pk_full variable for first version)
      result.k_pivot = *Dep::k_pivot; 
      // difference in k-space used when pivot-scale observables from mode equations are evaluated
      // Samples in uniform increments in log(k).
      result.dlnk = runOptions->getValue<double> ("dlnk");

    }

    /// Passes the inputs from the MultiModeCode initialisation function 
    /// and computes the outputs.
    /// TODO: split this up into primordial_ps and parametrised_ps versions.
    void get_multimode_primordial_ps(Primordial_ps &result)
    { 
      using namespace Pipes::get_multimode_primordial_ps;

      // Clear it all
      result = Primordial_ps();

      
      // Get the inflationary inputs
      Multimode_inputs inputs = *Dep::multimode_input_parameters;

      // TODO S.B shouldn't need these...
      // double N_pivot_prior_min = inputs.N_pivot;
      // double N_pivot_prior_max = inputs.N_pivot;

      // The parameters below are only used by multimode if the full Pk is requested. 
      int steps = inputs.numsteps;
      double kmin = inputs.k_min;
      double kmax = inputs.k_max;

      gambit_inflation_observables observables;

      //-------------------------------------------------------------
      // The function below calls the MultiModeCode backend
      //  for a given choice of inflationary model,
      //  which calculates the observables.
      //-------------------------------------------------------------

      try
      { 
        observables = BEreq::multimodecode_primordial_ps(inputs.num_inflaton,
                                                       inputs.potential_choice,
                                                       inputs.evaluate_modes,
                                                       inputs.get_runningofrunning,
                                                       byVal(&inputs.phi_init0[0]),
                                                       byVal(&inputs.dphi_init0[0]),
                                                       byVal(&inputs.vparams[0]),
                                                       inputs.N_pivot,
                                                       inputs.k_pivot,
                                                       inputs.dlnk,
                                                       steps,
                                                       kmin,
                                                       kmax,
                                                       inputs.vparam_rows);
      }
      catch(std::runtime_error &e)
      {
        logger() << e.what() << EOM;
        invalid_point().raise(e.what());
      }

      // If there's an error -> pass it to the helper function and invalidate the point
      if(observables.err != 0)
      {
        std::string message = multimode_error_handling(observables.err);
        logger() << message << EOM;
        invalid_point().raise(message);
      }

      // Fill up the GAMBIT prim. PS if we're good.
      result.fill_k(observables.k_array, inputs.numsteps);
      result.fill_P_s(observables.pks_array, inputs.numsteps);
      result.fill_P_s_iso(observables.pks_array_iso, inputs.numsteps);
      result.fill_P_t(observables.pkt_array, inputs.numsteps);

    }

    /// Passes the inputs from the MultiModeCode initialisation function 
    /// and computes the outputs.
    /// TODO: split this up into primordial_ps and parametrised_ps versions.
    void get_multimode_parametrised_ps(Parametrised_ps &result)
    { 
      using namespace Pipes::get_multimode_parametrised_ps;

      // Get the inflationary inputs
      Multimode_inputs inputs = *Dep::multimode_input_parameters;

      // TODO S.B shouldn't need these...
      // double N_pivot_prior_min = inputs.N_pivot;
      // double N_pivot_prior_max = inputs.N_pivot;

      gambit_inflation_observables observables;

      //-------------------------------------------------------------
      // The function below calls the MultiModeCode backend
      //  for a given choice of inflationary model,
      //  which calculates the observables.
      //-------------------------------------------------------------
      try
      {
        observables = BEreq::multimodecode_parametrised_ps(inputs.num_inflaton,
                                                         inputs.potential_choice,
                                                         inputs.evaluate_modes,
                                                         inputs.get_runningofrunning,
                                                         byVal(&inputs.phi_init0[0]),
                                                         byVal(&inputs.dphi_init0[0]),
                                                         byVal(&inputs.vparams[0]),
                                                         inputs.N_pivot,
                                                         inputs.k_pivot,
                                                         inputs.dlnk,
                                                         inputs.vparam_rows);
      }
      catch(std::runtime_error &e)
      {
        logger() << e.what() << EOM;
	invalid_point().raise(e.what());
      }
      // If there's an error -> pass it to the helper function and invalidate the point
      if(observables.err != 0)
      {
        std::string message = multimode_error_handling(observables.err);
        logger() << message << EOM;
        invalid_point().raise(message);
      }

      result.set_n_s(observables.ns);
      result.set_A_s(observables.As);
      result.set_r(observables.r);
 
    }

    void get_parametrised_ps_LCDM(Parametrised_ps &result)
    {
      using namespace Pipes::get_parametrised_ps_LCDM;

      // Check not using non-primordial version

      // (JR) got the error when the lines below were uncommented... todo check what's going on 
      //Problem with ModelInUse("LCDM_no_primordial").
      //This model is not known by CosmoBit::get_parametrised_ps_LCDM.  
      //Please make sure that it has been mentioned in some context in the
      //rollcall header declaration of this function.


      //if (ModelInUse("LCDM_no_primordial"))  
      //{
      //  CosmoBit_error().raise(LOCAL_INFO, "You cannot use the LCDM_no_primordial_ps model to get"
      //                      " a power spectrum!! Try the function get_multimode_parametrised_ps...");
      //}

      Parametrised_ps pps;
      pps.set_n_s(*Param["n_s"]);
      pps.set_A_s(*Param["ln10A_s"]);
      pps.set_r(0); 

      result = pps;
    }

    void print_parametrised_ps(map_str_dbl &result)
    {
      using namespace Pipes::print_parametrised_ps;

      Parametrised_ps pps  = *Dep::parametrised_power_spectrum;
      result = pps.get_parametrised_ps_map();
    }

    // still needs to be rewritten and dived into two functions. Also 
    // should fill the new 
    void get_parametrised_ps_SMASH(Parametrised_ps &result)
    {
      using namespace Pipes::get_parametrised_ps_SMASH;

      //---------------------------------------------------------------//
      // Below, SMASH inflationary potential is calculated and
      // the slow-roll apporximation is used to get the observables
      // A_s, n_s and r from model parameters and the inital value
      // of the inflaton potential at which infator starts its
      // slow-roll with zero-velocity.
      //---------------------------------------------------------------//

        std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

        vparams[0] = *Param["log10_xi"];
        vparams[1] = *Param["log10_beta"];
        vparams[2] = *Param["log10_lambda"];


        /* coefficients of:
         P(x) = -8+b*\phi**2+b*\xi*\phi**4+6*\xi**2*\phi**4
         */
        double smashp1[5] = { -8.0, 0, pow(10.0,vparams[1]), 0, (pow(10.0,vparams[1])*
                                                                 pow(10.0,vparams[0])+
                                                                 6.0*pow(10.0,vparams[0])*
                                                                 pow(10.0,vparams[0]))};

        double smashd1[8];

        gsl_poly_complex_workspace * wsmash = gsl_poly_complex_workspace_alloc (5);
        gsl_poly_complex_solve (smashp1, 5, wsmash, smashd1);
        gsl_poly_complex_workspace_free (wsmash);

        double badway[4] = {smashd1[0],smashd1[2],smashd1[4],smashd1[6]};
        double tempbw = 0;

        for(int i=0;i<4;i++)
        {
          if(badway[i]>tempbw) tempbw=badway[i];
        }

        vparams[3] = tempbw;

        //-------------------------------------------------------------
        //                  Sampling N_pivot
        //-------------------------------------------------------------
        double N_pivot = *Param["N_pivot"];

        // Solving the slow-roll equality to
        int status;
        int iter = 0, max_iter = 100000;
        const gsl_min_fminimizer_type *T;
        gsl_min_fminimizer *s;
        double m = tempbw*10.0; // Correct for.
        double a = 0.0, b = 100.0; // Correct for.
        gsl_function F;

        struct my_f_params params = { vparams[0], vparams[1], vparams[3] , N_pivot};

        F.function = &fn1;
        F.params = &params;

        T = gsl_min_fminimizer_brent;
        s = gsl_min_fminimizer_alloc (T);
        gsl_min_fminimizer_set (s, &F, m, a, b);

        do
        {
          iter++;
          status = gsl_min_fminimizer_iterate (s);

          m = gsl_min_fminimizer_x_minimum (s);
          a = gsl_min_fminimizer_x_lower (s);
          b = gsl_min_fminimizer_x_upper (s);

          status = gsl_min_test_interval (a, b, 0.01, 0.0);

        }
        while (status == GSL_CONTINUE && iter < max_iter);

        gsl_min_fminimizer_free (s);

        double smash_pot = 0.0;
        double smash_eps = 0.0;
        double smash_eta = 0.0;

        smash_pot = pot_SMASH(pow(10.0,*Param["log10_lambda"]),
                              a,
                              pow(10.0,*Param["log10_xi"]));
        smash_eps = SRparameters_epsilon_SMASH(a,
                                               pow(10.0,*Param["log10_beta"]),
                                               pow(10.0,*Param["log10_xi"]));
        smash_eta = SRparameters_eta_SMASH(a,
                                           pow(10.0,*Param["log10_beta"]),
                                           pow(10.0,*Param["log10_xi"]));

        double As_self = 0.0;
        double ns_self = 0.0;
        double r_self = 0.0;

        As_self = As_SR(smash_eps,smash_pot);
        ns_self = ns_SR(smash_eps,smash_eta);
        r_self  = r_SR(smash_eps);

        //-------------------------------------------------------------
        /* Have calculated the field value at N_pivot
         predicted by slow-roll approximation.
         Choosing a very close slightly higher
         value will be satisfactory for initial conditions.
         */
        //-------------------------------------------------------------
        Parametrised_ps pps;
        pps.set_n_s(ns_self);
        pps.set_A_s(As_self); // TODO check if we need to exponentiate
        pps.set_r(r_self); 

        result = pps;
    }


    /*
    /// have to think about how tho to that without causing a cricle in dep resolution
    //  -- need to pass the options below to class and 
    /// then class can return the vectors with the k, Pks and Pkt. todo
    void get_primordial_ps_SMASH(Primordial_ps &result)
    {
      using namespace Pipes::get_primordial_ps_SMASH;
      
      pybind11::dict classy_input;
      //-------------------------------------------------------------
      // Having CLASS primordial.h module calculate inflationary predictions of the SMASH potential.
      //-------------------------------------------------------------
      classy_input["P_k_ini type"] = "inflation_V";
      classy_input["potential"] = "smash_inflation";
      classy_input["V_0"] = *Param["log10_xi"];
      classy_input["V_1"] = *Param["log10_beta"];
      classy_input["V_3"] = *Param["log10_lambda"];
      classy_input["V_2"] = -100.; //Hard coding \chi inflaton potential boundaries
      
      // this does not do anything but and also does not work but I couldn't introduce a 
      // new compiler warning after Sanjay just removed all of them... 
      std::vector<double> k = result.get_k(); 
    }*/

    

// Getter functions for CL spectra from classy

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

    void class_get_unlensed_Cl_TE(std::vector<double>& result)
    {
      using namespace Pipes::class_get_unlensed_Cl_TE;
      result = BEreq::class_get_unlensed_cl("te");
    }

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

    void class_get_lensed_Cl_TE(std::vector<double>& result)
    {
      using namespace Pipes::class_get_lensed_Cl_TE;
      result = BEreq::class_get_lensed_cl("te");
    }

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

// ***************************************************************************************************************

    int diff_eq_rhs (double t, const double y[], double f[], void *params) // @Pat: this is the function for the differential equation that has to be solved in 'compute_dNeff_etaBBN_ALP' routine
    {
      /* RHS of differential equation
        dT/dt = 15/pi^2 (m_a n_a(t)/ tau_a) T^(-3) - H(T) T , where H(T) = 3.7978719e-7*T*T  // TODO: refer to Eq. number in paper when ready
        params: stores (m_a n_a(t)/ tau_a)
        y[0]: stores SM T[t0]
      */

      fast_interpolation injection_inter = *(static_cast<fast_interpolation*>(params));
      f[0] = (15.0/(4.0*pi*pi)) * injection_inter.interp(t)/pow(y[0], 3) - 3.7978719e-7*y[0]*y[0]*y[0];
      return GSL_SUCCESS;
    }

    void compute_dNeff_etaBBN_ALP(map_str_dbl &result)  // takes about 0.2 s with 3000 grid points and 1e-6 convergence criterion atm
    {
      using namespace Pipes::compute_dNeff_etaBBN_ALP;

      double dNeff, Neff_SM, Neff_ALP, eta_ratio;
      double temp_eta_ratio = 1; // temporary values to check convergence of iteration
      double temp_dNeff = 1;
      int ii=0;

      // --- precision parameters --
      double hstart = runOptions->getValueOrDef<double>(1e-6,"hstart"); // initial step size for GSL ODE solver
      double epsrel = runOptions->getValueOrDef<double>(1e-6,"epsrel"); // desired relative accuracy for GSL ODE solver and dNeff & etaBBN
      double max_iter = runOptions->getValueOrDef<int>(10,"max_iter"); // maximum number of iterations before error message is thrown if result not converged
      double grid_size = runOptions->getValueOrDef<int>(3000,"grid_size"); // number of (logarithmic) grid points in t

      double t0 = runOptions->getValueOrDef<double>(1e4,"t_initial"); // initial time in seconds
      double tf = runOptions->getValueOrDef<double>(1e12,"t_final"); // final time in seconds

      std::valarray<double> t_grid(grid_size), T_evo(grid_size), Tnu_grid(grid_size), na_grid(grid_size), injection_grid(grid_size); // arrays containing time dependence of variables

      SM_time_evo SM(t0,tf,grid_size);  // set time evolution of SM
      t_grid = SM.get_t_grid();         // has to be updated when solving the differential equation
      T_evo = SM.get_T_evo();

      // --- model parameters ----
      double Ya0 = *Dep::total_DM_abundance;
      double T0 = T_evo[0];
      double ssm_at_T0 = entropy_density_SM(T0, true);; // T0 in units of keV, set T_in_eV=True to interpret it correctly

      double na_t0 = Ya0 * ssm_at_T0;     // initial number density of a at t=t0, in units keV^3.
      double m_a = 1e-3*(*Param["ma0"]);  // mass of a in keV
      double tau_a = *Dep::lifetime;      // lifetime of a in seconds

      // loop over iterations until convergence or max_iter is reached
      while( ((fabs(temp_eta_ratio-eta_ratio) > epsrel) || fabs(temp_dNeff - dNeff) > epsrel) && (ii <= max_iter) )
      {
        temp_eta_ratio = eta_ratio;  // to check for convergence
        temp_dNeff = dNeff;

        SM.calc_H_int();
        na_grid = na_t0*exp(-3.0*SM.get_H_int())*exp(-t_grid/tau_a);    // na(t) in units keV^3
        injection_grid = (m_a/tau_a)*na_grid;                           // m_a*na(t)/tau_a in units keV^4/s
        Tnu_grid = SM.get_Tnu_evo()[0]*exp(-1.0*SM.get_H_int());        // T_nu(t) in units keV
        fast_interpolation injection_inter(t_grid, injection_grid);     // interpolating function

        gsl_odeiv2_system sys = {diff_eq_rhs, NULL, 1, &injection_inter};
        gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, hstart, 0.0, epsrel);

        double y[1] = {T_evo[0]}; // initial condition: T(t0) = T_SM(t0)
        double t=t0;

        for (int jj = 0; jj < grid_size; jj++)
        {
          double t_j = t_grid[jj];
          int status = gsl_odeiv2_driver_apply (d, &t, t_j, y);
          if (status != GSL_SUCCESS)
          {
            std::ostringstream err;
            err << "Failed to solve differential equation to compute dNeffCMB and etaBB for GeneralCosmoALP model at iteration step "<< ii <<". Invalidating point";
            invalid_point().raise(err.str());
          }
          T_evo[jj] = y[0];
        }

        gsl_odeiv2_driver_free (d);

        // update Hubble rate for next iteration
        SM.set_HT_evo(T_evo);

        eta_ratio = pow(T_evo[grid_size-1]/T_evo[0], 3) * exp(3.0*SM.get_H_int()[grid_size-1]);
        Neff_ALP = 3*pow(Tnu_grid[grid_size-1]/T_evo[grid_size-1], 4) * 3.85280407965; // (11/4)^(4/3) = 3.85280407965
        Neff_SM = 3*pow(Tnu_grid[0]/T_evo[0], 4) * 3.85280407965; // (11/4)^(4/3) = 3.85280407965
        dNeff = Neff_ALP - Neff_SM;
        ii ++;
      }

      // invalidate point if results not converged after 'max_iter' steps
      if( ((fabs(temp_eta_ratio-eta_ratio) > epsrel) || fabs(temp_dNeff - dNeff) > epsrel) && (ii >= max_iter) )
      {
        std::ostringstream err;
        err << "Computation of dNeffCMB and etaBBN for GeneralCosmoALP model did not converge after n = "<< ii <<" iterations. You can increase the maximum number of iterations with the run Option 'max_iter'. Invalidating point.";
        invalid_point().raise(err.str());
      }

      result["dNeff"] = dNeff;
      result["Neff_ratio"] = Neff_ALP/Neff_SM;
      result["eta_ratio"] = eta_ratio;
      logger() << "GeneralCosmoALP model: calculated Neff @BBN to be " << result["dNeff"] <<", and etaBB(ALP)/etaBBN(SM) = " << result["eta_ratio"] << ". Calculation converged after "<< ii <<" iterations." << EOM;

    }

/// -----------
/// Capabilities for general cosmological quantities
/// -----------

    void T_ncdm_SM(double &result)
    {
      using namespace Pipes::T_ncdm_SM;

      // set to 0.71611 in units of photon temperature, above the instantaneous decoupling value (4/11)^(1/3)
      // to recover Sum_i mNu_i/omega = 93.14 eV resulting from studies of active neutrino decoupling (hep-ph/0506164)
      result = 0.71611;
      // This standard values enters in many assumption entering class. Therefore changing this value in
      // the yaml file is disabled at the moment. If you still want to modify it uncomment the line below and
      // you can set is as a run option (as T_cmb).
      // result = runOptions->getValueOrDef<double>(0.71611,"T_ncdm");
    }

    void T_ncdm(double &result)
    {
      using namespace Pipes::T_ncdm;

      double rCMB = 1.0; // Default value if no energy injection is assumed.

      // If the "etaBBN_rBBN_rCMB_dNurBBN_dNurCMB" model is included in the scan,
      // we use rCMB of this model.
      if (ModelInUse("etaBBN_rBBN_rCMB_dNurBBN_dNurCMB"))
      {
        const ModelParameters& NP_params = *Dep::etaBBN_rBBN_rCMB_dNurBBN_dNurCMB_parameters;
        rCMB = NP_params.at("r_CMB");
      }

      // Take the SM value of T_ncdm (T_nu) and multiply it with the value of rCMB
      result = rCMB*(*Dep::T_ncdm_SM);
    }

    void calculate_eta0(double &result)
    {
      using namespace Pipes::calculate_eta0;

      double ngamma, nb;
      ngamma = 16*pi*zeta3*pow(*Param["T_cmb"]*_kB_eV_over_K_/_hc_eVcm_,3); // photon number density today
      nb = *Param["omega_b"]*3*100*1e3*100*1e3/_Mpc_SI_/_Mpc_SI_/(8*pi*_GN_cgs_* m_proton*1e9*eV2g); // baryon number density today

      result =  nb/ngamma;
      logger() << "Baryon to photon ratio (eta) today computed to be " << result << EOM;
    }

    void compute_Omega0_m(double &result)
    {
      using namespace Pipes::compute_Omega0_m;

      result =(*Dep::Omega0_b) + (*Dep::Omega0_cdm) + (*Dep::Omega0_ncdm);
    }

    void compute_Omega0_b(double &result)
    {
      using namespace Pipes::compute_Omega0_b;

      double h = *Param["H0"]/100;
      result =*Param["omega_b"]/h/h;
    }

    void compute_Omega0_cdm(double &result)
    {
      using namespace Pipes::compute_Omega0_cdm;

      double h = *Param["H0"]/100;
      result =*Param["omega_cdm"]/h/h;
    }

    void compute_Omega0_g(double &result)
    {
      using namespace Pipes::compute_Omega0_g;

      double h = *Param["H0"]/100.;
      result = (4.*_sigmaB_SI_/_c_SI_*pow(*Param["T_cmb"],4.)) / (3.*_c_SI_*_c_SI_*1.e10*h*h/_Mpc_SI_/_Mpc_SI_/8./pi/_GN_SI_);
    }

    void compute_n0_g(double &result)
    {
      using namespace Pipes::compute_n0_g;

      result = 2./pi/pi*zeta3 *pow(*Param["T_cmb"]*_kB_eV_over_K_,3.)/pow(_hP_eVs_*_c_SI_/2./pi,3)/100/100/100; // result per cm^3
    }

    void compute_Omega0_ur(double &result)
    {
      using namespace Pipes::compute_Omega0_ur;

      double N_ur = *Dep::N_ur;
      double Omega0_g = *Dep::Omega0_g;
      result = (N_ur)*7./8.*pow(4./11.,4./3.)* Omega0_g;
    }


    // (JR) delete when CLASS c interface is removed.
    void compute_Omega0_ncdm(double &result)
    {
      using namespace Pipes::compute_Omega0_ncdm;

      double mNu_tot_eV = *Dep::mNu_tot;
      double h = *Param["H0"]/100;

      result = mNu_tot_eV/(93.14*h*h);  // TODO: heads up: explicit assumption of T_ncdm = 0.71611 and T_cmb goes in here. Has to be generalised
    }

    void compute_Omega0_r(double &result)
    {
      using namespace Pipes::compute_Omega0_r;

      result = *Dep::Omega0_g + (*Dep::Omega0_ur);
      //std::cout << "omega0_g g " << *Dep::Omega0_g << " omega_ur " << *Dep::Omega0_ur <<  std::endl;
    }

/*
    (PS: Not sure if we really need this step in between. For all our casess so far we directly could map from eta0 to eta_BBN (up to constnat factors in front))
 * 
    void calculate_etaCMB_SM(double &result)
    {
      using namespace Pipes::calculate_etaCMB_SM;

      result =  *Dep::eta0; // in SM the baryon to photon ratio does not change between today an CMB release
      logger() << "Baryon to photon ratio (eta) @CMB computed to be " << result << EOM;
    }
*/


/// -----------
/// BBN related functions (call AlterBBN, BBN abundances & Likelihood)
/// -----------

    /// add parameters of relicparam structure that should be set to non-default values
    /// to the AlterBBN_input map.
    /// If you want to modify a parameter which has not been used in CosmoBit before simply
    /// add it to the function 'fill_cosmomodel' in 'AlterBBN_<version>.cpp' and to the
    /// set of 'known' parameters 'known_relicparam_options'
    void AlterBBN_Input(map_str_dbl &result)
    {
      using namespace Pipes::AlterBBN_Input;


      // If we are using some of the "non-standard energy content" models, set the
      // inputs for the AlterBBN_input map according to the parameters of that model.
      // In case we are not using one of these models, we use the default values
      // (i.e. eta inferred from LCDM, Nnu = 3.046 and dNnu = 0)
      if (ModelInUse("etaBBN_rBBN_rCMB_dNurBBN_dNurCMB"))
      {
        const ModelParameters& NP_params = *Dep::etaBBN_rBBN_rCMB_dNurBBN_dNurCMB_parameters;

        double dNurBBN =  NP_params.at("dNur_BBN");

        // Check if the input for dNeff is negative (unphysical)
        static bool allow_negative_dNeff = runOptions->getValueOrDef<bool>(false,"allow_negative_delta_N_ur");
        if ( (!allow_negative_dNeff) && (dNurBBN < 0.0) )
        {
          std::string err = "A negative value for \"dNur_BBN\" is unphysical and is not allowed in CosmoBit by default!\n\n";
          err += "If you want to proceed with megative values, please set the \"allow_negative_delta_N_ur\"-rule to \"true\" within the yaml-file.";
          CosmoBit_error().raise(LOCAL_INFO,err.c_str());
        }

        //If check is passed, set inputs.
        result["eta0"] = NP_params.at("eta_BBN");    // eta AFTER BBN (variable during)
        result["Nnu"]=3.046*pow(NP_params.at("r_BBN"),4); // contribution from SM neutrinos
        result["dNnu"]=dNurBBN;    // dNnu: within AlterBBN scenarios in which the sum Nnu+dNnu is the same are identical
      }
      else
      {
        result["eta0"] = *Dep::eta0;    // eta AFTER BBN equal to eta today calculated from omgea_b
        result["Nnu"]=3.046; // contribution from SM neutrinos
        result["dNnu"]=0.;    // no extra ur species in standard LCDM model
      }

      result["failsafe"] = runOptions->getValueOrDef<double>(3,"failsafe");
      result["err"] = runOptions->getValueOrDef<double>(3,"err");

      logger() << "Set AlterBBN with parameters eta = " << result["eta0"] << ", Nnu = " << result["Nnu"] << ", dNnu = " << result["dNnu"] << EOM;
      logger() << "     and error params: failsafe = " << result["failsafe"] << ", err = " << result["err"] << EOM;
    }


    /// (JR) this function is huge (200 lines..) -- considering that it is acually only there to call 
    /// a BE function from AlterBBN, can we split that up somehow or define utility functions
    /// to do some of the checks and error calculation things so it is not totally blown up? TODO
    void compute_BBN_abundances(CosmoBit::BBN_container &result)
    {
      using namespace Pipes::compute_BBN_abundances;

      // global variable of AlterBBN (# computed element abundances)
      const int NNUC = BEreq::get_NNUC();  
      result.set_abund_map(BEreq::get_abund_map_AlterBBN());
      // init arrays in BBN_container with right length
      result.init_arr_size(NNUC);         

      // in AlterBBN ratioH and cov_ratioH are arrays of fixed length
      // with certain compiler versions (gcc 5.4.0) we have seen memory corruption problems
      // if we initialise these as 
      // double ratioH[NNUC+1]
      // since the memory was not allocated properly. Fixed size arrays do not seem to be
      // properly supported even though there are no errors at compile time. 
      // using a unique pointer for ratioH and a 2d vector for cov_ratioH avoids 
      // these problems. 
      auto deleter = [&](double* ptr){delete [] ptr;};
      std::unique_ptr<double[], decltype(deleter)> ratioH(new double[NNUC+1](), deleter);
      std::unique_ptr<double[], decltype(deleter)> cov_ratioH(new double[(NNUC+1)*(NNUC+1)](), deleter);

      static bool first = true;
      const bool use_fudged_correlations = (runOptions->hasKey("correlation_matrix") && runOptions->hasKey("elements"));
      const bool use_relative_errors = runOptions->hasKey("relative_errors");
      static std::vector<double> relerr(NNUC+1, -1.0);
      static std::vector<std::vector<double>> corr(NNUC+1, std::vector<double>(NNUC+1, 0.0));
      if (use_fudged_correlations && first)
      {
        for (int ie = 1; ie < NNUC; ie++) corr.at(ie).at(ie) = 1.;
        std::vector<str> elements =  runOptions->getValue< std::vector<str> >("elements");
        std::vector<std::vector<double>> tmp_corr = runOptions->getValue< std::vector<std::vector<double>> >("correlation_matrix");
        std::vector<double> tmp_relerr;

        unsigned int nelements = elements.size();

        // Check if the size of the list of elements and the size of the correlation matrix agree
        if (nelements != tmp_corr.size())
        {
          std::ostringstream err;
          err << "The length of the list \'elements\' and the size of the correlation matrix \'correlation_matrix\' do not agree";
          CosmoBit_error().raise(LOCAL_INFO, err.str());
        }

        // If the relative errors are also given, then do also a check if the length of the list is correct and if the entries are positive.
        if (use_relative_errors)
        {
          tmp_relerr = runOptions->getValue< std::vector<double> >("relative_errors");
          if (nelements != tmp_relerr.size())
          {
            std::ostringstream err;
            err << "The length of the list \'elements\' and the length of \'relative_errors\' do not agree";
            CosmoBit_error().raise(LOCAL_INFO, err.str());
          }
          for (std::vector<double>::iterator it = tmp_relerr.begin(); it != tmp_relerr.end(); it++)
          {
            if (*it <= 0.0)
            {
              std::ostringstream err;
              err << "One entry for the relative error is not positive";
              CosmoBit_error().raise(LOCAL_INFO, err.str());
            }
          }
        }

        // Check if the correlation matrix is symmetric
        for (std::vector<std::vector<double>>::iterator it = tmp_corr.begin(); it != tmp_corr.end(); it++)
        {
          if (it->size() != nelements)
          {
            std::ostringstream err;
            err << "The correlation matrix is not a square matrix";
            CosmoBit_error().raise(LOCAL_INFO, err.str());
          }
        }

        // Check if the entries in the correlation matrix are reasonable
        for (unsigned int ie=0; ie<nelements; ie++)
        {
          //Check if the diagonal entries are equal to 1.
          if (std::abs(tmp_corr.at(ie).at(ie) - 1.) > 1e-6)
          {
            std::ostringstream err;
            err << "Not all diagonal elements of the correlation matirx are 1.";
            CosmoBit_error().raise(LOCAL_INFO, err.str());
          }
          for (unsigned int je=0; je<=ie; je++)
          {
            //Check for symmetry
            if (std::abs(tmp_corr.at(ie).at(je) - tmp_corr.at(je).at(ie)) > 1e-6)
            {
              std::ostringstream err;
              err << "The correlation matrix is not symmetric";
              CosmoBit_error().raise(LOCAL_INFO, err.str());
            }
            // Check if the off-diagonal elements are between -1 and 1.
            if (std::abs(tmp_corr.at(ie).at(je)) >= 1. && (ie != je))
            {
              std::ostringstream err;
              err << "The off-diagonal elements of the correlation matrix are not sensible (abs(..) > 1)";
              CosmoBit_error().raise(LOCAL_INFO, err.str());
            }
          }
        }

        // Check if the elements which are passed via runOptions are actually known.
        std::map<std::string, int> abund_map = result.get_abund_map();
        for (std::vector<str>::iterator it = elements.begin(); it != elements.end(); it++)
        {
          if (abund_map.count(*it) == 0)
          {
            std::ostringstream err;
            err << "I do not recognise the element \'" << *it << "\'";
            CosmoBit_error().raise(LOCAL_INFO, err.str());
          }
        }

        // Populate the correlation matrix
        for (std::vector<str>::iterator it1 = elements.begin(); it1 != elements.end(); it1++)
        {
          int ie  =  abund_map[*it1];
          int i = std::distance( elements.begin(), it1 );
          // If the relative errors are given fill it with the respective values (-1.0 refers to no errors given).
          if (use_relative_errors) relerr.at(ie) = tmp_relerr.at(i);
          for (std::vector<str>::iterator it2 = elements.begin(); it2 != elements.end(); it2++)
          {
            int je  =  abund_map[*it2];
            int j = std::distance( elements.begin(), it2 );
            corr.at(ie).at(je) = tmp_corr.at(i).at(j);
          }
        }

        first = false;
      }

      // fill AlterBBN_input map with the parameters for the model in consideration
      map_str_dbl AlterBBN_input = *Dep::AlterBBN_Input; 

      // call AlterBBN routine to calculate element abundances (& errors -- depending 
      // on error calculation settings made with parameters 'err' and failsafe set in 
      // 'AlterBBN_Input')
      int nucl_err = BEreq::call_nucl_err(AlterBBN_input, &ratioH[0], &cov_ratioH[0]);

      // TODO: replace .at() by [] to speed up
      std::vector<double> err_ratio(NNUC+1,0);
      if (use_fudged_correlations)
      {
        for (int ie=1;ie<=NNUC;ie++)
        {
          if (use_relative_errors && (relerr.at(ie) > 0.0))
          {
            err_ratio.at(ie) =  relerr.at(ie) * ratioH[ie];
          }
          else
          {
            // get every diagonal element (row and line 0not filled)
            err_ratio.at(ie) = sqrt(cov_ratioH[ie*(NNUC+1)+ie]);
          }
        }
      }

      // fill abundances and covariance matrix of BBN_container with results from AlterBBN
      if (nucl_err)
      {
        for(int ie=1;ie<=NNUC;ie++)
        {
          result.set_BBN_abund(ie, ratioH[ie]);
          for(int je=1;je<=NNUC;je++)
          {
            if (use_fudged_correlations)
              result.set_BBN_covmat(ie, je, corr.at(ie).at(je) * err_ratio.at(ie) * err_ratio.at(je));
            else
              result.set_BBN_covmat(ie, je, cov_ratioH[ie*(NNUC+1)+je]);
          }
        }
      }
      else
      {
        std::ostringstream err;
        err << "AlterBBN calculation for primordial element abundances failed. Invalidating Point.";
        invalid_point().raise(err.str());
      }
    }

    void get_Helium_abundance(std::vector<double> &result)
    {
      using namespace Pipes::get_Helium_abundance;

      CosmoBit::BBN_container BBN_res = *Dep::BBN_abundances;
      std::map<std::string, int> abund_map = BBN_res.get_abund_map();

      result.clear();
      result.push_back(BBN_res.get_BBN_abund(abund_map["Yp"]));
      result.push_back(sqrt(BBN_res.get_BBN_covmat(abund_map["Yp"], abund_map["Yp"])));

      //std::cout << "Helium Abundance from AlterBBN: " << result.at(0) << " +/- " << result.at(1) << std::endl;
      logger() << "Helium Abundance from AlterBBN: " << result.at(0) << " +/- " << result.at(1) << EOM;
    }

    void get_Helium3_abundance(std::vector<double> &result)
    {
      using namespace Pipes::get_Helium3_abundance;

      CosmoBit::BBN_container BBN_res = *Dep::BBN_abundances;
      std::map<std::string, int> abund_map = BBN_res.get_abund_map();

      result.clear();
      result.push_back(BBN_res.get_BBN_abund(abund_map["He3"]));
      result.push_back(sqrt(BBN_res.get_BBN_covmat(abund_map["He3"], abund_map["He3"])));

      logger() << "Helium 3 Abundance from AlterBBN: " << result.at(0) << " +/- " << result.at(1) << EOM;
    }

    void get_Deuterium_abundance(std::vector<double> &result)
    {
      using namespace Pipes::get_Deuterium_abundance;

      CosmoBit::BBN_container BBN_res = *Dep::BBN_abundances;
      std::map<std::string, int> abund_map = BBN_res.get_abund_map();

      result.clear();
      result.push_back(BBN_res.get_BBN_abund(abund_map["D"]));
      result.push_back(sqrt(BBN_res.get_BBN_covmat(abund_map["D"], abund_map["D"])));

      logger() << "Deuterium Abundance from AlterBBN: " << result.at(0) << " +/- " << result.at(1) << EOM;
    }

    void get_Lithium7_abundance(std::vector<double> &result)
    {
      using namespace Pipes::get_Lithium7_abundance;

      CosmoBit::BBN_container BBN_res = *Dep::BBN_abundances;
      std::map<std::string, int> abund_map = BBN_res.get_abund_map();

      result.clear();
      result.push_back(BBN_res.get_BBN_abund(abund_map["Li7"]));
      result.push_back(sqrt(BBN_res.get_BBN_covmat(abund_map["Li7"], abund_map["Li7"])));

      logger() << "Lithium Abundance from AlterBBN: " << result.at(0) << " +/- " << result.at(1) << EOM;
    }

    void get_Beryllium7_abundance(std::vector<double> &result)
    {
      using namespace Pipes::get_Beryllium7_abundance;

      CosmoBit::BBN_container BBN_res = *Dep::BBN_abundances;
      std::map<std::string, int> abund_map = BBN_res.get_abund_map();

      result.clear();
      result.push_back(BBN_res.get_BBN_abund(abund_map["Be7"]));
      result.push_back(sqrt(BBN_res.get_BBN_covmat(abund_map["Be7"], abund_map["Be7"])));

      logger() << "Beryllium Abundance from AlterBBN: " << result.at(0) << " +/- " << result.at(1) << EOM;
    }

    void compute_BBN_LogLike(double &result)
    {
      using namespace Pipes::compute_BBN_LogLike;

      double chi2 = 0;
      int ii = 0;
      int ie,je,s;

      CosmoBit::BBN_container BBN_res = *Dep::BBN_abundances; // fill BBN_container with abundance results from AlterBBN
      std::map<std::string, int> abund_map = BBN_res.get_abund_map();

      const str filename = runOptions->getValue<std::string>("DataFile"); // read in BBN data file
      const str path_to_file = GAMBIT_DIR "/CosmoBit/data/BBN/" + runOptions->getValue<std::string>("DataFile");
      static ASCIIdictReader data(path_to_file);
      static bool read_data = false;
      static std::map<std::string,std::vector<double>> dict;
      static int nobs;

      if(read_data == false)
      {
        logger() << "BBN data read from file '"<<filename<<"'." << EOM;
        nobs = data.nrow();
        dict = data.get_dict();
        if(data.duplicated_keys()==true) // check for double key entries in data file
        {
          std::ostringstream err;
          err << "Double entry for one element in BBN data file '"<< filename<<"'. \nYou can only enter one measurement per element.";
          CosmoBit_error().raise(LOCAL_INFO, err.str());
        }
        if((dict.count("Yp")>0 and dict.count("He4")>0) or (dict.count("D")>0 and dict.count("H2")>0))
        {
          std::ostringstream err;
          err << "Double entry for 'Yp' and 'He4' or 'D' and 'H2' in BBN data file '"<< filename<<"'. \nYou can only enter one measurement per element ('Yp' OR 'He4', 'D' OR 'H2').";
          CosmoBit_error().raise(LOCAL_INFO, err.str());
        }
        read_data = true;
      }

      // init vectors with observations, predictions and covmat
      double prediction[nobs],observed[nobs],sigmaobs[nobs],translate[nobs];
      gsl_matrix *cov = gsl_matrix_alloc(nobs, nobs);
      gsl_matrix *invcov = gsl_matrix_alloc(nobs, nobs);
      gsl_permutation *p = gsl_permutation_alloc(nobs);

      // iterate through observation dictionary to fill observed, sigmaobs and prediction arrays
      for(std::map<std::string,std::vector<double>>::iterator iter = dict.begin(); iter != dict.end(); ++iter)
      {
        std::string key = iter->first; // iter = ["element key", [mean, sigma]]
        if(abund_map.count(key)!=1)   // throw error if element not contained in abundance map was entered in data file
        {
          std::ostringstream err;
          err << "Unknown element '"<< key <<"' in BBN data file '"<< filename<<"'. \nYou can only enter 'Yp' or 'He4', 'D' or 'H2', 'He3', 'Li7'.";
          CosmoBit_error().raise(LOCAL_INFO, err.str());
        }
        translate[ii]=abund_map[key]; // to order observed and predicted abundances consistently
        observed[ii]=iter->second[0];
        sigmaobs[ii]=iter->second[1];
        prediction[ii]= BBN_res.get_BBN_abund(abund_map[key]);
        ii++;
      }

      // fill covmat
      for(ie=0;ie<nobs;ie++) {for(je=0;je<nobs;je++) gsl_matrix_set(cov, ie, je,pow(sigmaobs[ie],2.)*(ie==je)+BBN_res.get_BBN_covmat(translate[ie], translate[je]));}

      // Compute the LU decomposition and inverse of cov mat
      gsl_linalg_LU_decomp(cov,p,&s);
      gsl_linalg_LU_invert(cov,p,invcov);

      // Compute the determinant of the inverse of the covmat
      double det_cov = gsl_linalg_LU_det(cov,s);

      // compute chi2
      for(ie=0;ie<nobs;ie++) for(je=0;je<nobs;je++) chi2+=(prediction[ie]-observed[ie])*gsl_matrix_get(invcov,ie,je)*(prediction[je]-observed[je]);
      //std::cout << "    BBN Like: chi2 = " << chi2 << "  factor " <<  log(pow(2*pi,nobs)*det_cov) << "  det cov = " << det_cov << std::endl; 
      result = -0.5*(chi2 + log(pow(2*pi,nobs)*det_cov));

      //std::cout << "    BBN LogLike computed to be: " << result << std::endl;
      logger() << "BBN LogLike computed to be: " << result << EOM;

      gsl_matrix_free(cov);
      gsl_matrix_free(invcov);
      gsl_permutation_free(p);
    }

/// -----------
/// Background Likelihoods (SNe + H0 + BAO + Sigma8)
/// -----------

    void compute_H0_LogLike(double &result)
    {
      using namespace Pipes::compute_H0_LogLike;

      const str filename = runOptions->getValue<std::string>("DataFile");
      const str path_to_file = GAMBIT_DIR "/CosmoBit/data/H0/" + runOptions->getValue<std::string>("DataFile");
      static ASCIItableReader data(path_to_file);
      static bool read_data = false;

      if(read_data == false)
      {
        logger() << "H0 data read from file '"<<filename<<"'." << EOM;
        std::vector<std::string> colnames = initVector<std::string>("mean", "sigma");
        data.setcolnames(colnames);

        if(data.getnrow() != 1)
        {
          std::ostringstream err;
          err << data.getnrow() << " data points for H0 measurement data read from '"<< filename << "'.\n Only one expected.";
          CosmoBit_error().raise(LOCAL_INFO, err.str());
        }
        read_data = true;
      }
      result = -0.5 * pow(*Param["H0"] - data["mean"][0],2)/ pow(data["sigma"][0],2);
    }

    void compute_sigma8_LogLike(double &result)
    {
      using namespace Pipes::compute_sigma8_LogLike;

      double sigma8,Omega0_m; // cosmo distances and sound horizon at drag,Hubble, theoretical prediction
      double chi2 =0;

      const str filename = runOptions->getValue<std::string>("DataFile");
      const str path_to_file = GAMBIT_DIR "/CosmoBit/data/sigma8/" + runOptions->getValue<std::string>("DataFile");
      static ASCIItableReader data(path_to_file);
      static bool read_data = false;
      static int nrow;

      sigma8 = BEreq::class_get_sigma8();
      Omega0_m = *Dep::Omega0_m;

      if(read_data == false)
      {
        logger() << "Sigma8 data read from file '"<<filename<<"'." << EOM;
        std::vector<std::string> colnames = initVector<std::string>("Omega_m_ref","Omega_m_index", "bestfit","sigma");
        data.setcolnames(colnames);
        nrow=data.getnrow();

        read_data = true;
      }

      for(int ie=0;ie<nrow;ie++)
      {
        chi2 += pow(( sigma8*pow(Omega0_m/data["Omega_m_ref"][ie], data["Omega_m_index"][ie]) - data["bestfit"][ie])/ data["sigma"][ie],2);
      }
      result = -0.5*chi2;
    }

/***********************************************/
/* Interface to CLASS Python wrapper, 'classy' */
/***********************************************/

    /// Initialises the container within CosmoBit from classy. This holds
    /// an instance of the classy class Class() (Yep, I know...)
    /// which can be handed over to MontePython, or just used to compute
    /// some observables.
    void set_classy_input(CosmoBit::Classy_input& result)
    {
      using namespace Pipes::set_classy_input;

      result = *Dep::classy_primordial_parameters;
    }

    /// Initialises the container within CosmoBit from classy, but designed specifically
    /// to be used when MontePython is in use. This will ensure additional outputs are
    /// computed by classy CLASS to be passed to MontePython.
    void set_classy_input_with_MPLike(CosmoBit::Classy_input& result)
    {
      using namespace Pipes::set_classy_input_with_MPLike;

      // get extra cosmo_arguments from MP (gives a dictionary with output values that need
      // to be set for the class run)
      static pybind11::dict MP_cosmo_arguments = *Dep::cosmo_args_from_MPLike;
      logger() << LogTags::debug << "Extra cosmo_arguments needed from MP Likelihoods: ";
      logger() << pybind11::repr(MP_cosmo_arguments) << EOM;

      result = *Dep::classy_primordial_parameters;

      // add the arguments from Mp_cosmo_arguments which are not yet in cosmo_input_dict to it
      // also takes care of merging the "output" key values
      result.merge_input_dicts(MP_cosmo_arguments); // TODO (JR) use merge_pybind_dict routine instead -> don't need to duplicate code
      //merge_pybind_dicts(result, MP_cosmo_arguments);
    }


  void set_classy_parameters_EnergyInjection_AnnihilatingDM(pybind11::dict &result)
  {
    using namespace Pipes::set_classy_parameters_EnergyInjection_AnnihilatingDM;

    // make sure nothing from previous run is contained
    result.clear();

    // Set relevant inputs for the scenario of s-wave annihilating DM
    const ModelParameters& NP_params = *Dep::AnnihilatingDM_general_parameters;
    result["DM_annihilation_cross_section"] = NP_params.at("sigmav");
    result["DM_annihilation_mass"] = NP_params.at("mass");

    // Get the results from the DarkAges tables that hold extra information to be passed to the CLASS thermodynamics structure
    static DarkAges::fz_table fz;
    fz = *Dep::energy_injection_efficiency;
    bool f_eff_mode = fz.f_eff_mode;

    // flag passed to CLASS to signal that the energy_deposition_function is coming from GAMBIT
    // we patched exoclass to accept this. An alternative way without patching would be to write the tables to disk &
    // just have CLASS read in the file. To avoid the repeated file writing & deleting we pass pointers to the vector/arrays
    // to CLASS instead
    if (f_eff_mode)
    {
      result["f_eff_type"] = "pointer_to_fz_eff";
    }
    else
    {
      result["f_eff_type"] = "pointer_to_fz_channel";
    }

    // set the lengths of the input tables (since we are passing pointers to arrays CLASS has to know how long they are)
    result["energyinj_coef_num_lines"] = fz.redshift.size();

    // add the pointers to arrays class needs to know about to input dictionary
    // Note:
    //    - memory addresses are passed as strings (python wrapper for CLASS converts every entry to a string internally so
    //      we need to do that for the memory addresses before python casts them to something else)
    result["energyinj_coef_z"] = memaddress_to_uint(fz.redshift.data());
    if (f_eff_mode)
    {
      result["energyinj_coef_tot"] = memaddress_to_uint(fz.f_eff.data());
    }
    else
    {
      result["energyinj_coef_heat"] = memaddress_to_uint(fz.f_heat.data());
      result["energyinj_coef_lya"] = memaddress_to_uint(fz.f_lya.data());
      result["energyinj_coef_ionH"] = memaddress_to_uint(fz.f_hion.data());
      result["energyinj_coef_ionHe"] = memaddress_to_uint(fz.f_heion.data());
      result["energyinj_coef_lowE"] = memaddress_to_uint(fz.f_lowe.data());
    }
  }

  void set_classy_parameters_EnergyInjection_DecayingDM(pybind11::dict &result)
  {
    using namespace Pipes::set_classy_parameters_EnergyInjection_DecayingDM;

    // make sure nothing from previous run is contained
    result.clear();

    // Set relevant inputs for the scenario of decaying DM
    const ModelParameters& NP_params = *Dep::DecayingDM_general_parameters;
    result["DM_decay_tau"] = NP_params.at("lifetime");
    result["DM_decay_fraction"] = NP_params.at("fraction");

    // Get the results from the DarkAges tables that hold extra information to be passed to the CLASS thermodynamics structure
    static DarkAges::fz_table fz;
    fz = *Dep::energy_injection_efficiency;
    bool f_eff_mode = fz.f_eff_mode;

    // flag passed to CLASS to signal that the energy_deposition_function is coming from GAMBIT
    // we patched exoclass to accept this. An alternative way without patching would be to write the tables to disk &
    // just have CLASS read in the file. To avoid the repeated file writing & deleting we pass pointers to the vector/arrays
    // to CLASS instead
    if (f_eff_mode)
    {
      result["f_eff_type"] = "pointer_to_fz_eff";
    }
    else
    {
      result["f_eff_type"] = "pointer_to_fz_channel";
    }

    // set the lengths of the input tables (since we are passing pointers to arrays CLASS has to know how long they are)
    result["energyinj_coef_num_lines"] = fz.redshift.size();

    // add the pointers to arrays class needs to know about to input dictionary
    // Note:
    //    - memory addresses are passed as strings (python wrapper for CLASS converts every entry to a string internally so
    //      we need to do that for the memory addresses before python casts them to something else)
    result["energyinj_coef_z"] = memaddress_to_uint(fz.redshift.data());
    if (f_eff_mode)
    {
      result["energyinj_coef_tot"] = memaddress_to_uint(fz.f_eff.data());
    }
    else
    {
      result["energyinj_coef_heat"] = memaddress_to_uint(fz.f_heat.data());
      result["energyinj_coef_lya"] = memaddress_to_uint(fz.f_lya.data());
      result["energyinj_coef_ionH"] = memaddress_to_uint(fz.f_hion.data());
      result["energyinj_coef_ionHe"] = memaddress_to_uint(fz.f_heion.data());
      result["energyinj_coef_lowE"] = memaddress_to_uint(fz.f_lowe.data());
    }
  }

  /// add all inputs for CLASS needed to produce the correct output to be 
  /// able to compute the Planck CMB likelihoods
  void set_classy_PlanckLike_input(pybind11::dict &result)
  {
    using namespace Pipes::set_classy_PlanckLike_input;

    // make sure nothing from previous run is contained
    result.clear();

    static std::ostringstream output;
    static std::ostringstream l_max_scalars;

    static bool first = true;
    if(first)
    {
      int lmax = -1;
      bool needs_tCl = false;
      bool needs_pCl = false;

      // Get requirements of the loaded likelihoods in the plc backend
      BEreq::plc_required_Cl(lmax,needs_tCl,needs_pCl);

      // Prepare the classy input for "output"
      // -- The likelihoods need the lensed Cl such that lensing is required everytime
      output << "lCl";
      // -- Are additional Cl, other to Cl_phiphi, required?
      if (needs_tCl)
        output << ", tCl";
      if (needs_pCl)
        output << ", pCl";

      // Prepare the classy input for "l_max_scalars"
      l_max_scalars << lmax;

      first = false;
    }

    result["lensing"] = "yes";
    result["non linear"] = "halofit";
    result["output"] = output.str();
    result["l_max_scalars"] = l_max_scalars.str();
  }



    /* Classy getter functions */

    /// Energy densities *today* (Omega0)

    /// Matter
    void get_Omega0_m_classy(double& result)
    {
      using namespace Pipes::get_Omega0_m_classy;

      result = BEreq::class_get_Omega0_m();
    }

    /// Radiation
    void get_Omega0_r_classy(double& result)
    {
      using namespace Pipes::get_Omega0_r_classy;

      result = BEreq::class_get_Omega0_r();
    }

    /// Ultra-relativistic
    void get_Omega0_ur_classy(double& result)
    {
      using namespace Pipes::get_Omega0_ur_classy;

      result = BEreq::class_get_Omega0_ur();
    }

    /// Non-cold Dark Matter
    void get_Omega0_ncdm_tot_classy(double& result)
    {
      using namespace Pipes::get_Omega0_ncdm_tot_classy;

      result = BEreq::class_get_Omega0_ncdm_tot();
    }

    /// Sigma8
    void get_Sigma8_classy(double& result)
    {
      using namespace Pipes::get_Sigma8_classy;

      double sigma8 = BEreq::class_get_sigma8();
      double Omega0_m = *Dep::Omega0_m;

      result = sigma8*pow(Omega0_m/0.3, 0.5);
    }

    /// Effective number of neutrino species
    // (mostly for cross-checking!)
    void get_Neff_classy(double& result)
    {
      using namespace Pipes::get_Neff_classy;

      result = BEreq::class_get_Neff();
    }

    /// Comoving sound horizon at Baryon drag epoch
    void get_rs_drag_classy(double& result)
    {
      using namespace Pipes::get_rs_drag_classy;

      result = BEreq::class_get_rs();
    }


/***************/
/* MontePython */
/***************/

    /// function to fill the mcmc_parameters dictionary of MontePython's Data object with current 
    /// values of nuisance parameters
    void set_parameter_dict_for_MPLike(pybind11::dict & result)
    {
      using namespace Pipes::set_parameter_dict_for_MPLike;
      using namespace pybind11::literals;

      // The loop has to be executed for every parameter point. It takes about 0.00023 s -> ~4 minutes for 1e6 points 
      for (auto it=Param.begin(); it != Param.end(); it++)
      {
        std::string name = it->first;
        double value = *Param[name];

        // check if any models are scanned for which we had to rename the nuisance parameters due to 
        //    a) parameters having the same name  -> e.g. 'epsilon' and 'sigma_NL'
        //    b) parameter names containing symbols that can't be used in macros -> e.g. "^" in 'beta_0^Euclid'

        // a) have to rename parameters epsilon_ska, epsilon_euclid,.. to "epsilon" as they are implemented in MontePython
        if (name.find("epsilon") != std::string::npos){name="epsilon";}

        // a) have to rename parameters sigma_NL_ska, sigma_NL_euclid,.. to "sigma_NL" as they are implemented in MontePython
        else if (name.find("sigma_NL") != std::string::npos){name="sigma_NL";}

        // b) get the "^" characters back into the parameter names
        //   -> beta_x<experiment> has to be beta_x^<experiment> where x = 0 or 1 and <experiment> = Euclid, SKA1 or SKA2
        //if(ModelInUse("cosmo_nuisance_euclid_pk") or ModelInUse("cosmo_nuisance_ska"))
        else if (name.find("beta_")!= std::string::npos){name = name.insert(6,"^");}

        result[name.c_str()] = pybind11::dict("current"_a=value,"scale"_a=1.); // scale always 1 in GAMBIT
      }
    }

    /// function to fill the mcmc_parameters dictionary of MontePython's Data object with current 
    /// This version of the capability 'parameter_dict_for_MPLike' is used when no Likelihood with nuisance parameters is in use
    /// just passes an empty py dictionary
    void pass_empty_parameter_dict_for_MPLike(pybind11::dict & result)
    {
      using namespace Pipes::pass_empty_parameter_dict_for_MPLike;
      pybind11::dict r;
      result = r;
      // Nothing to do here.
    }

    /// Set the names of all experiments in use for scan
    /// basically reads in the runOptions 'Likelihoods' and/or 'Observables' of capability MP_experiment_names
    void set_MP_experiment_names(map_str_map_str_str & result)
    {
      using namespace Pipes::set_MP_experiment_names;

      static bool first = true;
      if (first)
      {
        // get list of likelihoods implemented in MP
        std::vector<str> avail_likes = BEreq::get_MP_available_likelihoods();

        // read in requested likelihoods
        YAML::Node MP_Likelihoods;
        if (runOptions->hasKey("Likelihoods")) MP_Likelihoods = runOptions->getValue<YAML::Node>("Likelihoods");

        // create string -> string map mapping the likelihood name to the '.data' file that will
        // be used to initialise likelihood object in MP. If you want to use the default one simply
        // set the the data-file string to "default"
        for (auto it=MP_Likelihoods.begin(); it != MP_Likelihoods.end(); it++)
        {
          std::string likelihood = it->first.as<std::string>();
          std::string data_file = it->second.as<std::string>();

          // check if requested likelihood is contained in vector of implemented likelihoods
          // if not throw error & list all available options
          if (std::find(avail_likes.begin(), avail_likes.end(), likelihood) == avail_likes.end())
          {

            str errmsg = "Likelihood '" + likelihood + "' is not implemented in MontePython. Check for typos or implement it.\nLikelihoods currently available are:\n"; 
            for(auto const& value: avail_likes)
            {
              errmsg += ("\t"+value+"\n");
            }
            CosmoBit_error().raise(LOCAL_INFO,errmsg);
          }
          else
          {
            result["Likelihoods"][likelihood] = data_file;
            logger() << LogTags::debug << "Read MontePythonLike option "<< likelihood << ", using data file " << data_file<< EOM;
          }
        }

        // & Observables
        YAML::Node MP_Observables;
        if (runOptions->hasKey("Observables")) MP_Observables = runOptions->getValue<YAML::Node>("Observables");

        for (auto it=MP_Observables.begin(); it != MP_Observables.end(); it++)
        {
          std::string obs = it->first.as<std::string>();
          std::string data_file = it->second.as<std::string>();

          // check if requested likelihood is contained in vector of implemented likelihoods
          // if not throw error & list all available options
          if (std::find(avail_likes.begin(), avail_likes.end(), obs) == avail_likes.end())
          {

            str errmsg = "Likelihood '" + obs + "' is not implemented in MontePython. Check for typos or implement it.\nLikelihoods currently available are:\n"; 
            for(auto const& value: avail_likes)
            {
              errmsg += ("\t"+value+"\n");
            }
            CosmoBit_error().raise(LOCAL_INFO,errmsg);
          }
          else
          {
            result["Observables"][obs] = data_file;
            logger() << LogTags::debug << "Read MontePythonLike option "<< obs << ", using data file " << data_file<< EOM;
          }
        }
        first = false;
      }
    }

    /// When initialising the MontePython Likelihood objects they add the output that needs to be computed by class
    /// to the input dictionary. We need to get these before starting the class run
    /// e.g. for Planck_SZ likelihood the entries {'output': ' mPk ', 'P_k_max_h/Mpc': '1.'} need to be added 
    /// to compute all needed observables
    void init_cosmo_args_from_MPLike(pybind11::dict &result)
    {
      using namespace Pipes::init_cosmo_args_from_MPLike;

      // CosmoBit::MPLike_data_container should only be created once when calculating the first point.
      // After that is has to be kept alive since it contains a vector with the initialised MPLike 
      // Likelihood objects.
      static bool first_run = true;
      static pybind11::object data;
      if(first_run)
      {
        // This is a map with the following structure:
        // { Likelihoods: {likelihood1: mode, likelihood2: mode...}
        //   Observables: {observable1: mode, observable2: mode...} }
        map_str_map_str_str experiment_names = *Dep::MP_experiment_names;

        // We need to pull the names of all Likelihoods AND observables
        // to initialise the data structures in MontePython
        map_str_str experiments;

        if (experiment_names.find("Likelihoods") != experiment_names.end())
        {
          for (auto it = experiment_names.at("Likelihoods").begin();
                    it != experiment_names.at("Likelihoods").end(); ++it)
          {
            // Add each experiment : datafile pair to the map_str_str.
            experiments[it->first] = it->second;
          }
        }
        if (experiment_names.find("Observables") != experiment_names.end())
        {
          for (auto it = experiment_names.at("Observables").begin();
                    it != experiment_names.at("Observables").end(); ++it)
          {
            // Add each experiment : datafile pair to the map_str_str.
            experiments[it->first] = it->second;
          }
        }

        logger() << LogTags::info << "(init_cosmo_args_from_MPLike) List of experiments: " << experiments << EOM;
        data = BEreq::create_MP_data_object(experiments);
        // add current parameters to data object to enable check if all nuisance parameters are 
        // scanned upon initialisation of likelihood objects
        data.attr("mcmc_parameters") = *Dep::parameter_dict_for_MPLike;
        map_str_pyobj likelihoods = BEreq::create_MP_likelihood_objects(data, experiments);
        first_run = false;
      }

      result.clear();
      pybind11::dict tmp_dict = data.attr("cosmo_arguments");
      // Stringify all values in the dictionary and strip off leading and trailing whitespaces
      for (auto it: tmp_dict)
      {
        std::string key = (pybind11::str(it.first)).cast<std::string>();
        std::string val = (pybind11::str(it.second)).cast<std::string>();
        boost::algorithm::trim(val);
        result[key.c_str()] = val.c_str();
      }
    }

    /// Computes lnL for each experiment initialised in MontePython
    void calc_MP_LogLikes(CosmoBit::MPLike_result_container & result)
    {
      using namespace Pipes::calc_MP_LogLikes;
      using namespace pybind11::literals;

      // A list of the experiments initialised in the YAML file
      // This is a map with the following structure:
      // { Likelihoods: {likelihood1: mode, likelihood2: mode...}
      //   Observables: {observable1: mode, observable2: mode...} }
      map_str_map_str_str experiment_names = *Dep::MP_experiment_names;

      // We need to pull the names of all Likelihoods AND observables
      // to initialise the data structures in MontePython
      map_str_str experiments;

      if (experiment_names.find("Likelihoods") != experiment_names.end())
      {
        for (auto it = experiment_names.at("Likelihoods").begin();
                  it != experiment_names.at("Likelihoods").end(); ++it)
        {
          // Add each experiment : datafile pair to the map_str_str.
          experiments[it->first] = it->second;
        }
      }
      if (experiment_names.find("Observables") != experiment_names.end())
      {
        for (auto it = experiment_names.at("Observables").begin();
                  it != experiment_names.at("Observables").end(); ++it)
        {
          // Add each experiment : datafile pair to the map_str_str.
          experiments[it->first] = it->second;
        }
      }
 
      // CosmoBit::MPLike_data_container should only be created once when calculating the first point.
      // After that is has to be kept alive since it contains a vector with the initialised MPLike 
      // Likelihood objects.
      pybind11::object data;
      map_str_pyobj likelihoods;
      static bool first_run = true;
      if(first_run)
      {
        data = BEreq::create_MP_data_object(experiments);
        // add current parameters to data object to enable check if all nuisance parameters are 
        // scanned upon initialisation of likelihood objects
        data.attr("mcmc_parameters") = *Dep::parameter_dict_for_MPLike;
        likelihoods = BEreq::create_MP_likelihood_objects(data, experiments);
        first_run = false;
      }

      static const CosmoBit::MPLike_data_container mplike_cont(data, likelihoods);

      // pass current values of nuisance parameters to data.mcmc_parameters dictionary for likelihood computation in MP
      mplike_cont.data.attr("mcmc_parameters") = *Dep::parameter_dict_for_MPLike;

      // Create instance of classy class Class
      pybind11::object cosmo = BEreq::get_classy_cosmo_object();

      // Loop through the list of experiments, and query the lnL from the
      // MontePython backend
      // Only if the experiment is requested as a Likelihood.
      // Separate loop below for observables.
      if (experiment_names.find("Likelihoods") != experiment_names.end())
      {
        logger() << LogTags::debug << "(calc_MP_LogLikes): Likelihoods\n\n";
        for (auto const& it : experiment_names.at("Likelihoods"))
        {
          // likelihood names are keys of experiment map (str, str map mapping likelihood name to .data file)
          std::string like_name = it.first;
          double logLike = BEreq::get_MP_loglike(mplike_cont, cosmo, like_name);
          result.add_logLike(like_name, logLike);
          logger() << "name: " << it.first << "\tvalue: " << logLike << "\n";
          //std::cout << "name: " << it.first << "\tvalue: " << logLike << std::endl;
        }
      }

      // Loop through the list of experiments, and query the lnL from the
      // MontePython backend.
      // ONLY if the experiment is requested as an observable.
      if (experiment_names.find("Observables") != experiment_names.end())
      {
        logger() << LogTags::debug << "(calc_MP_LogLikes): Observables\n\n";
        for (auto const& it : experiment_names.at("Observables"))
        {
          // likelihood names are keys of experiment map (str, str map mapping likelihood name to .data file)
          std::string like_name = it.first;
          double logLike = BEreq::get_MP_loglike(mplike_cont, cosmo, like_name);
          result.add_obs(like_name, logLike);
          logger() << "name: " << it.first << "\tvalue: " << logLike << "\n";
          //std::cout << "name: " << it.first << "\tvalue: " << logLike << std::endl;
        }
        logger() << EOM;
      }

    }

    /// Computes the combined lnL from the set of experiments
    /// given to MontePython.
    void calc_MP_combined_LogLike(double& result)
    {
      using namespace Pipes::calc_MP_combined_LogLike;

      // get string double map from MPlike result container mapping
      // experiment name to LogLike value for Likelihoods included 
      // into the scan
      CosmoBit::MPLike_result_container MPLike_results = *Dep::calc_MP_LogLikes;
      map_str_dbl MP_lnLs = MPLike_results.get_logLike_results();

      // Iterate through map of doubles and return one big fat double.
      double lnL = 0.;

      logger() << LogTags::debug << "(calc_MP_combined_LogLike):\n\n";
      for (const auto &p : MP_lnLs)
      {
          logger()  << "name: "  << p.first << "\tvalue: " << p.second << "\n";
        lnL += p.second;
      }
      logger() << EOM;

      result = lnL;
    }
    /// Computes lnL for each experiment initialised in MontePython
    /// but DOES NOT add it to the compound lnL in GAMBIT.
    /// This should be used when one wishes to grab a Likelihood from
    /// MontePython but does not want it to steer the scans
    /// (e.g. for forecasting; conflicting likelihoods; etc.) -- use with caution!
    void get_MP_LogLikes(map_str_dbl & result)
    {
      using namespace Pipes::get_MP_LogLikes;
      using namespace pybind11::literals;

      // get string double map from MPlike result container mapping
      // experiment name to LogLike value for Observables included 
      // into the scan
      CosmoBit::MPLike_result_container MPLike_results = *Dep::calc_MP_LogLikes;
      result = MPLike_results.get_logLike_results();
    }

    /// Computes lnL for each experiment initialised in MontePython
    /// but DOES NOT add it to the compound lnL in GAMBIT.
    /// This should be used when one wishes to grab a Likelihood from
    /// MontePython but does not want it to steer the scans
    /// (e.g. for forecasting; conflicting likelihoods; etc.) -- use with caution!
    void get_MP_Observables(map_str_dbl & result)
    {
      using namespace Pipes::get_MP_Observables;
      using namespace pybind11::literals;

      // get string double map from MPlike result container mapping
      // experiment name to LogLike value for Observables included 
      // into the scan
      CosmoBit::MPLike_result_container MPLike_results = *Dep::calc_MP_LogLikes;
      result = MPLike_results.get_obs_results();
    }
  }
}
