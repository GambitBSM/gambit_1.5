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
///  \date 2019 Jan - Feb
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 June
///  \date 2019 Mar,June
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date 2019 June
///
///  *********************************************
#include <cmath>
#include <functional>
#include <iostream>
#include <omp.h>
#include <stdlib.h>     /* malloc, free, rand */
#include <string>
#include <valarray>

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

    double r_SR(double eps,double eta)
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

    void injection_spectrum_ToyModel(DarkAges::injectionSpectrum& spectrum)
    {
      using namespace Pipes::injection_spectrum_ToyModel;

      double m = *Param["mass"];
      double BR_el = *Param["BR"];
      double BR_ph = 1.0 - BR_el;

      spectrum.E.clear();
      spectrum.spec_el.clear();
      spectrum.spec_ph.clear();

      spectrum.E.resize(1,m*0.5);
      spectrum.spec_el.resize(1,BR_el*2e9);
      spectrum.spec_ph.resize(1,BR_ph*2e9);
    }

    void lifetime_CosmoALP(double& result)
    {
      // lifetime in s
      using namespace Pipes::lifetime_CosmoALP;
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

    void minimum_abundance_CosmoALP(double& result)
    {
      using namespace Pipes::minimum_abundance_CosmoALP;

      // Read T_R in MeV and convert it to GeV
      static double T_R_in_GeV = 1e-3*(runOptions->getValueOrDef<double>(5.,"T_R")); // Read T_R in MeV and convert it to GeV

      // check for stupid input (T_R < m_e)
      // Throw an error if the user really pushed it that far.
      if (m_electron >= T_R_in_GeV)
        CosmoBit_error().raise(LOCAL_INFO,"The reheating temperature is below the electron mass.");

      double gagg = *Param["gagg"]; // in GeV^-1
      double Ya0_min = 1.56e-5 * pow(gagg,2) * m_planck * (T_R_in_GeV - m_electron);

      result = Ya0_min;
    }

    void DM_fraction_CosmoALP(double& result)
    {
      using namespace Pipes::DM_fraction_CosmoALP;
      double Ya0 = *Param["Ya0"]; // ALP abundance (after production)
      double ma0 = *Param["ma0"]; // non-thermal ALP mass in eV
      double T = *Dep::T_cmb; // CMB temperature in K
      double ssm0 = entropy_density_SM(T); // SM entropy density today in eV^3 (cf. footnote 24 of PDG2018-Astrophysical parameters)
      double rho0_a = Ya0 * ma0 * ssm0; // energy density of axions today in eV^4

      double omega_cdm = *Param["omega_cdm"]; // omega_cdm = Omega_cdm * h^2
      const double rho0_crit_by_h2 = 3.*pow(m_planck_red*1e9,2) * pow((1e5*1e9*hbar/_Mpc_SI_),2); // rho0_crit/(h^2)
      double rho0_cdm = omega_cdm * rho0_crit_by_h2; // rho0_cdm = Omega_cdm * rho0_crit;

      double xi_ini = rho0_a / rho0_cdm; // Initial fraction of decauying dark matter (rho_dcdm / rho_cdm) and the naive value if decay of the ALP would be neglected.
      double tau_a = *Dep::lifetime;
      const double t_rec = 1e12;
      double xi_at_rec = xi_ini * exp(-t_rec/tau_a );

      // invalidate if there are more ALPs than dark matter at the time of recombination (t ~ 1e12s)
      if (xi_at_rec  > 1.)
      {
        std::ostringstream err;
        err << "ALPs are over-abundant (n_a > n_cdm) at t = 10^12 s. (n_a/n_cdm = "<< xi_at_rec <<")";
        invalid_point().raise(err.str());
      }

      // TODO: Make the part which follows down below a proper likelihood

      // Consistency check if value for xi (rho0_a / rh0_cdm) is in contradiction with
      // T_R > 5 MeV (or user input). If so, invalidate the point.
      // (Reason: T_R predicts a minimum abundance Ya0_min through Primakoff processes. Any additional production mechanism would only increase Ya0)
      double Ya0_min = *Dep::minimum_abundance;
      if (Ya0 < Ya0_min)
      {
        std::ostringstream err;
        err << "The choice of Ya0 (";
        err << Ya0;
        err << ") is in contradiction with the minimum abundance Ya0_min (";
        err << Ya0_min;
        err << ") produced via Primakoff processes.";
        invalid_point().raise(err.str());
      }

      result = xi_ini;
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
      double z_eff = runOptions->getValueOrDef<double>(600.,"z_eff");

      DarkAges::fz_table fzt = *Dep::energy_injection_efficiency;
      std::vector<double> z = fzt.redshift;
      std::vector<double> fh = fzt.f_heat;
      std::vector<double> fly = fzt.f_lya;
      std::vector<double> fhi = fzt.f_hion;
      std::vector<double> fhei = fzt.f_heion;
      std::vector<double> flo = fzt.f_lowe;

      int npts = z.size();
      double ftot[npts];
      double red[npts];
      for (unsigned int i = 0; i < npts; i++)
      {
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
        std::cout << "tau = " << *Param["lifetime"] << std::endl;
        std::cout << "m = " << *Param["mass"] << std::endl;
        std::cout << "BR (electrom) = " << *Param["BR"] << std::endl;
        std::cout << "---------------" << std::endl;
        std::cout << "z\tf_heat\tf_lya\tf_hion\tf_heion\tf_lowe" << std::endl;
        for (unsigned int i = z.size() - last_steps; i < z.size(); i++)
        {
          std::cout << z.at(i) << "\t" << fh.at(i) << "\t" << fly.at(i) << "\t" << fhi.at(i) << "\t" << fhei.at(i) << "\t" << flo.at(i)  << std::endl;
        }
        std::cout << "f_eff (sum of all channels at z = "<< z_eff << ") = " << result << std::endl;
        std::cout << "################\n" << std::endl;
      }
    }

    void lnL_A_planck_gaussian(double& result)
    {
      using namespace Pipes::lnL_A_planck_gaussian;

      double obs_mean = runOptions->getValueOrDef<double>(1.,"mean");
      double obs_err = runOptions->getValueOrDef<double>(0.0025,"sigma") ;
      double pred_mean = *Param["A_planck"];
      double pred_err = 0.0;

      result = Stats::gaussian_loglikelihood(pred_mean, obs_mean, pred_err, obs_err, false);
    }


    void set_NuMasses_SM(map_str_dbl &result)
    {
      using namespace Pipes::set_NuMasses_SM;
      
      double mNu1, mNu2, mNu3;
      int N_ncdm = 0;
      if (ModelInUse("StandardModel_SLHA2"))
      {
        // (PS) Heads up! The units in StandardModel_SLHA2 are GeV
        // Here we are using eV
        mNu1 = 1e9*(*Param["mNu1"]);
        mNu2 = 1e9*(*Param["mNu2"]);
        mNu3 = 1e9*(*Param["mNu3"]);
      }
      else
      {
        mNu1 = 0.06;
        mNu2 = 0.;
        mNu3 = 0.;
      }

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

      switch (N_ncdm)
      {
        case 1:
          result["N_ur_SMnu"]= 2.0328;  // dNeff= 2.0328 for 1 massive neutrino at CMB release
          break;
        case 2:
          result["N_ur_SMnu"]= 1.0196;  // dNeff= 2.0328 for 1 massive neutrino at CMB release
          break;
        case 3:
          result["N_ur_SMnu"]= 0.00641;  // dNeff= 0.00641 for 3 massive neutrinos at CMB release
          break;
        case 0:
          {
            std::ostringstream err;
            err << "Warning: All your neutrino masses are zero. The Planck baseline LCDM model assumes at least one massive neutrino.\n";
            err << "A cosmological without massive neutrinos is not implemented in CosmoBit. If you want to consider this you can add it to the function 'set_NuMasses_SM' of the capability 'NuMasses_SM' ";
            CosmoBit_warning().raise(LOCAL_INFO, err.str());
          }
          break;
        default:
          {
            std::ostringstream err;
            err << "You are asking for more than three massive neutrino species.\n";
            err << "Such a case is not implemented here in CosmoBit. But if you know what you are doing just add the respective case in \'class_set_parameter_LCDM_family\' and you are done.";
            CosmoBit_error().raise(LOCAL_INFO, err.str());
          }
      }
    }

    void class_set_parameter_LCDM_family(Class_container& cosmo)
    {
      using namespace Pipes::class_set_parameter_LCDM_family;

      int l_max=cosmo.lmax,ii=1;

      cosmo.input.clear();

      map_str_dbl NuMasses_SM = *Dep::NuMasses_SM;
      int N_ncdm = (int)NuMasses_SM["N_ncdm"];

      if (N_ncdm > 0.)
      { 
        cosmo.input.addEntry("N_ncdm",N_ncdm);
        cosmo.input.addEntry("m_ncdm",m_ncdm_classInput(NuMasses_SM));

        std::vector<double> T_ncdm(N_ncdm,*Dep::T_ncdm);
        cosmo.input.addEntry("T_ncdm", T_ncdm);
      }
      else
      {
        cosmo.input.addEntry("T_ncdm", *Dep::T_ncdm); 
      }

      cosmo.input.addEntry("N_ur",*Dep::class_Nur);

      cosmo.input.addEntry("output","tCl pCl lCl");
      cosmo.input.addEntry("l_max_scalars",l_max);
      cosmo.input.addEntry("lensing","yes");

      cosmo.input.addEntry("T_cmb",*Dep::T_cmb);
      cosmo.input.addEntry("omega_b",*Param["omega_b"]);
      cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
      cosmo.input.addEntry("H0",*Param["H0"]);
      cosmo.input.addEntry("ln10^{10}A_s",*Param["ln10A_s"]);
      cosmo.input.addEntry("n_s",*Param["n_s"]);
      cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);

      if (ModelInUse("TestDecayingDM"))  // TODO: need to test if class or exo_class in use! does not work 
      {
        cosmo.input.addEntry("energy_deposition_function","GAMBIT");
        cosmo.input.addEntry("tau_dcdm",*Dep::lifetime);
        cosmo.input.addEntry("decay_fraction",*Dep::DM_fraction);
      }

      std::vector<double> Helium_abund = *Dep::Helium_abundance; // .at(0): mean, .at(1): uncertainty
      cosmo.input.addEntry("YHe",Helium_abund.at(0));

      YAML::Node class_dict;
      if (runOptions->hasKey("class_dict"))
      {
        class_dict = runOptions->getValue<YAML::Node>("class_dict");
        for (auto it=class_dict.begin(); it != class_dict.end(); it++)
        {
          std::string name = it->first.as<std::string>();
          std::string value = it->second.as<std::string>();
          cosmo.input.addEntry(name,value);
        }
      }

      std::string log_msg = cosmo.input.print_entries_to_logger(); // print all parameters and their values that are passed to class in debug.log
      logger() << LogTags::debug << log_msg << EOM;
    }

    void class_set_parameter_LCDM_SingletDM(Class_container& cosmo)
    {
      using namespace Pipes::class_set_parameter_LCDM_SingletDM;

      int l_max = cosmo.lmax;

      double sigmav = *Dep::sigmav; // in cm^3 s^-1
      double mass = *Dep::mwimp; // in GeV
      double feff = runOptions->getValueOrDef<double>(1.,"f_eff");
      double annihilation = (1.0/1.78e-21)*(sigmav/mass)*feff; // in m^3 s^-1 kg^-1

      cosmo.input.clear();

      cosmo.input.addEntry("output","tCl pCl lCl");
      cosmo.input.addEntry("l_max_scalars",l_max);
      cosmo.input.addEntry("lensing","yes");

      cosmo.input.addEntry("T_cmb",*Dep::T_cmb);
      cosmo.input.addEntry("omega_b",*Param["omega_b"]);
      cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
      cosmo.input.addEntry("H0",*Param["H0"]);
      cosmo.input.addEntry("ln10^{10}A_s",*Param["ln10A_s"]);
      cosmo.input.addEntry("n_s",*Param["n_s"]);
      cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
      cosmo.input.addEntry("annihilation",annihilation);

      YAML::Node class_dict;
      if (runOptions->hasKey("class_dict"))
      {
        class_dict = runOptions->getValue<YAML::Node>("class_dict");
        for (auto it=class_dict.begin(); it != class_dict.end(); it++)
        {
          std::string name = it->first.as<std::string>();
          std::string value = it->second.as<std::string>();
          cosmo.input.addEntry(name,value);
        }
      }
/*
      std::cout << "omega_b = " << *Param["omega_b"] << std::endl;
      std::cout << "omega_cdm = " << *Param["omega_cdm"] << std::endl;
      std::cout << "H0 = " << *Param["H0"] << std::endl;
      std::cout << "ln10A_s = " << *Param["ln10A_s"] << std::endl;
      std::cout << "n_s = " << *Param["n_s"] << std::endl;
      std::cout << "tau_reio = " << *Param["tau_reio"] << std::endl;
      std::cout << "annihilation = " << annihilation << "(Resulting from m = " << *Param["mS"] << " and lambda = "<< *Param["lambda_hS"] << ")" << std::endl;
*/
    }

    void class_set_parameter_LCDMtensor(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_LCDMtensor" << std::endl;
      using namespace Pipes::class_set_parameter_LCDMtensor;

      int l_max = cosmo.lmax;

      cosmo.input.clear();

      cosmo.input.addEntry("output","tCl pCl lCl");
      cosmo.input.addEntry("l_max_scalars",l_max);
      cosmo.input.addEntry("lensing","yes");
      cosmo.input.addEntry("modes","s,t");

      cosmo.input.addEntry("T_cmb",*Dep::T_cmb);
      cosmo.input.addEntry("omega_b",*Param["omega_b"]);
      cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
      cosmo.input.addEntry("H0",*Param["H0"]);
      cosmo.input.addEntry("ln10^{10}A_s",*Param["ln10A_s"]);
      cosmo.input.addEntry("n_s",*Param["n_s"]);
      cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
      cosmo.input.addEntry("r",*Param["r_tensor"]);

      YAML::Node class_dict;
      if (runOptions->hasKey("class_dict"))
      {
        class_dict = runOptions->getValue<YAML::Node>("class_dict");
        for (auto it=class_dict.begin(); it != class_dict.end(); it++)
        {
          std::string name = it->first.as<std::string>();
          std::string value = it->second.as<std::string>();
          cosmo.input.addEntry(name,value);
        }
      }
/*
      std::cout << "omega_b = " << *Param["omega_b"] << std::endl;
      std::cout << "omega_cdm = " << *Param["omega_cdm"] << std::endl;
      std::cout << "H0 = " << *Param["H0"] << std::endl;
      std::cout << "ln10A_s = " << *Param["ln10A_s"] << std::endl;
      std::cout << "n_s = " << *Param["n_s"] << std::endl;
      std::cout << "tau_reio = " << *Param["tau_reio"] << std::endl;
      std::cout << "r_tensor = " << *Param["r_tensor"] << std::endl;
*/
    }

    void class_set_parameter_inf_SR1quad_LCDMt(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_inf_SR1quad_LCDMt" << std::endl;
      using namespace Pipes::class_set_parameter_inf_SR1quad_LCDMt;

      int l_max=cosmo.lmax;

      cosmo.input.clear();

      cosmo.input.addEntry("T_cmb",*Dep::T_cmb);
      cosmo.input.addEntry("output","tCl pCl lCl");
      cosmo.input.addEntry("l_max_scalars",l_max);
      cosmo.input.addEntry("modes","s,t");
      cosmo.input.addEntry("lensing","yes");

      // Parameters to be passed to the potential                       //
      //---------------------------------------------------------------//
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");
      //---------------------------------------------------------------//
      // Ideally what's between this line and the line (*) indicated    //
      // below is a general enough example of needed to set a new       //
      // inflationary model with MultiModeCode, i.e. the model          //
      // parameters that's to be sampled.                               //
      //---------------------------------------------------------------//
      //--------------- Sampling m2_inflaton  -------------------------//
      //---------------------------------------------------------------//
      vparams[0] = *Param["m2_inflaton"];
      //---------------------------------------------------------------//
      //--------------- Sampling N_pivot ------------------------------//
      //---------------------------------------------------------------//
      double N_pivot = *Param["N_pivot"];
      //---------------------------------------------------------------//
      //------------------------------(*)------------------------------//



      //-------------------------------------------------------------
      //  Below we set settings for inflation solver MultiModecode.
      //-------------------------------------------------------------
      int silence = runOptions->getValue<int> ("is_not_silent");
      //-------------------------------------------------------------
      // Initialization parameters controlling main characteristics.
      //-------------------------------------------------------------
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");
      //-------------------------------------------------------------
      // Control the output of analytic approximations for comparison.
      //-------------------------------------------------------------
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");
      //-------------------------------------------------------------
      // Parameters to control how the ICs are sampled.
      //-------------------------------------------------------------
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");
      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this
      //-------------------------------------------------------------
      // Parameters to control how the vparams are sampled.
      //-------------------------------------------------------------
      int param_sampling = runOptions->getValue<int> ("param_sampling");
      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");
      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");
      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");
      //-------------------------------------------------------------
      // Priors on the IC and N_pivot ranges
      //-------------------------------------------------------------
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");
      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");
      //-------------------------------------------------------------
      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      //-------------------------------------------------------------
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");
      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");



      gambit_inflation_observables observs;

      //-------------------------------------------------------------
      // The function below calls the MultiModeCode backend
      //  for a given choice of inflationary model,
      //  which calculates the observables.
      //-------------------------------------------------------------
      BEreq::multimodecode_gambit_driver(&observs,
                                         num_inflaton,
                                         potential_choice,
                                         slowroll_infl_end,
                                         instreheat,
                                         vparam_rows,
                                         use_deltaN_SR,
                                         evaluate_modes,
                                         use_horiz_cross_approx,
                                         get_runningofrunning,
                                         ic_sampling,
                                         energy_scale,
                                         numb_samples,
                                         save_iso_N,
                                         N_iso_ref,
                                         param_sampling,
                                         byVal(&vp_prior_min[0]),
                                         byVal(&vp_prior_max[0]),
                                         varying_N_pivot,
                                         use_first_priorval,
                                         byVal(&phi_init0[0]),
                                         byVal(&dphi_init0[0]),
                                         byVal(&vparams[0]),
                                         N_pivot,
                                         k_pivot,
                                         dlnk,
                                         calc_full_pk,
                                         steps,
                                         kmin,
                                         byVal(&phi0_priors_min[0]),
                                         byVal(&phi0_priors_max[0]),
                                         byVal(&dphi0_priors_min[0]),
                                         byVal(&dphi0_priors_max[0]),
                                         N_pivot_prior_min,
                                         N_pivot_prior_max);

      if (calc_full_pk == 0)
      {
        if (silence == 1){
          std::cout << "A_s =" << observs.As << std::endl;
          // std::cout << "A_iso " << observs.A_iso << std::endl;
          // std::cout << "A_pnad " << observs.A_pnad << std::endl;
          // std::cout << "A_end " << observs.A_ent << std::endl;
          // std::cout << "A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
          std::cout << "n_s =" << observs.ns << std::endl;
          // std::cout << "n_t " << observs.nt << std::endl;
          // std::cout << "n_iso " << observs.n_iso << std::endl;
          // std::cout << "n_pnad " << observs.n_pnad << std::endl;
          // std::cout << "n_ent " << observs.n_ent << std::endl;
          std::cout << "r =" << observs.r << std::endl;
          std::cout << "alpha_s =" << observs.alpha_s << std::endl;
          // std::cout << "runofrun " << observs.runofrun << std::endl;
          // std::cout << "f_NL " << observs.f_NL << std::endl;
          // std::cout << "tau_NL " << observs.tau_NL << std::endl;
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("A_s",observs.As);
        cosmo.input.addEntry("n_s",observs.ns);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
        cosmo.input.addEntry("r",observs.r);
      }
      else
      {
        /* gambit_Pk type is set if full_spectra is asked. */
        cosmo.input.addEntry("P_k_ini type","gambit_Pk");

        cosmo.Pk_S.resize(steps+1, 0.);
        cosmo.Pk_T.resize(steps+1, 0.);
        cosmo.k_ar.resize(steps+1, 0.);

        for (int ii=0; ii < steps; ii++)
        {
          cosmo.k_ar.at(ii) = observs.k_array[ii];
          cosmo.Pk_S.at(ii) = observs.pks_array[ii];
          cosmo.Pk_T.at(ii) = observs.pkt_array[ii];
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
      }
      YAML::Node class_dict;
      if (runOptions->hasKey("class_dict"))
      {
        class_dict = runOptions->getValue<YAML::Node>("class_dict");
        for (auto it=class_dict.begin(); it != class_dict.end(); it++)
        {
          std::string name = it->first.as<std::string>();
          std::string value = it->second.as<std::string>();
          cosmo.input.addEntry(name,value);
        }
      }
    }

    void class_set_parameter_inf_1quarInf_LCDMt(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_inf_SR1quad_LCDMt" << std::endl;
      using namespace Pipes::class_set_parameter_inf_1quarInf_LCDMt;

      int l_max=cosmo.lmax;

      cosmo.input.clear();

      cosmo.input.addEntry("output","tCl pCl lCl");
      cosmo.input.addEntry("l_max_scalars",l_max);
      cosmo.input.addEntry("modes","s,t");
      cosmo.input.addEntry("lensing","yes");

      cosmo.input.addEntry("T_cmb",*Dep::T_cmb);
      // Parameters to be passed to the potential                       //
      //---------------------------------------------------------------//
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");
      //---------------------------------------------------------------//
      // Ideally what's between this line and the line (*) indicated    //
      // below is a general enough example of needed to set a new       //
      // inflationary model with MultiModeCode, i.e. the model          //
      // parameters that's to be sampled.                               //
      //---------------------------------------------------------------//
      //--------------- Sampling \lambda ------------------------------//
      //---------------------------------------------------------------//
      vparams[0] = *Param["lambda"];
      //---------------------------------------------------------------//
      //--------------- Sampling N_pivot ------------------------------//
      //---------------------------------------------------------------//
      double N_pivot = *Param["N_pivot"];
      //---------------------------------------------------------------//
      //------------------------------(*)------------------------------//



      //-------------------------------------------------------------
      //  Below we set settings for inflation solver MultiModecode.
      //-------------------------------------------------------------
      int silence = runOptions->getValue<int> ("is_not_silent");
      //-------------------------------------------------------------
      // Initialization parameters controlling main characteristics.
      //-------------------------------------------------------------
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");
      //-------------------------------------------------------------
      // Control the output of analytic approximations for comparison.
      //-------------------------------------------------------------
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");
      //-------------------------------------------------------------
      // Parameters to control how the ICs are sampled.
      //-------------------------------------------------------------
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");
      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this
      //-------------------------------------------------------------
      // Parameters to control how the vparams are sampled.
      //-------------------------------------------------------------
      int param_sampling = runOptions->getValue<int> ("param_sampling");
      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");
      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");
      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");
      //-------------------------------------------------------------
      // Priors on the IC and N_pivot ranges
      //-------------------------------------------------------------
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");
      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");
      //-------------------------------------------------------------
      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      //-------------------------------------------------------------
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");
      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");



      gambit_inflation_observables observs;

      //-------------------------------------------------------------
      // The function below calls the MultiModeCode backend
      //  for a given choice of inflationary model,
      //  which calculates the observables.
      //-------------------------------------------------------------
      BEreq::multimodecode_gambit_driver(&observs,
                                         num_inflaton,
                                         potential_choice,
                                         slowroll_infl_end,
                                         instreheat,
                                         vparam_rows,
                                         use_deltaN_SR,
                                         evaluate_modes,
                                         use_horiz_cross_approx,
                                         get_runningofrunning,
                                         ic_sampling,
                                         energy_scale,
                                         numb_samples,
                                         save_iso_N,
                                         N_iso_ref,
                                         param_sampling,
                                         byVal(&vp_prior_min[0]),
                                         byVal(&vp_prior_max[0]),
                                         varying_N_pivot,
                                         use_first_priorval,
                                         byVal(&phi_init0[0]),
                                         byVal(&dphi_init0[0]),
                                         byVal(&vparams[0]),
                                         N_pivot,
                                         k_pivot,
                                         dlnk,
                                         calc_full_pk,
                                         steps,
                                         kmin,
                                         byVal(&phi0_priors_min[0]),
                                         byVal(&phi0_priors_max[0]),
                                         byVal(&dphi0_priors_min[0]),
                                         byVal(&dphi0_priors_max[0]),
                                         N_pivot_prior_min,
                                         N_pivot_prior_max);

      if (calc_full_pk == 0)
      {
        if (silence == 1){
          std::cout << "A_s =" << observs.As << std::endl;
          // std::cout << "A_iso " << observs.A_iso << std::endl;
          // std::cout << "A_pnad " << observs.A_pnad << std::endl;
          // std::cout << "A_end " << observs.A_ent << std::endl;
          // std::cout << "A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
          std::cout << "n_s =" << observs.ns << std::endl;
          // std::cout << "n_t " << observs.nt << std::endl;
          // std::cout << "n_iso " << observs.n_iso << std::endl;
          // std::cout << "n_pnad " << observs.n_pnad << std::endl;
          // std::cout << "n_ent " << observs.n_ent << std::endl;
          std::cout << "r =" << observs.r << std::endl;
          std::cout << "alpha_s =" << observs.alpha_s << std::endl;
          // std::cout << "runofrun " << observs.runofrun << std::endl;
          // std::cout << "f_NL " << observs.f_NL << std::endl;
          // std::cout << "tau_NL " << observs.tau_NL << std::endl;
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("A_s",observs.As);
        cosmo.input.addEntry("n_s",observs.ns);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
        cosmo.input.addEntry("r",observs.r);
      }
      else
      {
        /* gambit_Pk type is set if full_spectra is asked. */
        cosmo.input.addEntry("P_k_ini type","gambit_Pk");

        cosmo.Pk_S.resize(steps+1, 0.);
        cosmo.Pk_T.resize(steps+1, 0.);
        cosmo.k_ar.resize(steps+1, 0.);

        for (int ii=0; ii < steps; ii++)
        {
          cosmo.k_ar.at(ii) = observs.k_array[ii];
          cosmo.Pk_S.at(ii) = observs.pks_array[ii];
          cosmo.Pk_T.at(ii) = observs.pkt_array[ii];
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
      }

      YAML::Node class_dict;
      if (runOptions->hasKey("class_dict"))
      {
        class_dict = runOptions->getValue<YAML::Node>("class_dict");
        for (auto it=class_dict.begin(); it != class_dict.end(); it++)
        {
          std::string name = it->first.as<std::string>();
          std::string value = it->second.as<std::string>();
          cosmo.input.addEntry(name,value);
        }
      }
    }

    void class_set_parameter_inf_1mono32Inf_LCDMt(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_inf_SR1quad_LCDMt" << std::endl;
      using namespace Pipes::class_set_parameter_inf_1mono32Inf_LCDMt;

      int l_max=cosmo.lmax;

      cosmo.input.clear();

      cosmo.input.addEntry("output","tCl pCl lCl");
      cosmo.input.addEntry("l_max_scalars",l_max);
      cosmo.input.addEntry("modes","s,t");
      cosmo.input.addEntry("lensing","yes");

      cosmo.input.addEntry("T_cmb",*Dep::T_cmb);
      // Parameters to be passed to the potential                       //
      //---------------------------------------------------------------//
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");
      //---------------------------------------------------------------//
      // Ideally what's between this line and the line (*) indicated    //
      // below is a general enough example of needed to set a new       //
      // inflationary model with MultiModeCode, i.e. the model          //
      // parameters that's to be sampled.                               //
      //---------------------------------------------------------------//
      //---------------------------------------------------------------//
      //--------------- Sampling \lambda ------------------------------//
      //---------------------------------------------------------------//
      vparams[0] = *Param["lambda"];
      //---------------------------------------------------------------//
      //--------------- Sampling N_pivot ------------------------------//
      //---------------------------------------------------------------//
      double N_pivot = *Param["N_pivot"];
      //---------------------------------------------------------------//
      //------------------------------(*)------------------------------//



      //-------------------------------------------------------------
      //  Below we set settings for inflation solver MultiModecode.
      //-------------------------------------------------------------
      int silence = runOptions->getValue<int> ("is_not_silent");
      //-------------------------------------------------------------
      // Initialization parameters controlling main characteristics.
      //-------------------------------------------------------------
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");
      //-------------------------------------------------------------
      // Control the output of analytic approximations for comparison.
      //-------------------------------------------------------------
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");
      //-------------------------------------------------------------
      // Parameters to control how the ICs are sampled.
      //-------------------------------------------------------------
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");
      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this
      //-------------------------------------------------------------
      // Parameters to control how the vparams are sampled.
      //-------------------------------------------------------------
      int param_sampling = runOptions->getValue<int> ("param_sampling");
      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");
      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");
      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");
      //-------------------------------------------------------------
      // Priors on the IC and N_pivot ranges
      //-------------------------------------------------------------
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");
      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");
      //-------------------------------------------------------------
      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      //-------------------------------------------------------------
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");
      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");



      gambit_inflation_observables observs;

      //-------------------------------------------------------------
      // The function below calls the MultiModeCode backend
      //  for a given choice of inflationary model,
      //  which calculates the observables.
      //-------------------------------------------------------------
      BEreq::multimodecode_gambit_driver(&observs,
                                         num_inflaton,
                                         potential_choice,
                                         slowroll_infl_end,
                                         instreheat,
                                         vparam_rows,
                                         use_deltaN_SR,
                                         evaluate_modes,
                                         use_horiz_cross_approx,
                                         get_runningofrunning,
                                         ic_sampling,
                                         energy_scale,
                                         numb_samples,
                                         save_iso_N,
                                         N_iso_ref,
                                         param_sampling,
                                         byVal(&vp_prior_min[0]),
                                         byVal(&vp_prior_max[0]),
                                         varying_N_pivot,
                                         use_first_priorval,
                                         byVal(&phi_init0[0]),
                                         byVal(&dphi_init0[0]),
                                         byVal(&vparams[0]),
                                         N_pivot,
                                         k_pivot,
                                         dlnk,
                                         calc_full_pk,
                                         steps,
                                         kmin,
                                         byVal(&phi0_priors_min[0]),
                                         byVal(&phi0_priors_max[0]),
                                         byVal(&dphi0_priors_min[0]),
                                         byVal(&dphi0_priors_max[0]),
                                         N_pivot_prior_min,
                                         N_pivot_prior_max);

      if (calc_full_pk == 0)
      {
        if (silence == 1){
          std::cout << "A_s =" << observs.As << std::endl;
          // std::cout << "A_iso " << observs.A_iso << std::endl;
          // std::cout << "A_pnad " << observs.A_pnad << std::endl;
          // std::cout << "A_end " << observs.A_ent << std::endl;
          // std::cout << "A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
          std::cout << "n_s =" << observs.ns << std::endl;
          // std::cout << "n_t " << observs.nt << std::endl;
          // std::cout << "n_iso " << observs.n_iso << std::endl;
          // std::cout << "n_pnad " << observs.n_pnad << std::endl;
          // std::cout << "n_ent " << observs.n_ent << std::endl;
          std::cout << "r =" << observs.r << std::endl;
          std::cout << "alpha_s =" << observs.alpha_s << std::endl;
          // std::cout << "runofrun " << observs.runofrun << std::endl;
          // std::cout << "f_NL " << observs.f_NL << std::endl;
          // std::cout << "tau_NL " << observs.tau_NL << std::endl;
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("A_s",observs.As);
        cosmo.input.addEntry("n_s",observs.ns);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
        cosmo.input.addEntry("r",observs.r);
      }
      else
      {
        /* gambit_Pk type is set if full_spectra is asked. */
        cosmo.input.addEntry("P_k_ini type","gambit_Pk");

        cosmo.Pk_S.resize(steps+1, 0.);
        cosmo.Pk_T.resize(steps+1, 0.);
        cosmo.k_ar.resize(steps+1, 0.);

        for (int ii=0; ii < steps; ii++)
        {
          cosmo.k_ar.at(ii) = observs.k_array[ii];
          cosmo.Pk_S.at(ii) = observs.pks_array[ii];
          cosmo.Pk_T.at(ii) = observs.pkt_array[ii];
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
      }

      YAML::Node class_dict;
      if (runOptions->hasKey("class_dict"))
      {
        class_dict = runOptions->getValue<YAML::Node>("class_dict");
        for (auto it=class_dict.begin(); it != class_dict.end(); it++)
        {
          std::string name = it->first.as<std::string>();
          std::string value = it->second.as<std::string>();
          cosmo.input.addEntry(name,value);
        }
      }
    }

    void class_set_parameter_inf_1linearInf_LCDMt(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_inf_SR1quad_LCDMt" << std::endl;
      using namespace Pipes::class_set_parameter_inf_1linearInf_LCDMt;

      int l_max=cosmo.lmax;

      cosmo.input.clear();

      cosmo.input.addEntry("output","tCl pCl lCl");
      cosmo.input.addEntry("l_max_scalars",l_max);
      cosmo.input.addEntry("modes","s,t");
      cosmo.input.addEntry("lensing","yes");

      cosmo.input.addEntry("T_cmb",*Dep::T_cmb);
      // Parameters to be passed to the potential                       //
      //---------------------------------------------------------------//
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");
      //---------------------------------------------------------------//
      // Ideally what's between this line and the line (*) indicated    //
      // below is a general enough example of needed to set a new       //
      // inflationary model with MultiModeCode, i.e. the model          //
      // parameters that's to be sampled.                               //
      //---------------------------------------------------------------//
      //--------------- Sampling \lambda ------------------------------//
      //---------------------------------------------------------------//
      vparams[0] = *Param["lambda"];
      //---------------------------------------------------------------//
      //--------------- Sampling N_pivot ------------------------------//
      //---------------------------------------------------------------//
      double N_pivot = *Param["N_pivot"];
      //---------------------------------------------------------------//
      //------------------------------(*)------------------------------//



      //-------------------------------------------------------------
      //  Below we set settings for inflation solver MultiModecode.
      //-------------------------------------------------------------
      int silence = runOptions->getValue<int> ("is_not_silent");
      //-------------------------------------------------------------
      // Initialization parameters controlling main characteristics.
      //-------------------------------------------------------------
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");
      //-------------------------------------------------------------
      // Control the output of analytic approximations for comparison.
      //-------------------------------------------------------------
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");
      //-------------------------------------------------------------
      // Parameters to control how the ICs are sampled.
      //-------------------------------------------------------------
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");
      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this
      //-------------------------------------------------------------
      // Parameters to control how the vparams are sampled.
      //-------------------------------------------------------------
      int param_sampling = runOptions->getValue<int> ("param_sampling");
      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");
      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");
      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");
      //-------------------------------------------------------------
      // Priors on the IC and N_pivot ranges
      //-------------------------------------------------------------
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");
      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");
      //-------------------------------------------------------------
      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      //-------------------------------------------------------------
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");
      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");



      gambit_inflation_observables observs;

      //-------------------------------------------------------------
      // The function below calls the MultiModeCode backend
      //  for a given choice of inflationary model,
      //  which calculates the observables.
      //-------------------------------------------------------------
      BEreq::multimodecode_gambit_driver(&observs,
                                         num_inflaton,
                                         potential_choice,
                                         slowroll_infl_end,
                                         instreheat,
                                         vparam_rows,
                                         use_deltaN_SR,
                                         evaluate_modes,
                                         use_horiz_cross_approx,
                                         get_runningofrunning,
                                         ic_sampling,
                                         energy_scale,
                                         numb_samples,
                                         save_iso_N,
                                         N_iso_ref,
                                         param_sampling,
                                         byVal(&vp_prior_min[0]),
                                         byVal(&vp_prior_max[0]),
                                         varying_N_pivot,
                                         use_first_priorval,
                                         byVal(&phi_init0[0]),
                                         byVal(&dphi_init0[0]),
                                         byVal(&vparams[0]),
                                         N_pivot,
                                         k_pivot,
                                         dlnk,
                                         calc_full_pk,
                                         steps,
                                         kmin,
                                         byVal(&phi0_priors_min[0]),
                                         byVal(&phi0_priors_max[0]),
                                         byVal(&dphi0_priors_min[0]),
                                         byVal(&dphi0_priors_max[0]),
                                         N_pivot_prior_min,
                                         N_pivot_prior_max);

      if (calc_full_pk == 0)
      {
        if (silence == 1){
          std::cout << "A_s =" << observs.As << std::endl;
          // std::cout << "A_iso " << observs.A_iso << std::endl;
          // std::cout << "A_pnad " << observs.A_pnad << std::endl;
          // std::cout << "A_end " << observs.A_ent << std::endl;
          // std::cout << "A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
          std::cout << "n_s =" << observs.ns << std::endl;
          // std::cout << "n_t " << observs.nt << std::endl;
          // std::cout << "n_iso " << observs.n_iso << std::endl;
          // std::cout << "n_pnad " << observs.n_pnad << std::endl;
          // std::cout << "n_ent " << observs.n_ent << std::endl;
          std::cout << "r =" << observs.r << std::endl;
          std::cout << "alpha_s =" << observs.alpha_s << std::endl;
          // std::cout << "runofrun " << observs.runofrun << std::endl;
          // std::cout << "f_NL " << observs.f_NL << std::endl;
          // std::cout << "tau_NL " << observs.tau_NL << std::endl;
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("A_s",observs.As);
        cosmo.input.addEntry("n_s",observs.ns);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
        cosmo.input.addEntry("r",observs.r);
      }
      else
      {
        /* gambit_Pk type is set if full_spectra is asked. */
        cosmo.input.addEntry("P_k_ini type","gambit_Pk");

        cosmo.Pk_S.resize(steps+1, 0.);
        cosmo.Pk_T.resize(steps+1, 0.);
        cosmo.k_ar.resize(steps+1, 0.);

        for (int ii=0; ii < steps; ii++)
        {
          cosmo.k_ar.at(ii) = observs.k_array[ii];
          cosmo.Pk_S.at(ii) = observs.pks_array[ii];
          cosmo.Pk_T.at(ii) = observs.pkt_array[ii];
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
      }

      YAML::Node class_dict;
      if (runOptions->hasKey("class_dict"))
      {
        class_dict = runOptions->getValue<YAML::Node>("class_dict");
        for (auto it=class_dict.begin(); it != class_dict.end(); it++)
        {
          std::string name = it->first.as<std::string>();
          std::string value = it->second.as<std::string>();
          cosmo.input.addEntry(name,value);
        }
      }
    }

    void class_set_parameter_inf_1naturalInf_LCDMt(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_inf_SR1quad_LCDMt" << std::endl;
      using namespace Pipes::class_set_parameter_inf_1naturalInf_LCDMt;

      int l_max=cosmo.lmax;

      cosmo.input.clear();

      cosmo.input.addEntry("output","tCl pCl lCl");
      cosmo.input.addEntry("l_max_scalars",l_max);
      cosmo.input.addEntry("modes","s,t");
      cosmo.input.addEntry("lensing","yes");

      cosmo.input.addEntry("T_cmb",*Dep::T_cmb);
      // Parameters to be passed to the potential                       //
      //---------------------------------------------------------------//
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");
      //---------------------------------------------------------------//
      // Ideally what's between this line and the line (*) indicated    //
      // below is a general enough example of needed to set a new       //
      // inflationary model with MultiModeCode, i.e. the model          //
      // parameters that's to be sampled.                               //
      //---------------------------------------------------------------//
      //--------------- Sampling \lambda ------------------------------//
      //---------------------------------------------------------------//
      vparams[0] = *Param["lambda"];
      //---------------------------------------------------------------//
      //--------------- Sampling f       ------------------------------//
      //---------------------------------------------------------------//
      vparams[1] = *Param["faxion"];
      //---------------------------------------------------------------//
      //--------------- Sampling N_pivot ------------------------------//
      //---------------------------------------------------------------//
      double N_pivot = *Param["N_pivot"];
      //---------------------------------------------------------------//
      //------------------------------(*)------------------------------//



      //-------------------------------------------------------------
      //  Below we set settings for inflation solver MultiModecode.
      //-------------------------------------------------------------
      int silence = runOptions->getValue<int> ("is_not_silent");
      //-------------------------------------------------------------
      // Initialization parameters controlling main characteristics.
      //-------------------------------------------------------------
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");
      //-------------------------------------------------------------
      // Control the output of analytic approximations for comparison.
      //-------------------------------------------------------------
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");
      //-------------------------------------------------------------
      // Parameters to control how the ICs are sampled.
      //-------------------------------------------------------------
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");
      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this
      //-------------------------------------------------------------
      // Parameters to control how the vparams are sampled.
      //-------------------------------------------------------------
      int param_sampling = runOptions->getValue<int> ("param_sampling");
      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");
      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");
      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");
      //-------------------------------------------------------------
      // Priors on the IC and N_pivot ranges
      //-------------------------------------------------------------
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");
      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");
      //-------------------------------------------------------------
      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      //-------------------------------------------------------------
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");
      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");



      gambit_inflation_observables observs;

      //-------------------------------------------------------------
      // The function below calls the MultiModeCode backend
      //  for a given choice of inflationary model,
      //  which calculates the observables.
      //-------------------------------------------------------------
      BEreq::multimodecode_gambit_driver(&observs,
                                         num_inflaton,
                                         potential_choice,
                                         slowroll_infl_end,
                                         instreheat,
                                         vparam_rows,
                                         use_deltaN_SR,
                                         evaluate_modes,
                                         use_horiz_cross_approx,
                                         get_runningofrunning,
                                         ic_sampling,
                                         energy_scale,
                                         numb_samples,
                                         save_iso_N,
                                         N_iso_ref,
                                         param_sampling,
                                         byVal(&vp_prior_min[0]),
                                         byVal(&vp_prior_max[0]),
                                         varying_N_pivot,
                                         use_first_priorval,
                                         byVal(&phi_init0[0]),
                                         byVal(&dphi_init0[0]),
                                         byVal(&vparams[0]),
                                         N_pivot,
                                         k_pivot,
                                         dlnk,
                                         calc_full_pk,
                                         steps,
                                         kmin,
                                         byVal(&phi0_priors_min[0]),
                                         byVal(&phi0_priors_max[0]),
                                         byVal(&dphi0_priors_min[0]),
                                         byVal(&dphi0_priors_max[0]),
                                         N_pivot_prior_min,
                                         N_pivot_prior_max);

      if (calc_full_pk == 0)
      {
        if (silence == 1){
          std::cout << "A_s =" << observs.As << std::endl;
          // std::cout << "A_iso " << observs.A_iso << std::endl;
          // std::cout << "A_pnad " << observs.A_pnad << std::endl;
          // std::cout << "A_end " << observs.A_ent << std::endl;
          // std::cout << "A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
          std::cout << "n_s =" << observs.ns << std::endl;
          // std::cout << "n_t " << observs.nt << std::endl;
          // std::cout << "n_iso " << observs.n_iso << std::endl;
          // std::cout << "n_pnad " << observs.n_pnad << std::endl;
          // std::cout << "n_ent " << observs.n_ent << std::endl;
          std::cout << "r =" << observs.r << std::endl;
          std::cout << "alpha_s =" << observs.alpha_s << std::endl;
          // std::cout << "runofrun " << observs.runofrun << std::endl;
          // std::cout << "f_NL " << observs.f_NL << std::endl;
          // std::cout << "tau_NL " << observs.tau_NL << std::endl;
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("A_s",observs.As);
        cosmo.input.addEntry("n_s",observs.ns);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
        cosmo.input.addEntry("r",observs.r);
      }
      else
      {
        /* gambit_Pk type is set if full_spectra is asked. */
        cosmo.input.addEntry("P_k_ini type","gambit_Pk");

        cosmo.Pk_S.resize(steps+1, 0.);
        cosmo.Pk_T.resize(steps+1, 0.);
        cosmo.k_ar.resize(steps+1, 0.);

        for (int ii=0; ii < steps; ii++)
        {
          cosmo.k_ar.at(ii) = observs.k_array[ii];
          cosmo.Pk_S.at(ii) = observs.pks_array[ii];
          cosmo.Pk_T.at(ii) = observs.pkt_array[ii];
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
      }

      YAML::Node class_dict;
      if (runOptions->hasKey("class_dict"))
      {
        class_dict = runOptions->getValue<YAML::Node>("class_dict");
        for (auto it=class_dict.begin(); it != class_dict.end(); it++)
        {
          std::string name = it->first.as<std::string>();
          std::string value = it->second.as<std::string>();
          cosmo.input.addEntry(name,value);
        }
      }
    }

    void class_set_parameter_inf_1hilltopInf_LCDMt(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_inf_SR1quad_LCDMt" << std::endl;
      using namespace Pipes::class_set_parameter_inf_1hilltopInf_LCDMt;

      int l_max=cosmo.lmax;

      cosmo.input.clear();

      cosmo.input.addEntry("output","tCl pCl lCl");
      cosmo.input.addEntry("l_max_scalars",l_max);
      cosmo.input.addEntry("modes","s,t");
      cosmo.input.addEntry("lensing","yes");

      cosmo.input.addEntry("T_cmb",*Dep::T_cmb);
      // Parameters to be passed to the potential                       //
      //---------------------------------------------------------------//
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");
      //---------------------------------------------------------------//
      // Ideally what's between this line and the line (*) indicated    //
      // below is a general enough example of needed to set a new       //
      // inflationary model with MultiModeCode, i.e. the model          //
      // parameters that's to be sampled.                               //
      //---------------------------------------------------------------//
      //--------------- Sampling \lambda ------------------------------//
      //---------------------------------------------------------------//
      vparams[0] = *Param["lambda"];
      //---------------------------------------------------------------//
      //--------------- Sampling mu       ------------------------------//
      //---------------------------------------------------------------//
      vparams[1] = *Param["mu"];
      //---------------------------------------------------------------//
      //--------------- Sampling N_pivot ------------------------------//
      //---------------------------------------------------------------//
      double N_pivot = *Param["N_pivot"];
      //---------------------------------------------------------------//
      //------------------------------(*)------------------------------//
      std::cout << "lambda =" << vparams[0] << std::endl;
      std::cout << "mu =" << vparams[1] << std::endl;
      std::cout << "N_pivot =" << N_pivot << std::endl;


      //-------------------------------------------------------------
      //  Below we set settings for inflation solver MultiModecode.
      //-------------------------------------------------------------
      int silence = runOptions->getValue<int> ("is_not_silent");
      //-------------------------------------------------------------
      // Initialization parameters controlling main characteristics.
      //-------------------------------------------------------------
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");
      //-------------------------------------------------------------
      // Control the output of analytic approximations for comparison.
      //-------------------------------------------------------------
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");
      //-------------------------------------------------------------
      // Parameters to control how the ICs are sampled.
      //-------------------------------------------------------------
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");
      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this
      //-------------------------------------------------------------
      // Parameters to control how the vparams are sampled.
      //-------------------------------------------------------------
      int param_sampling = runOptions->getValue<int> ("param_sampling");
      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");
      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");
      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");
      //-------------------------------------------------------------
      // Priors on the IC and N_pivot ranges
      //-------------------------------------------------------------
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");
      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");
      //-------------------------------------------------------------
      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      //-------------------------------------------------------------
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");
      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");



      gambit_inflation_observables observs;

      //-------------------------------------------------------------
      // The function below calls the MultiModeCode backend
      //  for a given choice of inflationary model,
      //  which calculates the observables.
      //-------------------------------------------------------------
      BEreq::multimodecode_gambit_driver(&observs,
                                         num_inflaton,
                                         potential_choice,
                                         slowroll_infl_end,
                                         instreheat,
                                         vparam_rows,
                                         use_deltaN_SR,
                                         evaluate_modes,
                                         use_horiz_cross_approx,
                                         get_runningofrunning,
                                         ic_sampling,
                                         energy_scale,
                                         numb_samples,
                                         save_iso_N,
                                         N_iso_ref,
                                         param_sampling,
                                         byVal(&vp_prior_min[0]),
                                         byVal(&vp_prior_max[0]),
                                         varying_N_pivot,
                                         use_first_priorval,
                                         byVal(&phi_init0[0]),
                                         byVal(&dphi_init0[0]),
                                         byVal(&vparams[0]),
                                         N_pivot,
                                         k_pivot,
                                         dlnk,
                                         calc_full_pk,
                                         steps,
                                         kmin,
                                         byVal(&phi0_priors_min[0]),
                                         byVal(&phi0_priors_max[0]),
                                         byVal(&dphi0_priors_min[0]),
                                         byVal(&dphi0_priors_max[0]),
                                         N_pivot_prior_min,
                                         N_pivot_prior_max);

      if (calc_full_pk == 0)
      {
        if (silence == 1){
          std::cout << "A_s =" << observs.As << std::endl;
          // std::cout << "A_iso " << observs.A_iso << std::endl;
          // std::cout << "A_pnad " << observs.A_pnad << std::endl;
          // std::cout << "A_end " << observs.A_ent << std::endl;
          // std::cout << "A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
          std::cout << "n_s =" << observs.ns << std::endl;
          // std::cout << "n_t " << observs.nt << std::endl;
          // std::cout << "n_iso " << observs.n_iso << std::endl;
          // std::cout << "n_pnad " << observs.n_pnad << std::endl;
          // std::cout << "n_ent " << observs.n_ent << std::endl;
          std::cout << "r =" << observs.r << std::endl;
          std::cout << "alpha_s =" << observs.alpha_s << std::endl;
          // std::cout << "runofrun " << observs.runofrun << std::endl;
          // std::cout << "f_NL " << observs.f_NL << std::endl;
          // std::cout << "tau_NL " << observs.tau_NL << std::endl;
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("A_s",observs.As);
        cosmo.input.addEntry("n_s",observs.ns);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
        cosmo.input.addEntry("r",observs.r);
      }
      else
      {
        /* gambit_Pk type is set if full_spectra is asked. */
        cosmo.input.addEntry("P_k_ini type","gambit_Pk");

        cosmo.Pk_S.resize(steps+1, 0.);
        cosmo.Pk_T.resize(steps+1, 0.);
        cosmo.k_ar.resize(steps+1, 0.);

        for (int ii=0; ii < steps; ii++)
        {
          cosmo.k_ar.at(ii) = observs.k_array[ii];
          cosmo.Pk_S.at(ii) = observs.pks_array[ii];
          cosmo.Pk_T.at(ii) = observs.pkt_array[ii];
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
      }

      YAML::Node class_dict;
      if (runOptions->hasKey("class_dict"))
      {
        class_dict = runOptions->getValue<YAML::Node>("class_dict");
        for (auto it=class_dict.begin(); it != class_dict.end(); it++)
        {
          std::string name = it->first.as<std::string>();
          std::string value = it->second.as<std::string>();
          cosmo.input.addEntry(name,value);
        }
      }
    }

    void class_set_parameter_inf_smashInf_LCDMt(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_inf_smashInf_LCDMt" << std::endl;
      using namespace Pipes::class_set_parameter_inf_smashInf_LCDMt;

      //---------------------------------------------------------------//
      // Below, SMASH inflationary potential is calculated and
      // the slow-roll apporximation is used to get the observables
      // A_s, n_s and r from model parameters and the inital value
      // of the inflaton potential at which infator starts its
      // slow-roll with zero-velocity.
      //---------------------------------------------------------------//

      int calc_tech = runOptions->getValue<int> ("calc_technique");
      int silence = runOptions->getValue<int> ("is_not_silent");

      int l_max=cosmo.lmax;

      cosmo.input.clear();

      cosmo.input.addEntry("output","tCl pCl lCl");
      cosmo.input.addEntry("l_max_scalars",l_max);
      cosmo.input.addEntry("modes","s,t");
      cosmo.input.addEntry("lensing","yes");

      cosmo.input.addEntry("T_cmb",*Dep::T_cmb);


      if (calc_tech==0){
        //-------------------------------------------------------------
        // Having GAMBIT-defined functions solve for the slow-roll parameters which is then given to CLASS.
        //-------------------------------------------------------------


        //-------------------------------------------------------------
        // Parameters to be passed to the potential
        //-------------------------------------------------------------
        std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

        vparams[0] = *Param["log10_xi"];
        vparams[1] = *Param["log10_beta"];
        vparams[2] = *Param["log10_lambda"];

        if (silence==1){
          std::cout << "log10[xi] = " << vparams[0] << std::endl;
          std::cout << "log10[beta] = " << vparams[1] << std::endl;
          std::cout << "log10[lambda] = " << vparams[2] << std::endl;
        }
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
        double m = tempbw*10.0, m_expected = M_PI; // Correct for.
        double a = 0.0, b = 100.0; // Correct for.
        gsl_function F;

        struct my_f_params params = { vparams[0], vparams[1], vparams[3] , N_pivot};

        F.function = &fn1;
        F.params = &params;

        T = gsl_min_fminimizer_brent;
        s = gsl_min_fminimizer_alloc (T);
        gsl_min_fminimizer_set (s, &F, m, a, b);

        /*
         printf ("using %s method\n",gsl_min_fminimizer_name (s));
         printf ("%5s [%9s, %9s] %9s %10s %9s\n","iter", "lower", "upper", "min","err", "err(est)");
         printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",iter, a, b, m, m - m_expected, b - a);
         */

        do
        {
          iter++;
          status = gsl_min_fminimizer_iterate (s);

          m = gsl_min_fminimizer_x_minimum (s);
          a = gsl_min_fminimizer_x_lower (s);
          b = gsl_min_fminimizer_x_upper (s);

          status = gsl_min_test_interval (a, b, 0.01, 0.0);

          /*
           if (status == GSL_SUCCESS) printf ("Converged:\n");
           printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
           iter, a, b,
           m, m - m_expected, b - a);
           */
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
        r_self  = r_SR(smash_eps,smash_eta);

        //-------------------------------------------------------------
        /* Have calculated the field value at N_pivot
         predicted by slow-roll approximation.
         Choosing a very close slightly higher
         value will be satisfactory for initial conditions.
         */
        //-------------------------------------------------------------

        if (silence == 1)
        {
          printf("n_s from smash = %e\n",ns_self);
          printf("A_s from smash = %e\n",As_self);
          printf("r from smash = %e\n",r_self);
        }

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("A_s",As_self);
        cosmo.input.addEntry("n_s",ns_self);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);
        cosmo.input.addEntry("r",r_self);

      }
      else if (calc_tech==1)
      {
        //-------------------------------------------------------------
        // Having CLASS primordial.h module calculate inflationary predictions of the SMASH potential.
        //-------------------------------------------------------------

        cosmo.input.addEntry("P_k_ini type","inflation_V");
        cosmo.input.addEntry("potential","smash_inflation");

        cosmo.input.addEntry("V_0",*Param["log10_xi"]);
        cosmo.input.addEntry("V_1",*Param["log10_beta"]);
        cosmo.input.addEntry("V_3",*Param["log10_lambda"]);

        cosmo.input.addEntry("V_2",-100.); //Hard coding \chi inflaton potential boundaries.
//        cosmo.input.addEntry("N_star",*Param["N_pivot"]); //Hard coding \chi inflaton potential boundaries.

        cosmo.input.addEntry("omega_b",*Param["omega_b"]);
        cosmo.input.addEntry("omega_cdm",*Param["omega_cdm"]);
        cosmo.input.addEntry("H0",*Param["H0"]);
        cosmo.input.addEntry("tau_reio",*Param["tau_reio"]);

      }

      YAML::Node class_dict;
      if (runOptions->hasKey("class_dict"))
      {
        class_dict = runOptions->getValue<YAML::Node>("class_dict");
        for (auto it=class_dict.begin(); it != class_dict.end(); it++)
        {
          std::string name = it->first.as<std::string>();
          std::string value = it->second.as<std::string>();
          cosmo.input.addEntry(name,value);
        }
      }
    }

    void class_get_spectra_func(Class_container& cosmo)
    {
//      std::cout << "Last seen alive in: class_get_spectra_func" << std::endl;
      using namespace Pipes::class_get_spectra_func;
      cosmo = BEreq::get_ptr_to_class();

      cosmo.Cl_TT = BEreq::class_get_cl("tt");
      cosmo.Cl_TE = BEreq::class_get_cl("te");
      cosmo.Cl_EE = BEreq::class_get_cl("ee");
      cosmo.Cl_BB = BEreq::class_get_cl("bb");
      cosmo.Cl_PhiPhi = BEreq::class_get_cl("pp");
    }

    void function_Planck_high_TT_loglike(double& result)
    {
//      std::cout << "Last seen alive in: function_Planck_high_TT_loglike" << std::endl;
      using namespace Pipes::function_Planck_high_TT_loglike;

      double  cl_and_pars[2525];
      int idx_tt;
      Class_container cosmo = *Dep::class_get_spectra;

      //--------------------------------------------------------------------------
      //------addition of the Cl for TT, TE, EE and BB to Cl array----------------
      //--------------------------------------------------------------------------

      for(int ii = 0; ii < 2509 ; ii++)
      {
        idx_tt = ii;

        if (ii >= 2)
        {
          cl_and_pars[idx_tt] = cosmo.Cl_TT.at(ii);
        }
        else
        {
          cl_and_pars[idx_tt] = 0.;
        }
      }

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      cl_and_pars[2509] = *Param["A_cib_217"];
      cl_and_pars[2510] = *Param["cib_index"];
      cl_and_pars[2511] = *Param["xi_sz_cib"];
      cl_and_pars[2512] = *Param["A_sz"];
      cl_and_pars[2513] = *Param["ps_A_100_100"];
      cl_and_pars[2514] = *Param["ps_A_143_143"];
      cl_and_pars[2515] = *Param["ps_A_143_217"];
      cl_and_pars[2516] = *Param["ps_A_217_217"];
      cl_and_pars[2517] = *Param["ksz_norm"];
      cl_and_pars[2518] = *Param["gal545_A_100"];
      cl_and_pars[2519] = *Param["gal545_A_143"];
      cl_and_pars[2520] = *Param["gal545_A_143_217"];
      cl_and_pars[2521] = *Param["gal545_A_217"];
      cl_and_pars[2522] = *Param["calib_100T"];
      cl_and_pars[2523] = *Param["calib_217T"];
      cl_and_pars[2524] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_object* high_clikid;
      clik_error *_err;

      high_clikid = BEreq::return_high_TT();
      _err = BEreq::clik_initialize_error();
      result = BEreq::clik_compute_loglike(byVal(high_clikid),
               byVal(cl_and_pars),
               &_err);

//      std::cout << "Log likelihood (of high_TT) is : " << result << std::endl;
    }

    void function_Planck_high_TTTEEE_loglike(double& result)
    {
      //std::cout << "Last seen alive in: function_Planck_high_TTTEEE_loglike" << std::endl;
      using namespace Pipes::function_Planck_high_TTTEEE_loglike;

      double  cl_and_pars[7621];
      int idx_tt, idx_te, idx_ee;
      Class_container cosmo = *Dep::class_get_spectra;

      //--------------------------------------------------------------------------
      //------addition of the Cl for TT, TE and EE to Cl array--------------------
      //--------------------------------------------------------------------------
      for(int ii = 0; ii < 2509 ; ii++)
      {
        idx_tt = ii;
        idx_ee = ii + 2509;
        idx_te = ii + (2 * 2509);
        if (ii >= 2)
        {
          cl_and_pars[idx_tt] = cosmo.Cl_TT.at(ii);
          cl_and_pars[idx_ee] = cosmo.Cl_EE.at(ii);
          cl_and_pars[idx_te] = cosmo.Cl_TE.at(ii);
        }
        else
        {
          cl_and_pars[idx_tt] = 0.;
          cl_and_pars[idx_ee] = 0.;
          cl_and_pars[idx_te] = 0.;
        }
      }

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      cl_and_pars[7527] = *Param["A_cib_217"];
      cl_and_pars[7528] = *Param["cib_index"];
      cl_and_pars[7529] = *Param["xi_sz_cib"];
      cl_and_pars[7530] = *Param["A_sz"];
      cl_and_pars[7531] = *Param["ps_A_100_100"];
      cl_and_pars[7532] = *Param["ps_A_143_143"];
      cl_and_pars[7533] = *Param["ps_A_143_217"];
      cl_and_pars[7534] = *Param["ps_A_217_217"];
      cl_and_pars[7535] = *Param["ksz_norm"];
      cl_and_pars[7536] = *Param["gal545_A_100"];
      cl_and_pars[7537] = *Param["gal545_A_143"];
      cl_and_pars[7538] = *Param["gal545_A_143_217"];
      cl_and_pars[7539] = *Param["gal545_A_217"];
      cl_and_pars[7540] = *Param["galf_EE_A_100"];
      cl_and_pars[7541] = *Param["galf_EE_A_100_143"];
      cl_and_pars[7542] = *Param["galf_EE_A_100_217"];
      cl_and_pars[7543] = *Param["galf_EE_A_143"];
      cl_and_pars[7544] = *Param["galf_EE_A_143_217"];
      cl_and_pars[7545] = *Param["galf_EE_A_217"];
      cl_and_pars[7546] = *Param["galf_EE_index"];
      cl_and_pars[7547] = *Param["galf_TE_A_100"];
      cl_and_pars[7548] = *Param["galf_TE_A_100_143"];
      cl_and_pars[7549] = *Param["galf_TE_A_100_217"];
      cl_and_pars[7550] = *Param["galf_TE_A_143"];
      cl_and_pars[7551] = *Param["galf_TE_A_143_217"];
      cl_and_pars[7552] = *Param["galf_TE_A_217"];
      cl_and_pars[7553] = *Param["galf_TE_index"];
      // set beam-leakage to zero (60 nusissance parameter)
      for (int i = 0; i < 60; i++) cl_and_pars[(i+7554)] = 0.;
      cl_and_pars[7614] = *Param["calib_100T"];
      cl_and_pars[7615] = *Param["calib_217T"];
      cl_and_pars[7616] = *Param["calib_100P"];
      cl_and_pars[7617] = *Param["calib_143P"];
      cl_and_pars[7618] = *Param["calib_217P"];
      cl_and_pars[7619] = *Param["A_pol"];
      cl_and_pars[7620] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_object* high_clikid;
      clik_error *_err;

      high_clikid = BEreq::return_high_TTTEEE();
      _err = BEreq::clik_initialize_error();

      result = BEreq::clik_compute_loglike(byVal(high_clikid),
               byVal(cl_and_pars),
               &_err);

//      std::cout << "Log likelihood (of high_TTTEEE) is : " << result << std::endl;
    }

    void function_Planck_high_TT_lite_loglike(double& result)
    {
      //std::cout << "Last seen alive in: function_Planck_high_TT_lite_loglike" << std::endl;
      using namespace Pipes::function_Planck_high_TT_lite_loglike;

      double cl_and_pars[2510];
      int idx_tt;
      Class_container cosmo = *Dep::class_get_spectra;

      //--------------------------------------------------------------------------
      //------addition of the Cl for TT, TE, EE and BB to Cl array----------------
      //--------------------------------------------------------------------------
      for(int ii = 0; ii < 2509 ; ii++)
      {
        idx_tt = ii;
        if (ii >= 2)
        {
          cl_and_pars[idx_tt] = cosmo.Cl_TT.at(ii);
        }
        else
        {
          cl_and_pars[idx_tt] = 0.;
        }
      }

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      cl_and_pars[2509] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_object* high_clikid;
      clik_error *_err;

      high_clikid = BEreq::return_high_TT_lite();
      _err = BEreq::clik_initialize_error();
      result = BEreq::clik_compute_loglike(byVal(high_clikid),
               byVal(cl_and_pars),
               &_err);

//      std::cout << "Log likelihood (of high_TT_lite) is : " << result << std::endl;
    }

    void function_Planck_high_TTTEEE_lite_loglike(double& result)
    {
      //std::cout << "Last seen alive in: function_Planck_high_TTTEEE_lite_loglike" << std::endl;
      using namespace Pipes::function_Planck_high_TTTEEE_lite_loglike;

      double  cl_and_pars[7528];
      int idx_tt, idx_te, idx_ee;
      Class_container cosmo = *Dep::class_get_spectra;

      //--------------------------------------------------------------------------
      //------addition of the Cl for TT, TE and EE to Cl array--------------------
      //--------------------------------------------------------------------------
      for(int ii = 0; ii < 2509 ; ii++)
      {
        idx_tt = ii;
        idx_ee = ii + 2509;
        idx_te = ii + (2 * 2509);
        if (ii >= 2)
        {
          cl_and_pars[idx_tt] = cosmo.Cl_TT.at(ii);
          cl_and_pars[idx_ee] = cosmo.Cl_EE.at(ii);
          cl_and_pars[idx_te] = cosmo.Cl_TE.at(ii);
        }
        else
        {
          cl_and_pars[idx_tt] = 0.;
          cl_and_pars[idx_ee] = 0.;
          cl_and_pars[idx_te] = 0.;
        }
      }

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      cl_and_pars[7527] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_object* high_clikid;
      clik_error *_err;

      high_clikid = BEreq::return_high_TTTEEE_lite();
      _err = BEreq::clik_initialize_error();

      result = BEreq::clik_compute_loglike(byVal(high_clikid),
               byVal(cl_and_pars),
               &_err);

      //std::cout << "Log likelihood (of high_TTTEEE_lite) is : " << result << std::endl;
    }

    void function_Planck_lensing_loglike(double& result)
    {
      //std::cout << "Last seen alive in: function_Planck_lensing_loglike" << std::endl;
      using namespace Pipes::function_Planck_lensing_loglike;

      double  cl_and_pars[8197];
      int idx_pp, idx_tt, idx_te, idx_ee;
      Class_container cosmo = *Dep::class_get_spectra;

      //--------------------------------------------------------------------------
      //------addition of the Cl for PhiPhi,  TT, TE and EE to Cl array-----------
      //--------------------------------------------------------------------------
      for(int ii = 0; ii < 2049 ; ii++)
      {
        idx_pp = ii;
        idx_tt = ii + 2049;
        idx_ee = ii + (2 * 2049);
        idx_te = ii + (3 * 2049);
        if (ii >= 2)
        {
          cl_and_pars[idx_pp] = cosmo.Cl_PhiPhi.at(ii);
          cl_and_pars[idx_tt] = cosmo.Cl_TT.at(ii);
          cl_and_pars[idx_ee] = cosmo.Cl_EE.at(ii);
          cl_and_pars[idx_te] = cosmo.Cl_TE.at(ii);
        }
        else
        {
          cl_and_pars[idx_pp] = 0.;
          cl_and_pars[idx_tt] = 0.;
          cl_and_pars[idx_ee] = 0.;
          cl_and_pars[idx_te] = 0.;
        }
      }

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      cl_and_pars[8196] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_lensing_object* lensing_clikid;
      clik_error *_err;

      lensing_clikid = BEreq::return_lensing();
      _err = BEreq::clik_initialize_error();
      result = BEreq::clik_lensing_compute_loglike(byVal(lensing_clikid),
               byVal(cl_and_pars),
               &_err);

//      std::cout << "Log likelihood (of lensing) is : " << result << std::endl;
    }

    void function_Planck_lowp_TT_loglike(double& result)
    {
      //std::cout << "Last seen alive in: function_Planck_lowp_TT_loglike" << std::endl;
      using namespace Pipes::function_Planck_lowp_TT_loglike;

      double  cl_and_pars[121];
      int idx_tt, idx_te, idx_ee, idx_bb;
      Class_container cosmo = *Dep::class_get_spectra;

      //--------------------------------------------------------------------------
      //------addition of the Cl for TT, TE, EE and BB to Cl array----------------
      //--------------------------------------------------------------------------
      for(int ii = 0; ii < 30 ; ii++)
      {
        idx_tt = ii;
        idx_ee = ii + 30;
        idx_bb = ii + (2 * 30);
        idx_te = ii + (3 * 30);
        if (ii >= 2)
        {
          cl_and_pars[idx_tt] = cosmo.Cl_TT.at(ii);
          cl_and_pars[idx_ee] = cosmo.Cl_EE.at(ii);
          cl_and_pars[idx_bb] = cosmo.Cl_BB.at(ii);
          cl_and_pars[idx_te] = cosmo.Cl_TE.at(ii);
        }
        else
        {
          cl_and_pars[idx_tt] = 0.;
          cl_and_pars[idx_te] = 0.;
          cl_and_pars[idx_bb] = 0.;
          cl_and_pars[idx_te] = 0.;
        }
      }

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      cl_and_pars[120] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_object* lowl_clikid;
      clik_error *_err;

      lowl_clikid = BEreq::return_lowp_TT();
      _err = BEreq::clik_initialize_error();

      result = BEreq::clik_compute_loglike(byVal(lowl_clikid),
               byVal(cl_and_pars),
               &_err);

//      std::cout << "Log likelihood (of lowp_TT) is : " << result << std::endl;
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
      double Ya0 = *Param["Ya0"];
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
            err << "Failed to solve differential equation to compute dNeffCMB and etaBB for Cosmo_ALP model at iteration step "<< ii <<". Invalidating point";
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
        err << "Computation of dNeffCMB and etaBBN for Cosmo_ALP model did not converge after n = "<< ii <<" iterations. You can increase the maximum number of iterations with the run Option 'max_iter'. Invalidating point.";
        invalid_point().raise(err.str());
      }

      result["dNeff"] = dNeff;
      result["Neff_ratio"] = Neff_ALP/Neff_SM;
      result["eta_ratio"] = eta_ratio;
      logger() << "CosmoALP model: calculated Neff @BBN to be " << result["dNeff"] <<", and etaBB(ALP)/etaBBN(SM) = " << result["eta_ratio"] << ". Calculation converged after "<< ii <<" iterations." << EOM;

    }

/// -----------
/// Capabilities for general cosmological quantities
/// -----------


    void set_T_cmb(double &result)
    {
      using namespace Pipes::set_T_cmb;

      // use COBE/FIRAS (0911.1955) mes. if not specified otherwise
      result = runOptions->getValueOrDef<double>(2.7255,"T_cmb");
    }

    void set_T_ncdm_SM(double &result)
    {
      using namespace Pipes::set_T_ncdm_SM;

      // set to 0.71611 in units of photon temperature, above the instantaneous decoupling value (4/11)^(1/3)
      // to recover Sum_i mNu_i/omega = 93.14 eV resulting from studies of active neutrino decoupling (hep-ph/0506164)
      result = 0.71611;
      // This standard values enters in many assumption entering class. Therefore changing this value in 
      // the yaml file is disabled at the moment. If you still want to modify it uncomment the line below and 
      // you can set is as a run option (as T_cmb).
      // result = runOptions->getValueOrDef<double>(0.71611,"T_ncdm");
    }

    void set_T_ncdm(double &result)
    {
      using namespace Pipes::set_T_ncdm;

      result = *Dep::T_ncdm_SM;
    }

    void set_T_ncdm_CosmoALP(double &result)
    {
      using namespace Pipes::set_T_ncdm_CosmoALP;

      map_str_dbl dNeff_CosmoALP = *Dep::external_dNeff_etaBBN;
      // SM value of 0.71611 (see set_T_ncdm) multiplied by (Neff_ALP/Neff_SM)^(1/4) = T_nu_ALP/T_nu_SM
      // to account for modified Neutrino temperature today (ALPs heat gamma w.r.t. to photons -> Neutrino today are cooler)
      result = 0.71611*pow(dNeff_CosmoALP["Neff_ratio"], 1./4.);

    }

    void calculate_eta0(double &result)
    {
      using namespace Pipes::calculate_eta0;

      double ngamma, nb;
      ngamma = 16*pi*zeta3*pow(*Dep::T_cmb*_kB_eV_over_K_/_hc_eVcm_,3); // photon number density today
      nb = *Param["omega_b"]*3*100*1e3*100*1e3/_Mpc_SI_/_Mpc_SI_/(8*pi*_GN_cgs_*_m_proton_g_); // baryon number density today

      result =  nb/ngamma;
      logger() << "Baryon to photon ratio (eta) today computed to be " << result << EOM;
    }

    void compute_Sigma8(double &result)
    {
      using namespace Pipes::compute_Sigma8;

      double sigma8 = BEreq::class_get_sigma8(0.);
      double Omega0_m = *Dep::Omega0_m;

      result = sigma8*pow(Omega0_m/0.3, 0.5);
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
      result = (4.*_sigmaB_SI_/_c_SI_*pow(*Dep::T_cmb,4.)) / (3.*_c_SI_*_c_SI_*1.e10*h*h/_Mpc_SI_/_Mpc_SI_/8./pi/_GN_SI_);
    }

    void compute_n0_g(double &result)
    {
      using namespace Pipes::compute_n0_g;

      result = 2./pi/pi*zeta3 *pow(*Dep::T_cmb*_kB_eV_over_K_,3.)/pow(_hP_eVs_*_c_SI_/2./pi,3)/100/100/100; // result per cm^3  
    }

    void compute_Omega0_ur(double &result)
    {
      using namespace Pipes::compute_Omega0_ur;

      // Capability class_Nur takes care of setting N_ur right for considered number of 
      // massive neutrinos & extra model dependent contributions
      result = (*Dep::class_Nur)*7./8.*pow(4./11.,4./3.)*(*Dep::Omega0_g); 
    }

    void compute_Omega0_ncdm(double &result)
    {
      using namespace Pipes::compute_Omega0_ncdm;

      double mNu_tot_eV; // sum of neutrino masses
      if (ModelInUse("StandardModel_SLHA2"))
      {
        // (PS) Heads up! The units in StandardModel_SLHA2 are GeV
        // Here we are using eV
        mNu_tot_eV = 1e9*(*Param["mNu1"] + (*Param["mNu2"]) + (*Param["mNu3"]) );
      }

      else
      {
        mNu_tot_eV = 0.06; // default to standard Planck case of 1 massive neutrino with m = 0.06 eV
      }

      double h = *Param["H0"]/100;
      result = mNu_tot_eV/(93.14*h*h);  // TODO: heads up: explicit assumption of T_ncdm = 0.71611 and T_cmb goes in here. Has to be generalised  
    }

    void compute_Omega0_r(double &result)
    {
      using namespace Pipes::compute_Omega0_r;

      result = *Dep::Omega0_g + (*Dep::Omega0_ur);
      //std::cout << "omega0_g g " << *Dep::Omega0_g << " omega_ur " << *Dep::Omega0_ur <<  std::endl; 
    }


    void calculate_etaCMB_SM(double &result)
    {
      using namespace Pipes::calculate_etaCMB_SM;

      result =  *Dep::eta0; // in SM the baryon to photon ratio does not change between today an CMB release
      logger() << "Baryon to photon ratio (eta) @CMB computed to be " << result << EOM;
    }

    void set_etaBBN(double &result)
    {
      using namespace Pipes::set_etaBBN;

      result =  *Param["eta_BBN"]; // in SM the baryon to photon ratio does not change between today an CMB release
      logger() << "Baryon to photon ratio (eta) @BBN set to be " << result << EOM;
    }

    void calculate_etaBBN_SM(double &result)
    {
      using namespace Pipes::calculate_etaBBN_SM;

      result =  *Dep::etaCMB; // in SM the baryon to photon ratio does not change between today an CMB release
      logger() << "Baryon to photon ratio (eta) @BBN computed to be " << result << EOM;
    }

    void calculate_etaBBN_ALP(double &result)
    {
      using namespace Pipes::calculate_etaBBN_ALP;

       map_str_dbl dNeff_etaBBN = *Dep::external_dNeff_etaBBN;
       result = (*Dep::etaCMB) * dNeff_etaBBN["eta_ratio"];
       logger() << "Baryon to photon ratio (eta) @BBN computed to be " << result << EOM;
    }

    void calculate_dNeffCMB_ALP(double &result)
    {
      using namespace Pipes::calculate_dNeffCMB_ALP;
     
      map_str_dbl dNeff_etaBBN = *Dep::external_dNeff_etaBBN;
      result = dNeff_etaBBN["dNeff"];
      logger() << "dNeff for ALP calculated to be " << result << EOM;

    }

     void set_class_Nur_LCDM_family(double &result)
    {
      using namespace Pipes::set_class_Nur_LCDM_family;

      map_str_dbl NuMasses_SM = *Dep::NuMasses_SM;
      result = *Param["dNeff"] + NuMasses_SM["N_ur_SMnu"];
      logger() << "N_ur for LCDM family calculated to be " << result << EOM;
    }


    void set_class_Nur_CosmoALP(double &result)
    {
      using namespace Pipes::set_class_Nur_CosmoALP;

      map_str_dbl dNeff_etaBBN = *Dep::external_dNeff_etaBBN;
      map_str_dbl NuMasses_SM = *Dep::NuMasses_SM;
      // N_ur = (NeffCMB_CosmoALP/NeffCMB_SM) * Nur_SM , to rescale the value for N_ur from non-instantaneous 
      // decoupling with the standard value of 3.046 to the according contribution for a modified Neutrino temperature
      // (CosmoALP model heats CMB photons w.r.t. to neutrinos)
      // TODO: refer to relevant section in paper!
      result = dNeff_etaBBN["Neff_ratio"]*NuMasses_SM["N_ur_SMnu"];
      logger() << "N_ur for ALP calculated to be " << result << EOM;
    }


/// -----------
/// BBN related functions (call AlterBBN, BBN abundances & Likelihood)
/// -----------

    void AlterBBN_Input_LCDM_dNeffCMB_dNeffBBN_etaBBN(map_str_dbl &result)
    {
      // add parameters of relicparam structure that should be set to non-default values
      // to the AlterBBN_input map. 
      // If you want to modify a parameter which has not been used in CosmoBit before simply
      // add it to the function 'fill_cosmomodel' in 'AlterBBN_<version>.cpp' and to the 
      // set of 'known' parameters 'known_relicparam_options'
      using namespace Pipes::AlterBBN_Input_LCDM_dNeffCMB_dNeffBBN_etaBBN;

      result["eta0"] = *Param["eta_BBN"];  // eta AFTER BBN (variable during)
      result["Nnu"]=3.046;                 // contribution from SM neutrinos
      result["dNnu"]=*Param["dNeff_BBN"];  // dNnu: within AlterBBN scenarios in which the sum Nnu+dNnu is the same are identical
      result["failsafe"] = runOptions->getValueOrDef<double>(3,"failsafe");
      result["err"] = runOptions->getValueOrDef<double>(3,"err");

      std::cout << "Set AlterBBN with parameters eta = " << result["eta0"] << ", Nnu = " << result["Nnu"] << ", dNeffBBN = " << result["dNnu"] << std::endl;
      logger() << "Set AlterBBN with parameters eta = " << result["eta0"] << ", Nnu = " << result["Nnu"] << ", dNeffBBN = " << result["dNnu"] << EOM;
      logger() << "     and error params: failsafe = " << result["failsafe"] << ", err = " << result["err"] << EOM;
    }

    void compute_BBN_abundances(CosmoBit::BBN_container &result)
    {
      using namespace Pipes::compute_BBN_abundances;

      int NNUC = BEreq::get_NNUC();  // global variable of AlterBBN (# computed element abundances)
      result.init_arr(NNUC);         // init arrays in BBN_container with right length
      double ratioH [NNUC+1], cov_ratioH [NNUC+1][NNUC+1];

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
            err << "The correlation matirx is not a square matrix";
            CosmoBit_error().raise(LOCAL_INFO, err.str());
          }
        }

        // Check if the entries in the correlation matrix are reasonable
        for (int ie=0; ie<nelements; ie++)
        {
          //Check if the diagonal entries are equal to 1.
          if (std::abs(tmp_corr.at(ie).at(ie) - 1.) > 1e-6)
          {
            std::ostringstream err;
            err << "Not all diagonal elements of the correlation matirx are 1.";
            CosmoBit_error().raise(LOCAL_INFO, err.str());
          }
          for (int je=0; je<=ie; je++)
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
        std::map<std::string, int> abund_map = result.get_map();
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

      map_str_dbl AlterBBN_input = *Dep::AlterBBN_setInput; // fill AlterBBN_input map with the parameters for the model in consideration
      int nucl_err = BEreq::call_nucl_err(AlterBBN_input, byVal(ratioH), byVal(cov_ratioH[0]));

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
            err_ratio.at(ie) = sqrt(cov_ratioH[ie][ie]);
          }
        }
      }

      // fill abundances and covariance matrix of BBN_container with results from AlterBBN
      if (nucl_err)
      {
        for(int ie=1;ie<=NNUC;ie++)
        {
          result.BBN_abund.at(ie) = ratioH[ie];
          for(int je=1;je<=NNUC;je++)
          {
            if (use_fudged_correlations)
              result.BBN_covmat.at(ie).at(je) = corr.at(ie).at(je) * err_ratio.at(ie) * err_ratio.at(je);
            else
              result.BBN_covmat.at(ie).at(je) = cov_ratioH[ie][je];
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
      std::map<std::string, int> abund_map = BBN_res.get_map();

      result.clear();
      result.push_back(BBN_res.BBN_abund.at(abund_map["Yp"]));
      result.push_back(sqrt(BBN_res.BBN_covmat.at(abund_map["Yp"]).at(abund_map["Yp"])));

      logger() << "Helium Abundance from AlterBBN: " << result.at(0) << " +/- " << result.at(1) << EOM;
    }

    void get_Helium3_abundance(std::vector<double> &result)
    {
      using namespace Pipes::get_Helium3_abundance;

      CosmoBit::BBN_container BBN_res = *Dep::BBN_abundances;
      std::map<std::string, int> abund_map = BBN_res.get_map();

      result.clear();
      result.push_back(BBN_res.BBN_abund.at(abund_map["He3"]));
      result.push_back(sqrt(BBN_res.BBN_covmat.at(abund_map["He3"]).at(abund_map["He3"])));

      logger() << "Helium 3 Abundance from AlterBBN: " << result.at(0) << " +/- " << result.at(1) << EOM;
    }

    void get_Deuterium_abundance(std::vector<double> &result)
    {
      using namespace Pipes::get_Deuterium_abundance;

      CosmoBit::BBN_container BBN_res = *Dep::BBN_abundances;
      std::map<std::string, int> abund_map = BBN_res.get_map();

      result.clear();
      result.push_back(BBN_res.BBN_abund.at(abund_map["D"]));
      result.push_back(sqrt(BBN_res.BBN_covmat.at(abund_map["D"]).at(abund_map["D"])));

      logger() << "Deuterium Abundance from AlterBBN: " << result.at(0) << " +/- " << result.at(1) << EOM;
    }

    void get_Lithium7_abundance(std::vector<double> &result)
    {
      using namespace Pipes::get_Lithium7_abundance;

      CosmoBit::BBN_container BBN_res = *Dep::BBN_abundances;
      std::map<std::string, int> abund_map = BBN_res.get_map();

      result.clear();
      result.push_back(BBN_res.BBN_abund.at(abund_map["Li7"]));
      result.push_back(sqrt(BBN_res.BBN_covmat.at(abund_map["Li7"]).at(abund_map["Li7"])));

      logger() << "Lithium Abundance from AlterBBN: " << result.at(0) << " +/- " << result.at(1) << EOM;
    }

    void get_Beryllium7_abundance(std::vector<double> &result)
    {
      using namespace Pipes::get_Beryllium7_abundance;

      CosmoBit::BBN_container BBN_res = *Dep::BBN_abundances;
      std::map<std::string, int> abund_map = BBN_res.get_map();

      result.clear();
      result.push_back(BBN_res.BBN_abund.at(abund_map["Be7"]));
      result.push_back(sqrt(BBN_res.BBN_covmat.at(abund_map["Be7"]).at(abund_map["Be7"])));

      logger() << "Beryllium Abundance from AlterBBN: " << result.at(0) << " +/- " << result.at(1) << EOM;
    }

    void compute_BBN_LogLike(double &result)
    {
      using namespace Pipes::compute_BBN_LogLike;

      double chi2 = 0;
      int ii = 0;
      int ie,je,s;

      CosmoBit::BBN_container BBN_res = *Dep::BBN_abundances; // fill BBN_container with abundance results from AlterBBN
      std::map<std::string, int> abund_map = BBN_res.get_map();

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
        prediction[ii]= BBN_res.BBN_abund.at(abund_map[key]);
        ii++;
      }

      // fill covmat
      for(ie=0;ie<nobs;ie++) {for(je=0;je<nobs;je++) gsl_matrix_set(cov, ie, je,pow(sigmaobs[ie],2.)*(ie==je)+BBN_res.BBN_covmat.at(translate[ie]).at(translate[je]));}

      // Compute the LU decomposition and inverse of cov mat
      gsl_linalg_LU_decomp(cov,p,&s);
      gsl_linalg_LU_invert(cov,p,invcov);

      // Compute the determinant of the inverse of the covmat
      double det_cov = gsl_linalg_LU_det(cov,s);

      // compute chi2
      for(ie=0;ie<nobs;ie++) for(je=0;je<nobs;je++) chi2+=(prediction[ie]-observed[ie])*gsl_matrix_get(invcov,ie,je)*(prediction[je]-observed[je]);
      result = -0.5*(chi2 + log(pow(2*pi,nobs)*det_cov));

      gsl_matrix_free(cov);
      gsl_matrix_free(invcov);
      gsl_permutation_free(p);
    }

/// -----------
/// Background Likelihoods (SNe + H0 + BAO + Sigma8)
/// -----------

    void compute_Pantheon_LogLike(double &result)
    {
      using namespace Pipes::compute_Pantheon_LogLike;

      const str path_to_file = GAMBIT_DIR "/CosmoBit/data/SNeIa/Pantheon/lcparam_full_long.dat";
      const str path_to_covmat = GAMBIT_DIR "/CosmoBit/data/SNeIa/Pantheon/sys_full_long.dat";
      static ASCIItableReader data(path_to_file);
      static ASCIItableReader covmat_data(path_to_covmat);
      static int nobs = covmat_data[0][0];
      static bool read_data = false;

      static gsl_matrix *invcov = gsl_matrix_alloc(nobs, nobs);
      gsl_vector *residuals = gsl_vector_alloc(nobs);
      gsl_vector *product = gsl_vector_alloc(nobs);

      int ie,je;
      double dl,chi2;
      double M = *Param["M"]; 

      if(read_data == false)
      {
        read_data = true;
        logger() << "Pantheon data read from file '"<<path_to_file<<"'." << EOM;
        std::vector<std::string> colnames = initVector<std::string>("zcmb", "zhel", "dz", "mb", "dmb", "x1", "dx1", "color", "dcolor", "3rdvar", "d3rdvar", "cov_m_s", "cov_m_c", "cov_s_c", "set", "ra", "dec");
        data.setcolnames(colnames);

        int covmat_line = 1; // 0-th entry stores dimension of covmat, start with 1

        for(ie=0;ie<nobs;ie++)  // fill covmat and invert after that
        {
          for(je=0;je<nobs;je++)
          {
            if (ie == je)   // add statistical error of mb to diagonal elements
            {
              gsl_matrix_set(invcov, ie, je, covmat_data[0][covmat_line]+data["dmb"][ie]*data["dmb"][ie]);
              covmat_line+=1;
            }
            else
            {
              gsl_matrix_set(invcov, ie, je, covmat_data[0][covmat_line]);
              covmat_line+=1;
            }
          }
        }

        gsl_linalg_cholesky_decomp(invcov);
        gsl_linalg_cholesky_invert(invcov); // now the inverse of the covariance matrix is stored in invcov
        // delete read covmat ASCIItableReader covmat_data here?
      }

      for(ie=0;ie<nobs;ie++)
      {
        dl = BEreq::class_get_Dl(byVal(data["zcmb"][ie]));
        gsl_vector_set(residuals, ie, data["mb"][ie] - (5 * log10(dl) + 25 + M));
      }

      gsl_blas_dgemv(CblasNoTrans, 1.0 , invcov, residuals, 0.,  product); // solves product = 1 x invcov residuals + 0 x product
      gsl_blas_ddot(residuals, product, &chi2); // computes chi2 = residuals^T product = residuals^T invcov residuals

      result = -0.5* chi2;
      // do not free invcov since it is needed for all points
      gsl_vector_free(residuals);
      gsl_vector_free(product);
    }

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

    void compute_BAO_LogLike(double &result)
    {
      using namespace Pipes::compute_BAO_LogLike;

      int type;
      double da,dr,dv,rs,Hz,z,theo; // cosmo distances and sound horizon at drag,Hubble, theoretical prediction
      double chi2 =0;

      const str filename = runOptions->getValue<std::string>("DataFile");
      const str path_to_file = GAMBIT_DIR "/CosmoBit/data/BAO/" + runOptions->getValue<std::string>("DataFile");
      static ASCIItableReader data(path_to_file);
      static bool read_data = false;
      static int nrow;

      if(read_data == false)
      {
        logger() << "BAO data read from file '"<<filename<<"'." << EOM;
        std::vector<std::string> colnames = initVector<std::string>("z", "mean","sigma","type");
        data.setcolnames(colnames);
        nrow=data.getnrow();

        read_data = true;
      }

      for(int ie=0;ie<nrow;ie++)
      {
        z = data["z"][ie];
        type = data["type"][ie];
        da = BEreq::class_get_Da(byVal(z)); // angular diam. dist as function of redshift
        Hz = BEreq::class_get_Hz(byVal(z)); // Hubble parameter as function of redshift (in 1/Mpc)
        dr = z/Hz;
        dv = pow(da*da*(1.+z)*(1.+z)*dr, 1./3.);
        rs = BEreq::class_get_rs();

        switch(type)
        {
          case 3: theo = dv /rs; break;
          case 4: theo = dv; break;
          case 5: theo = da / rs; break;
          case 6: theo = 1./Hz/rs; break;
          case 7: theo = rs/dv; break;
          default:
            std::ostringstream err;
            err << "Type " << type << " in "<< ie+1 <<". data point in BAO data file '"<<filename <<"' not recognised.";
            CosmoBit_error().raise(LOCAL_INFO, err.str());
        }
        //pybind11::print("         (MOntePython internal) compute BAO for z = ", z, ": ",theo," ",data["mean"][ie]," ",data["sigma"][ie]);
        chi2 += pow((theo - data["mean"][ie]) / data["sigma"][ie],2);
      }
      result = -0.5*chi2;
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

      sigma8 = BEreq::class_get_sigma8(0.);
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
    void init_Classy_cosmo_container(CosmoBit::Classy_cosmo_container& ccc)
    {
      using namespace Pipes::init_Classy_cosmo_container;

      // create instance of classy class Class()
      BEreq::classy_create_class_instance(ccc.cosmo);

      pybind11::dict cosmo_input_dict = *Dep::set_classy_parameters;
      ccc.set_input_dict(cosmo_input_dict);

      BEreq::classy_compute(ccc);

      std::cout << "(CosmoBit): initialised cosmo object & computed vals"<< std::endl;
    }
    
    /// Initialises the container within CosmoBit from classy, but designed specifically
    /// to be used when MontePython is in use. This will ensure additional outputs are
    /// computed by classy CLASS to be passed to MontePython.
    void init_Classy_cosmo_container_with_MPLike(CosmoBit::Classy_cosmo_container& ccc)
    {
      using namespace Pipes::init_Classy_cosmo_container_with_MPLike;

      // create instance of classy class Class()
      BEreq::classy_create_class_instance(ccc.cosmo);

      // get extra cosmo_arguments from MP (gives a dictionary with output values that need
      // to be set for the class run)
      static pybind11::dict MP_cosmo_arguments = *Dep::cosmo_args_from_MPLike;
      pybind11::print("  ----|> Extra cosmo_arguments needed from MP Likelihoods: ",MP_cosmo_arguments);

      pybind11::dict cosmo_input_dict = *Dep::set_classy_parameters;

      // add the arguments from Mp_cosmo_arguments which are not yet in cosmo_input_dict to it
      // also takes care of merging the "output" key values
      merge_pybind_dicts(cosmo_input_dict,MP_cosmo_arguments);

      pybind11::print("  ----|> Combined dict", cosmo_input_dict);

      // set the input dict of Classy cosmo container
      ccc.set_input_dict(cosmo_input_dict);

      // actual CLASS run through classy
      BEreq::classy_compute(ccc);

      std::cout << "(CosmoBit): initialised cosmo object & computed vals"<< std::endl;
    }

    /// Set the LCDM parameters in classy. Looks at the parameters used in a run,
    /// and passes them to classy in the form of a Python dictionary.
    void set_classy_parameters_LCDM(pybind11::dict& result)
    {
      using namespace Pipes::set_classy_parameters_LCDM;
      using namespace pybind11::literals;


      map_str_dbl NuMasses_SM = *Dep::NuMasses_SM;
      int N_ncdm = (int)NuMasses_SM["N_ncdm"];

      // Number of non-cold DM species
      if (N_ncdm > 0.)
      {
        result["N_ncdm"] = N_ncdm;   

        std::vector<double> m_ncdm = m_ncdm_classInput(NuMasses_SM);
        std::vector<double> T_ncdm(N_ncdm,*Dep::T_ncdm);

        std::ostringstream ss1, ss2;
        std::string separator;
        for (auto x : m_ncdm) {
          ss1 << separator << x;
          separator = ",";
        }
        separator = "";
        for (auto x : T_ncdm) {
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
      result["N_ur"] = *Dep::class_Nur; // Number of ultra relativistic species
      result["output"] = runOptions->getValueOrDef<str>("", "output");
      // JR we should not set this by default when using the python wrapper: if only bg needed 
      // classy will complain that it did not read the input parameter
      // if we use Planck through MP this should be taken care of by itself anyway
      //result["l_max_scalars"] = runOptions->getValueOrDef<int>(2508., "l_max"); 
      result["lensing"] = runOptions->getValueOrDef<str>("no","lensing");
      result["T_cmb"] = *Dep::T_cmb;
      result["omega_b"] = *Param["omega_b"];
      result["omega_cdm"] = *Param["omega_cdm"];
      result["H0"] = *Param["H0"];
      result["ln10^{10}A_s"] = *Param["ln10A_s"];
      result["n_s"] = *Param["n_s"];
      result["tau_reio"] = *Param["tau_reio"];
      std::vector<double> Helium_abundance = *Dep::Helium_abundance;
      result["YHe"] = Helium_abundance.at(0); // .at(0): mean, .at(1): uncertainty


      
      // Other Class input direct from the YAML file
      YAML::Node classy_dict;
      if (runOptions->hasKey("classy_dict"))
      {
        classy_dict = runOptions->getValue<YAML::Node>("classy_dict");
        for (auto it=classy_dict.begin(); it != classy_dict.end(); it++)
        {
          std::string name = it->first.as<std::string>();
          std::string value = it->second.as<std::string>();
          // Check if the key exists in the dictionary
          if (!result.contains(name.c_str()))
          {
            result[name.c_str()] = value;
          }
          // If it does, throw an error, there's some YAML conflict going on.
          else
          {
            CosmoBit_error().raise(LOCAL_INFO, "The key " + name + " already "
              "exists in the CLASSY dictionary. Please check your YAML file.");
          }

        }
      }

      std::cout << "(CosmoBit): set ccc.cosmo_input_dict values to "<< std::endl;
      pybind11::print("      ",result);

    }


    /* Classy getter functions */

    /// Energy densities *today* (Omega0)

    /// Matter
    void get_Omega0_m_classy(double& result)
    {
      using namespace Pipes::get_Omega0_m_classy;

      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      result = ccc.cosmo.attr("Omega0_m")().cast<double>();
    }

    /// Baryons
    void get_Omega0_b_classy(double& result)
    {
      using namespace Pipes::get_Omega0_b_classy;

      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      result = ccc.cosmo.attr("Omega_b")().cast<double>();
    }

    /// Cold Dark Matter
    void get_Omega0_cdm_classy(double& result)
    {
      using namespace Pipes::get_Omega0_cdm_classy;

      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      // TODO -- does not currently work
      //result = ccc.cosmo.attr("Omega0_cdm")().cast<double>();
      result = 0.;    }

    /// Radiation
    void get_Omega0_r_classy(double& result)
    {
      using namespace Pipes::get_Omega0_r_classy;

      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      // TODO -- does not currently work
      //result = ccc.cosmo.attr("Omega0_r")().cast<double>();
      result = 0.;
    }

    /// Photons
    void get_Omega0_g_classy(double& result)
    {
      using namespace Pipes::get_Omega0_g_classy;

      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      result = ccc.cosmo.attr("Omega_g")().cast<double>();
    }

    /// Ultra-relativistic
    void get_Omega0_ur_classy(double& result)
    {
      using namespace Pipes::get_Omega0_ur_classy;

      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      // TODO -- does not currently work
      //result = ccc.cosmo.attr("Omega0_ur")().cast<double>();
      result = 0.;
    }

    /// Non-cold Dark Matter
    void get_Omega0_ncdm_classy(double& result)
    {
      using namespace Pipes::get_Omega0_ncdm_classy;

      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      // TODO -- does not currently work
      // result = ccc.cosmo.attr("Omega0_ncdm")().cast<double>();
      result = 0.;
    }

    /// Other observables.

    /// Sigma8
    void get_Sigma8_classy(double& result)
    {
      using namespace Pipes::get_Sigma8_classy;

      // We need to check if the matter power spectrum is
      // being computed or not. If not, throw an error.
      pybind11::dict classy_params = *Dep::set_classy_parameters;
      pybind11::str output = classy_params["output"];
      if (output.attr("find")("mPk").cast<int>() == -1) {
        CosmoBit_error().raise(LOCAL_INFO, "CLASSY can not compute "
          "σ8 unless you ask for the matter power spectrum.\n"
          "Try adding the following entry to your YAML file:\n"
          "  - capability: set_classy_parameters\n"
          "    function: set_classy_parameters_LCDM\n"
          "    options:\n"
          "        output: \"lCl pCl tCl mPk\"");
      }

      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      result = ccc.cosmo.attr("sigma8")().cast<double>();
    }

    /// Effective number of neutrino species
    // (mostly for cross-checking!)
    void get_Neff_classy(double& result)
    {
      using namespace Pipes::get_Neff_classy;

      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      result = ccc.cosmo.attr("Neff")().cast<double>();
    }

    /// Comoving sound horizon at Baryon drag epoch
    void get_rs_drag_classy(double& result)
    {
      using namespace Pipes::get_rs_drag_classy;

      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      result = ccc.cosmo.attr("rs_drag")().cast<double>();
    }

    /// for the future, maybe -- or to live in CB_utils/frontend...
    /// Angular distance
    /*
    void get_angular_distance_classy(double& result)
    {
      using namespace Pipes::get_angular_distance_classy;
      double z = ...?
      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      result = ccc.cosmo.attr("angular_distance")(z).cast<double>();
    }    
    void get_luminosity_distance_classy(double& result)
    {
      using namespace Pipes::get_luminosity_distance_classy;
      double z = ...?
      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      result = ccc.cosmo.attr("luminosity_distance")(z).cast<double>();
    }
    void get_scale_invariant_growth_factor(double& result)
    {
      using namespace Pipes::get_scale_invariant_growth_factor;
      double z = ...?
      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;
      result = ccc.cosmo.attr("scale_invariant_growth_factor")(z).cast<double>();
    }
    */

/***************/
/* MontePython */
/***************/

    /// function to fill the mcmc_parameters dictionary of MontePython's Data object with current 
    /// values of nuisance parameters
    void set_parameter_dict_for_MPLike(pybind11::dict & result)
    {
      using namespace Pipes::set_parameter_dict_for_MPLike;
      using namespace pybind11::literals;

      for (auto it=Param.begin(); it != Param.end(); it++)
      {
        std::string name = it->first;
        double value = *Param[name];
        result[name.c_str()] = pybind11::dict("current"_a=*Param[name],"scale"_a=1.); // scale always 1 in GAMBIT
        std::cout<<"    (CosmoBit) Nuisance parameters are "<< name << " with value " << byVal(value)<<std::endl;
      }
      std::cout<<"    (CosmoBit) Using Likelihood with cosmological nuisance parameters.. setting values for MPLike data.mcmc_parameters dict"<<std::endl;      
    }

    /// function to fill the mcmc_parameters dictionary of MontePython's Data object with current 
    /// This version of the capability 'parameter_dict_for_MPLike' is used when no Likelihood with nuisance parameters is in use
    /// just passes an empty py dictionary
    void pass_empty_parameter_dict_for_MPLike(pybind11::dict & result)
    {
      using namespace Pipes::pass_empty_parameter_dict_for_MPLike;
      using namespace pybind11::literals;

      std::cout<<"    (CosmoBit) No Likelihood with cosmological nuisance parameters in use.. doing nothing in parameter_dict_for_MPLike "<<std::endl;      
    }

    /// Set the names of all experiments in use for scan 
    /// basically reads in the runOption 'Likelihoods' of capability MP_experiment_names
    void set_MP_experiment_names(std::vector<std::string> & result)
    {
      using namespace Pipes::set_MP_experiment_names;

      static bool first = true;
      if (first)
      {

        // Read analysis names from the yaml file
        std::vector<std::string> default_analyses;  // The default is empty lists of analyses TODO : deal with what happens if no Like is specified, throw error
        result = runOptions->getValueOrDef<std::vector<std::string> >(default_analyses, "Likelihoods");

        // Check that the analysis names listed in the yaml file all correspond to actual ColliderBit analyses
        for (std::string & name : result)
        {
          std::cout << "(CosmoBit) Read MontePythonLike option  "<< name <<std::endl;
            
          /*if (!checkAnalysis(analysis_name)) // TODO add check if likelihood implemented (similar to what we did in MPLike patch file)
          {
            str errmsg = "The analysis " + analysis_name + " is not a known ColliderBit analysis.";
            ColliderBit_error().raise(LOCAL_INFO, errmsg);
          }*/
          
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

      std::cout << "(CosmoBit) entered init_cosmo_args_from_MPLike"<< std::endl;

      // Need to know where CLASSY lives
      std::string classyDir = BEreq::path_to_classy();

      // CosmoBit::MPLike_data_container should only be created once when calculating the first point.
      // After that is has to be kept alive since it contains a vector with the initialised MPLike 
      // Likelihood objects.
      clock_t tStart = clock();
      static bool first_run = true;
      static pybind11::object data;
      if(first_run)
      {
        std::vector<std::string> experiments = *Dep::MP_experiment_names;
        data = BEreq::create_data_object(experiments,classyDir);
        map_str_pyobj likelihoods = BEreq::create_likelihood_objects(data, experiments);
        first_run = false;
      }
      
      result = data.attr("cosmo_arguments");

      std::cout << "(CosmoBit) time took to create data & likelihood objects "<< (double)(clock() - tStart)/CLOCKS_PER_SEC<<std::endl;
      pybind11::print("     cosmo_arguments are ", result);

    }

    /// Computes lnL for each experiment initialised in MontePython
    void calc_MP_LogLikes(map_str_dbl & result)
    {
      using namespace Pipes::calc_MP_LogLikes;
      using namespace pybind11::literals;

      // A list of the experiments initialised in the YAML file
      std::vector<std::string> experiments = *Dep::MP_experiment_names;

      std::cout << "(CosmoBit): init_MontePythonLike start"<< std::endl;

      // Need to know where CLASSY lives
      std::string classyDir = BEreq::path_to_classy();
      
      // CosmoBit::MPLike_data_container should only be created once when calculating the first point.
      // After that is has to be kept alive since it contains a vector with the initialised MPLike 
      // Likelihood objects.
      pybind11::object data;
      map_str_pyobj likelihoods;
      static bool first_run = true;
      if(first_run)
      {
        data = BEreq::create_data_object(experiments, classyDir);
        likelihoods = BEreq::create_likelihood_objects(data, experiments);
        first_run = false;
      }
      
      static const CosmoBit::MPLike_data_container mplike_cont(data, likelihoods);
      
      // pass current values of nuisance parameters to data.mcmc_parameters dictionary for likelihood computation in MP
      mplike_cont.data.attr("mcmc_parameters") = *Dep::parameter_dict_for_MPLike;

      // Test if cosmo_arguments & mcmc_parameters dictionaries contain (nuisance) parameters that need to go into ccc.cosmo_input_dict 
      // and likelihood calculation (values depend on used likelihoods)  --- Works! 
      pybind11::print("(CosmoBit) data.cosmo_arguments contains ", mplike_cont.data.attr("cosmo_arguments"));
      pybind11::print("(CosmoBit) data.mcmc_parameters contains ", mplike_cont.data.attr("mcmc_parameters"));

      // Create instance of classy class Class
      CosmoBit::Classy_cosmo_container ccc = *Dep::get_Classy_cosmo_container;

      // Loop through the list of experiments, and query the lnL from the
      // MontePython backend
      for (std::vector<std::string>::const_iterator it = experiments.begin();
        it != experiments.end(); ++it)
      {
        std::string like_name = *it;
        result[like_name] = BEreq::get_MP_loglike(mplike_cont, ccc.cosmo, like_name);
      }

      std::cout << "(CosmoBit): get_MP_loglike end with result \n"<< result << std::endl;
    }

    /// Computes the combined lnL from the set of experiments 
    /// given to MontePython.
    void calc_MP_combined_LogLike(double& result)
    {
      using namespace Pipes::calc_MP_combined_LogLike;

      map_str_dbl MP_lnLs = *Dep::MP_LogLikes;

      // Iterate through map of doubles and return one big fat double.
      double lnL = 0.;
      for (const auto &p : MP_lnLs)
      {
        std::cout << "Likelihood name: "  << p.first << std::endl;
        std::cout << "Likelihood value: " << p.second<< std::endl;
        lnL += p.second;
      }
      result = lnL;
    }
  }
}

