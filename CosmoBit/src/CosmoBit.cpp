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
///
///  *********************************************
#include <string>
#include <iostream>
#include <cmath>
#include <functional>

#include <omp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>


#include <stdlib.h>     /* malloc, free, rand */
#include "gambit/Utils/yaml_options.hpp"
#include "gambit/Utils/statistics.hpp"
#include "gambit/Utils/ascii_table_reader.hpp"
#include "gambit/Utils/ascii_dict_reader.hpp"
#include "gambit/Utils/numerical_constants.hpp"
#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/CosmoBit/CosmoBit_rollcall.hpp"
#include "gambit/CosmoBit/CosmoBit_types.hpp"


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

    void DM_mass_ToyModel(double& result)
    {
      using namespace Pipes::DM_mass_ToyModel;
      result = *Param["mass"];
    }

    void lifetime_ToyModel(double& result)
    {
      using namespace Pipes::lifetime_ToyModel;
      result = *Param["lifetime"];
    }

    void DM_fraction_ToyModel(double& result)
    {
      using namespace Pipes::DM_fraction_ToyModel;
      result = *Param["fraction"];
    }

    void f_effective_func(double& result)
    {
      using namespace Pipes::f_effective_func;

      bool silent = runOptions->getValueOrDef<bool>(false,"silent_mode");
      int last_steps = runOptions->getValueOrDef<int>(4,"show_last_steps");
      double z_eff = runOptions->getValueOrDef<double>(600.,"z_eff");

      DarkAges::fz_table fzt = BEreq::DA_efficiency_function();
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

    void PyTest_func_1(double& result)
    {
      using namespace Pipes::PyTest_func_1;

      std::vector<double> vec = BEreq::PyArrayTest_Py_to_cpp();
      double sum = 0.0;
      for (std::vector<double>::iterator it = vec.begin(); it != vec.end(); ++it)
      {
        sum += *it;
      }
      result = sum;
    }

    void PyTest_func_2(double& result)
    {
      using namespace Pipes::PyTest_func_2;

      int len = runOptions->getValueOrDef<int>(10,"array_len");
      result = BEreq::PyArrayTest_cpp_to_Py(byVal(len));
    }

    void lnL_A_planck_gaussian(double& result)
    {
      using namespace Pipes::lnL_A_planck_gaussian;

      double obs_mean = 1.0;
      double obs_err = 0.0025;
      double pred_mean = *Param["A_planck"];
      double pred_err = 0.0;

      result = Stats::gaussian_loglikelihood(pred_mean, obs_mean, pred_err, obs_err, false);
    }

    void class_set_parameter_LCDM(Class_container& cosmo)
    {
      using namespace Pipes::class_set_parameter_LCDM;

      int l_max=cosmo.lmax;

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

      std::vector<double> Helium_abund = *Dep::Helium_abundance; // .at(0): mean, .at(1): uncertainty
      cosmo.input.addEntry("YHe",Helium_abund.at(0)); 

      // JR: heads-up! Neutrinos + dNeff also need to be set. accordingly for base model we want to scan
      // Planck baseline should be 
      cosmo.input.addEntry("N_ur",2.0328);  //1 massive neutrinos 
      cosmo.input.addEntry("N_ncdm",1);
      cosmo.input.addEntry("m_ncdm","0.06");


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

    void class_set_parameter_LCDM_dNeff_Smu(Class_container& cosmo)
    {
      using namespace Pipes::class_set_parameter_LCDM_dNeff_Smu;
      // on the level of class the model LCDM_dNeff_Smu and LCDM_dNeff_Smu_etaBBN are identical

      int l_max=cosmo.lmax;

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

      cosmo.input.addEntry("N_ur",*Param["dNeff"]+0.00641);  // dNeff= 0.00641 for 3 massive neutrinos at CMB release
      cosmo.input.addEntry("N_ncdm",3);
      std::stringstream sstream;
      double numass = *Param["Smu"]/3.;
      sstream << numass << ", " << numass << ", " << numass;
      cosmo.input.addEntry("m_ncdm",sstream.str());

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

      // Maximal value of l.
      int l_max = cosmo.lmax;

      // Number of Cl-spectra (columns of the Cl-table).
      // The order of the spectra is [TT, EE, TE, BB, PhiPhi, TPhi, EPhi]
      int num_ct_max=7;

      // Define an array which takes the values of Cl of the different
      // spectra at a given value of l
      double* cl = new double[num_ct_max];

      // Loop through all l from 0 to l_max (including) and ask for the cl-spectra.
      for (int l=0; l <= l_max; l++)
      {
        if (l < 2)
        {
          // The entries for l=0 and l=1 are zero per defintion
          cosmo.Cl_TT.at(l) = 0.;
          cosmo.Cl_TE.at(l) = 0.;
          cosmo.Cl_EE.at(l) = 0.;
          cosmo.Cl_BB.at(l) = 0.;
          cosmo.Cl_PhiPhi.at(l) = 0.;
        }
        else
        {
          // For l >= 2 ask for the cl-spectra.
          if (BEreq::class_output_total_cl_at_l(&cosmo.sp,&cosmo.le,&cosmo.op,byVal(l),byVal(cl)) == _SUCCESS_)
          {
            cosmo.Cl_TT.at(l) = cl[cosmo.sp.index_ct_tt]*pow(cosmo.ba.T_cmb*1.e6,2);
            cosmo.Cl_TE.at(l) = cl[cosmo.sp.index_ct_te]*pow(cosmo.ba.T_cmb*1.e6,2);
            cosmo.Cl_EE.at(l) = cl[cosmo.sp.index_ct_ee]*pow(cosmo.ba.T_cmb*1.e6,2);
            cosmo.Cl_BB.at(l) = cl[cosmo.sp.index_ct_bb]*pow(cosmo.ba.T_cmb*1.e6,2);
            cosmo.Cl_PhiPhi.at(l) = cl[cosmo.sp.index_ct_pp];
          }
          else
          {
            // Failsafe for unexpected behaviour of "class_outpout_at_cl"
            cosmo.Cl_TT.at(l) = 0.;
            cosmo.Cl_TE.at(l) = 0.;
            cosmo.Cl_EE.at(l) = 0.;
            cosmo.Cl_BB.at(l) = 0.;
            cosmo.Cl_PhiPhi.at(l) = 0.;
          }
        }
      }

      // We do not need "cl" anymore
      delete cl;
    }

/* Begin of outdated functions (might be removed) */
/* end of outdated functions */

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

/// -----------
/// BBN related functions
/// -----------

    void set_T_cmb(double &result)
    {
      using namespace Pipes::set_T_cmb;

      // use COBE/FIRAS (0911.1955) mes. if not specified otherwise
      result = runOptions->getValueOrDef<double>(2.7255,"T_cmb"); 

    }

    void calculate_eta(double &result)
    {
      using namespace Pipes::calculate_eta;

      double ngamma, nb;
      ngamma = 16*pi*zeta3*pow(*Dep::T_cmb*kb/hc,3); // photon number density today
      nb = *Param["omega_b"]*3*100*1e3*100*1e3/Mpc/Mpc/(8*pi*Gn*m_proton_g); // baryon number density today

      result =  nb/ngamma;
      logger() << "Baryon to photon ratio (eta) computed to be " << result << EOM;
    }

    void AlterBBN_fill(relicparam &result) 
    {
      // fill AlterBBN structure for LCDM
      using namespace Pipes::AlterBBN_fill;

      BEreq::Init_cosmomodel(&result);

      result.eta0=*Dep::eta; // eta0 = eta_BBN = eta_CMB
      result.Nnu=3.046;     // 3 massive neutrinos
      result.dNnu=0;        // no extra rel. d.o.f. in  base LCDM
      result.failsafe = 3;  // set precision parameters for AlterBBN
      result.err = 3;
    }

    void AlterBBN_fill_dNeff(relicparam &result)
    {
      // fill AlterBBN structure for LCDM_dNeff_Smu
      using namespace Pipes::AlterBBN_fill_dNeff;

      BEreq::Init_cosmomodel(&result);

      result.eta0=*Dep::eta; // eta0 = eta_BBN = eta_CMB
      result.Nnu=3.046;            // 3 massive neutrinos
      result.dNnu=*Param["dNeff"]; // no extra rel. d.o.f. in  base LCDM
      result.failsafe = 3;         // set precision parameters for AlterBBN
      result.err = 3;
    }


    void AlterBBN_fill_etaBBN(relicparam &result)
    {
      // fill AlterBBN structure for LCDM__dNeff_Smu_etaBBN
      using namespace Pipes::AlterBBN_fill_etaBBN;

      BEreq::Init_cosmomodel(&result);
      
      result.eta0 = *Param["eta_BBN"];  // eta AFTER BBN (variable during)
      result.Nnu=3.046;                 // 3 massive neutrinos
      result.dNnu=*Param["dNeff_BBN"];
      result.failsafe = 3;              // set precision parameters for AlterBBN
      result.err = 3; 
    }


    void compute_BBN_abundances(CosmoBit::BBN_container &result)
    {
      using namespace Pipes::compute_BBN_abundances;

      int NNUC = BEreq::get_NNUC();  // global variable of AlterBBN (# computed element abundances)  
      result.init_arr(NNUC);         // init arrays in BBN_container with right length
      double ratioH [NNUC+1], cov_ratioH [NNUC+1][NNUC+1];
      
      relicparam const& paramrelic = *Dep::AlterBBN_modelinfo;
      int nucl_err = BEreq::nucl_err(&paramrelic, byVal(ratioH), byVal(cov_ratioH[0]));
      
      // fill abundances and covariance matrix of BBN_container with results from AlterBBN
      if (nucl_err)
      {
        for(int ie=1;ie<=NNUC;ie++)
        {
          result.BBN_abund.at(ie) = ratioH[ie];
          for(int je=1;je<=NNUC;je++) result.BBN_covmat.at(ie).at(je) =cov_ratioH[ie][je];
        }
      }
      else
      {
        std::ostringstream err;
        err << "AlterBBN calculation for primordial element abundances failed.";
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

      // compute chi2
      for(ie=0;ie<nobs;ie++) for(je=0;je<nobs;je++) chi2+=(prediction[ie]-observed[ie])*gsl_matrix_get(invcov,ie,je)*(prediction[je]-observed[je]);
      result = -0.5*chi2;  
      
      gsl_matrix_free(cov);
      gsl_matrix_free(invcov);
      gsl_permutation_free(p);
    }

/// -----------
/// Background Likelihoods (H0 + BAO)
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
      
      Class_container cosmo = *Dep::class_get_spectra;
    
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
        rs = cosmo.th.rs_d;

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
        chi2 += pow((theo - data["mean"][ie]) / data["sigma"][ie],2);
      }
      result = -0.5*chi2;
    }

    void compute_sigma8_LogLike(double &result)
    {
      using namespace Pipes::compute_sigma8_LogLike;

      double sigma8,Omega_m; // cosmo distances and sound horizon at drag,Hubble, theoretical prediction
      double chi2 =0;
      
      const str filename = runOptions->getValue<std::string>("DataFile");
      const str path_to_file = GAMBIT_DIR "/CosmoBit/data/sigma8/" + runOptions->getValue<std::string>("DataFile");
      static ASCIItableReader data(path_to_file);
      static bool read_data = false;
      static int nrow;
      
      Class_container cosmo = *Dep::class_get_spectra;
      sigma8 = BEreq::class_get_sigma8(0.);
      Omega_m = cosmo.ba.Omega0_cdm+cosmo.ba.Omega0_dcdmdr+cosmo.ba.Omega0_b+cosmo.ba.Omega0_ncdm_tot;
    
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
        chi2 += pow(( sigma8*pow(Omega_m/data["Omega_m_ref"][ie], data["Omega_m_index"][ie]) - data["bestfit"][ie])/ data["sigma"][ie],2);
      }
      result = -0.5*chi2;
    }

  }

}
