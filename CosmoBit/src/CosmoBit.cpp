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
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2017 Nov
///  \date 2018 Jan,Feb, Mar
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

#include <stdlib.h>     /* malloc, free, rand */
#include "gambit/Utils/yaml_options.hpp"
#include "gambit/Utils/statistics.hpp"
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

      std::cout << "pot = " << pot << std::endl;

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
      //std::cout << "Last seen alive in: class_set_parameter_LCDM" << std::endl;
      using namespace Pipes::class_set_parameter_LCDM;

      int l_max=cosmo.lmax;

      // What follows is a loop to read all generic additional parameters for CLASS
      // needed to (hopefully) reproduce the LCDM values of Planck, which we pass
      // as options a YAML-node named "class_dict" within the 'Rules' section of the input,
      // and which will be added if there are present.
      int len_of_input = 9;
      Options class_dict;
      std::vector<str> names;
      if (runOptions->hasKey("class_dict"))
      {
        class_dict = Options(runOptions->getValue<YAML::Node>("class_dict"));
        names = class_dict.getNames();
        len_of_input += names.size();
      }

      BEreq::class_parser_initialize(&cosmo.fc,byVal(len_of_input),"",cosmo.class_errmsg);

      strcpy(cosmo.fc.name[0],"output");
      strcpy(cosmo.fc.value[0],"tCl pCl lCl");
      strcpy(cosmo.fc.name[7],"l_max_scalars");
      sprintf(cosmo.fc.value[7],"%d",l_max);
      strcpy(cosmo.fc.name[8],"lensing");
      strcpy(cosmo.fc.value[8],"yes");

      int i = 9;

      if (runOptions->hasKey("class_dict"))
      {
        for (std::vector<str>::iterator name_it = names.begin(); name_it != names.end(); name_it++)
        {
          //std::cout << "Key = " << *name_it << " and Value = " << class_dict.getValue<str>(*name_it) << std::endl;
          str key, value;
          key = *name_it;
          value = class_dict.getValue<str>(*name_it);
          sprintf(cosmo.fc.name[i],"%s",key.c_str());
          sprintf(cosmo.fc.value[i],"%s",value.c_str());
          i++;
        }
      }

      strcpy(cosmo.fc.name[1],"omega_b");
      strcpy(cosmo.fc.name[2],"omega_cdm");
      strcpy(cosmo.fc.name[3],"H0");
      strcpy(cosmo.fc.name[4],"ln10^{10}A_s");
      strcpy(cosmo.fc.name[5],"n_s");
      strcpy(cosmo.fc.name[6],"tau_reio");

      sprintf(cosmo.fc.value[1],"%e",*Param["omega_b"]);
      sprintf(cosmo.fc.value[2],"%e",*Param["omega_cdm"]);
      sprintf(cosmo.fc.value[3],"%e",*Param["H0"]);
      sprintf(cosmo.fc.value[4],"%e",*Param["ln10A_s"]);
      sprintf(cosmo.fc.value[5],"%e",*Param["n_s"]);
      sprintf(cosmo.fc.value[6],"%e",*Param["tau_reio"]);
      /*
      std::cout << "omega_b = " << *Param["omega_b"] << std::endl;
      std::cout << "omega_cdm = " << *Param["omega_cdm"] << std::endl;
      std::cout << "H0 = " << *Param["H0"] << std::endl;
      std::cout << "ln10A_s = " << *Param["ln10A_s"] << std::endl;
      std::cout << "n_s = " << *Param["n_s"] << std::endl;
      std::cout << "tau_reio = " << *Param["tau_reio"] << std::endl;
      */
    }

    void class_set_parameter_LCDM_SingletDM(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_LCDM_SingletDM" << std::endl;
      using namespace Pipes::class_set_parameter_LCDM_SingletDM;

      int l_max=cosmo.lmax;
      double sigmav = *Dep::sigmav; // in cm^3 s^-1
      double mass = *Dep::mwimp; // in GeV
      double feff = runOptions->getValueOrDef<double>(1.,"f_eff");
      double annihilation = (1.0/1.78e-21)*(sigmav/mass)*feff; // in m^3 s^-1 kg^-1

      // What follows is a loop to read all generic additional parameters for CLASS
      // needed to (hopefully) reproduce the LCDM values of Planck, which we pass
      // as options a YAML-node named "class_dict" within the 'Rules' section of the input,
      // and which will be added if there are present.

      int len_of_input = 10;
      Options class_dict;
      std::vector<str> names;

      if (runOptions->hasKey("class_dict"))
      {
        class_dict = Options(runOptions->getValue<YAML::Node>("class_dict"));
        names = class_dict.getNames();
        len_of_input += names.size();
      }

      BEreq::class_parser_initialize(&cosmo.fc,byVal(len_of_input),"",cosmo.class_errmsg);

      strcpy(cosmo.fc.name[0],"output");
      strcpy(cosmo.fc.value[0],"tCl pCl lCl");
      strcpy(cosmo.fc.name[8],"l_max_scalars");
      sprintf(cosmo.fc.value[8],"%d",l_max);
      strcpy(cosmo.fc.name[9],"lensing");
      strcpy(cosmo.fc.value[9],"yes");

      int i = 10;

      if (runOptions->hasKey("class_dict"))
      {
        for (std::vector<str>::iterator name_it = names.begin(); name_it != names.end(); name_it++)
        {
          //std::cout << "Key = " << *name_it << " and Value = " << class_dict.getValue<str>(*name_it) << std::endl;
          str key, value;
          key = *name_it;
          value = class_dict.getValue<str>(*name_it);
          sprintf(cosmo.fc.name[i],"%s",key.c_str());
          sprintf(cosmo.fc.value[i],"%s",value.c_str());
          i++;
        }
      }

      strcpy(cosmo.fc.name[1],"omega_b");
      strcpy(cosmo.fc.name[2],"omega_cdm");
      strcpy(cosmo.fc.name[3],"H0");
      strcpy(cosmo.fc.name[4],"ln10^{10}A_s");
      strcpy(cosmo.fc.name[5],"n_s");
      strcpy(cosmo.fc.name[6],"tau_reio");
      strcpy(cosmo.fc.name[7],"annihilation");

      sprintf(cosmo.fc.value[1],"%e",*Param["omega_b"]);
      sprintf(cosmo.fc.value[2],"%e",*Param["omega_cdm"]);
      sprintf(cosmo.fc.value[3],"%e",*Param["H0"]);
      sprintf(cosmo.fc.value[4],"%e",*Param["ln10A_s"]);
      sprintf(cosmo.fc.value[5],"%e",*Param["n_s"]);
      sprintf(cosmo.fc.value[6],"%e",*Param["tau_reio"]);
      sprintf(cosmo.fc.value[7],"%e",annihilation);
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

      int l_max=cosmo.lmax;;

      // What follows is a loop to read all generic additional parameters for CLASS
      // needed to (hopefully) reproduce the LCDM values of Planck, which we pass
      // as options a YAML-node named "class_dict" within the 'Rules' section of the input,
      // and which will be added if there are present.

      int len_of_input = 11;
      Options class_dict;
      std::vector<str> names;

      if (runOptions->hasKey("class_dict"))
      {
        class_dict = Options(runOptions->getValue<YAML::Node>("class_dict"));
        names = class_dict.getNames();
        len_of_input += names.size();
      }

      BEreq::class_parser_initialize(&cosmo.fc,byVal(len_of_input),"",cosmo.class_errmsg);

      strcpy(cosmo.fc.name[0],"output");
      strcpy(cosmo.fc.value[0],"tCl pCl lCl");
      strcpy(cosmo.fc.name[8],"l_max_scalars");
      sprintf(cosmo.fc.value[8],"%d",l_max);
      strcpy(cosmo.fc.name[9],"modes");
      strcpy(cosmo.fc.value[9],"s,t");
      strcpy(cosmo.fc.name[10],"lensing");
      strcpy(cosmo.fc.value[10],"yes");

      int i = 11;

      if (runOptions->hasKey("class_dict"))
      {
        for (std::vector<str>::iterator name_it = names.begin(); name_it != names.end(); name_it++)
        {
          //std::cout << "Key = " << *name_it << " and Value = " << class_dict.getValue<str>(*name_it) << std::endl;
          str key, value;
          key = *name_it;
          value = class_dict.getValue<str>(*name_it);
          sprintf(cosmo.fc.name[i],"%s",key.c_str());
          sprintf(cosmo.fc.value[i],"%s",value.c_str());
          i++;
        }
      }

      strcpy(cosmo.fc.name[1],"omega_b");
      strcpy(cosmo.fc.name[2],"omega_cdm");
      strcpy(cosmo.fc.name[3],"H0");
      strcpy(cosmo.fc.name[4],"ln10^{10}A_s");
      strcpy(cosmo.fc.name[5],"n_s");
      strcpy(cosmo.fc.name[6],"tau_reio");
      strcpy(cosmo.fc.name[7],"r");

      sprintf(cosmo.fc.value[1],"%e",*Param["omega_b"]);
      sprintf(cosmo.fc.value[2],"%e",*Param["omega_cdm"]);
      sprintf(cosmo.fc.value[3],"%e",*Param["H0"]);
      sprintf(cosmo.fc.value[4],"%e",*Param["ln10A_s"]);
      sprintf(cosmo.fc.value[5],"%e",*Param["n_s"]);
      sprintf(cosmo.fc.value[6],"%e",*Param["tau_reio"]);
      sprintf(cosmo.fc.value[7],"%e",*Param["r_tensor"]);
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

      // What follows is a loop to read all generic additional parameters for CLASS
      // needed to (hopefully) reproduce the LCDM values of Planck, which we pass
      // as options a YAML-node named "class_dict" within the 'Rules' section of the input,
      // and which will be added if there are present.

      int len_of_input = 11;
      Options class_dict;
      std::vector<str> names;

      if (runOptions->hasKey("class_dict"))
      {
        class_dict = Options(runOptions->getValue<YAML::Node>("class_dict"));
        names = class_dict.getNames();
        len_of_input += names.size();
      }

      BEreq::class_parser_initialize(&cosmo.fc,byVal(len_of_input),"",cosmo.class_errmsg);

      strcpy(cosmo.fc.name[0],"output");
      strcpy(cosmo.fc.value[0],"tCl pCl lCl");
      strcpy(cosmo.fc.name[8],"l_max_scalars");
      sprintf(cosmo.fc.value[8],"%d",l_max);
      strcpy(cosmo.fc.name[9],"modes");
      strcpy(cosmo.fc.value[9],"s,t");
      strcpy(cosmo.fc.name[10],"lensing");
      strcpy(cosmo.fc.value[10],"yes");

      int i = 11;

      if (runOptions->hasKey("class_dict"))
      {
        for (std::vector<str>::iterator name_it = names.begin(); name_it != names.end(); name_it++)
        {
          //std::cout << "Key = " << *name_it << " and Value = " << class_dict.getValue<str>(*name_it) << std::endl;
          str key, value;
          key = *name_it;
          value = class_dict.getValue<str>(*name_it);
          sprintf(cosmo.fc.name[i],"%s",key.c_str());
          sprintf(cosmo.fc.value[i],"%s",value.c_str());
          i++;
        }
      }

      // Initialization parameters controlling main characteristics.
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");

      // Control the output of analytic approximations for comparison.
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");

      // Parameters to control how the ICs are sampled.
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");

      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this

      // Parameters to control how the vparams are sampled.
      int param_sampling = runOptions->getValue<int> ("param_sampling");

      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");

      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");

      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");

      // Parameters to be passed to the potential
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

      //-------------------------------------------------
      //----------------- Sampling m^2 ------------------
      //-------------------------------------------------
      vparams[0] = *Param["m2_inflaton"];

      //-------------------------------------------------
      //--------------- Sampling N_pivot ----------------
      //-------------------------------------------------
      double N_pivot = *Param["N_pivot"];

      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");

      // Priors on the IC and N_pivot ranges
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");

      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");

      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");

      std::cout << "calc_full_pk = " << calc_full_pk << std::endl;

      calc_full_pk = 0; // lets hack for now as I decide how to implement full spectra to Patrick's version.

      std::cout << "steps = " << steps << std::endl;

      gambit_inflation_observables observs;

      std::cout << "num_inflaton = " << num_inflaton << std::endl;
      std::cout << "save_iso_N = " << save_iso_N << std::endl;
      std::cout << "N_pivot = " << N_pivot << std::endl;
      std::cout << "kmin = " << kmin << std::endl;
      std::cout << "dlnk = " << dlnk << std::endl;
      std::cout << "lambda = " << vparams[0] << std::endl;

      // The function below calls the multimodecode backend for a given choice of inflationary model,
      // which calculates the observables.
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
           byVal(&dphi_init0[0]), // dphi correct at the multimodecode
           byVal(&vparams[0]),
           N_pivot,
           k_pivot,
           dlnk,
           calc_full_pk,
           steps,
           kmin,
           //kmax,
           byVal(&phi0_priors_min[0]),
           byVal(&phi0_priors_max[0]),
           byVal(&dphi0_priors_min[0]),
           byVal(&dphi0_priors_max[0]),
           N_pivot_prior_min,
           N_pivot_prior_max);

      std::cout << " we are out of the multimodecode_gambit_driver function" << std::endl;

      std::cout << "observs.k_array = " << sizeof(observs.k_array) << std::endl;
      std::cout << "observs.k_array[0] = " << observs.k_array[0] << std::endl;
      //std::cout << "observs.pks_array = " << observs.pks_array << std::endl;
      //std::cout << "observs.pkt_array = " << observs.pkt_array << std::endl;

      std::cout << "calc_full_pk = " << calc_full_pk << std::endl;
      if (calc_full_pk == 0)
      {
        std::cout << "observs.As " << observs.As << std::endl;
        std::cout << "observs.A_iso " << observs.A_iso << std::endl;
        std::cout << "observs.A_pnad " << observs.A_pnad << std::endl;
        std::cout << "observs.A_ent " << observs.A_ent << std::endl;
        std::cout << "observs.A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
        //std::cout << "observs.A_bundle " << observs.A_bundle << std::endl;
        std::cout << "observs.ns " << observs.ns << std::endl;
        std::cout << "observs.nt " << observs.nt << std::endl;
        std::cout << "observs.n_iso " << observs.n_iso << std::endl;
        std::cout << "observs.n_pnad " << observs.n_pnad << std::endl;
        std::cout << "observs.n_ent " << observs.n_ent << std::endl;
        std::cout << "observs.r " << observs.r << std::endl;
        std::cout << "observs.alpha_s " << observs.alpha_s << std::endl;
        std::cout << "observs.runofrun " << observs.runofrun << std::endl;
        std::cout << "observs.f_NL " << observs.f_NL << std::endl;
        std::cout << "observs.tau_NL " << observs.tau_NL << std::endl;

        strcpy(cosmo.fc.name[1],"omega_b");
        strcpy(cosmo.fc.name[2],"omega_cdm");
        strcpy(cosmo.fc.name[3],"H0");
        strcpy(cosmo.fc.name[4],"ln10^{10}A_s");
        strcpy(cosmo.fc.name[5],"n_s");
        strcpy(cosmo.fc.name[6],"tau_reio");
        strcpy(cosmo.fc.name[7],"r");

        sprintf(cosmo.fc.value[1],"%e",*Param["omega_b"]);
        sprintf(cosmo.fc.value[2],"%e",*Param["omega_cdm"]);
        sprintf(cosmo.fc.value[3],"%e",*Param["H0"]);
        sprintf(cosmo.fc.value[4],"%e",observs.As);
        sprintf(cosmo.fc.value[5],"%e",observs.ns);
        sprintf(cosmo.fc.value[6],"%e",*Param["tau_reio"]);
        sprintf(cosmo.fc.value[7],"%e",observs.r);
      }
    }

    void class_set_parameter_inf_1quarInf_LCDMt(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_inf_1quarInf_LCDMt" << std::endl;
      using namespace Pipes::class_set_parameter_inf_1quarInf_LCDMt;

      int l_max=cosmo.lmax;

      // What follows is a loop to read all generic additional parameters for CLASS
      // needed to (hopefully) reproduce the LCDM values of Planck, which we pass
      // as options a YAML-node named "class_dict" within the 'Rules' section of the input,
      // and which will be added if there are present.

      int len_of_input = 11;
      Options class_dict;
      std::vector<str> names;

      if (runOptions->hasKey("class_dict"))
      {
  class_dict = Options(runOptions->getValue<YAML::Node>("class_dict"));
  names = class_dict.getNames();
  len_of_input += names.size();
      }

      BEreq::class_parser_initialize(&cosmo.fc,byVal(len_of_input),"",cosmo.class_errmsg);

      strcpy(cosmo.fc.name[0],"output");
      strcpy(cosmo.fc.value[0],"tCl pCl lCl");
      strcpy(cosmo.fc.name[8],"l_max_scalars");
      sprintf(cosmo.fc.value[8],"%d",l_max);
      strcpy(cosmo.fc.name[9],"modes");
      strcpy(cosmo.fc.value[9],"s,t");
      strcpy(cosmo.fc.name[10],"lensing");
      strcpy(cosmo.fc.value[10],"yes");

      int i = 11;

      if (runOptions->hasKey("class_dict"))
      {
        for (std::vector<str>::iterator name_it = names.begin(); name_it != names.end(); name_it++)
        {
          //std::cout << "Key = " << *name_it << " and Value = " << class_dict.getValue<str>(*name_it) << std::endl;
          str key, value;
          key = *name_it;
          value = class_dict.getValue<str>(*name_it);
          sprintf(cosmo.fc.name[i],"%s",key.c_str());
          sprintf(cosmo.fc.value[i],"%s",value.c_str());
          i++;
        }
      }

      // Initialization parameters controlling main characteristics.
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");

      // Control the output of analytic approximations for comparison.
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");

      // Parameters to control how the ICs are sampled.
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");

      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this

      // Parameters to control how the vparams are sampled.
      int param_sampling = runOptions->getValue<int> ("param_sampling");

      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");

      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");

      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");

      // Parameters to be passed to the potential
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

      //-------------------------------------------------
      //----------------- Sampling \lambda ------------------
      //-------------------------------------------------
      vparams[0] = *Param["lambda"];

      //-------------------------------------------------
      //--------------- Sampling N_pivot ----------------
      //-------------------------------------------------
      double N_pivot = *Param["N_pivot"];

      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");

      // Priors on the IC and N_pivot ranges
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");

      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");

      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");

      std::cout << "calc_full_pk = " << calc_full_pk << std::endl;

      calc_full_pk = 0; // lets hack for now as I decide how to implement full spectra to Patrick's version.

      std::cout << "steps = " << steps << std::endl;

      gambit_inflation_observables observs;

      std::cout << "num_inflaton = " << num_inflaton << std::endl;
      std::cout << "save_iso_N = " << save_iso_N << std::endl;
      std::cout << "N_pivot = " << N_pivot << std::endl;
      std::cout << "kmin = " << kmin << std::endl;
      std::cout << "dlnk = " << dlnk << std::endl;
      std::cout << "lambda = " << vparams[0] << std::endl;

      // The function below calls the multimodecode backend for a given choice of inflationary model,
      // which calculates the observables.
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
           byVal(&dphi_init0[0]), // dphi correct at the multimodecode
           byVal(&vparams[0]),
           N_pivot,
           k_pivot,
           dlnk,
           calc_full_pk,
           steps,
           kmin,
           //kmax,
           byVal(&phi0_priors_min[0]),
           byVal(&phi0_priors_max[0]),
           byVal(&dphi0_priors_min[0]),
           byVal(&dphi0_priors_max[0]),
           N_pivot_prior_min,
           N_pivot_prior_max);

      std::cout << " we are out of the multimodecode_gambit_driver function" << std::endl;

      std::cout << "observs.k_array = " << sizeof(observs.k_array) << std::endl;
      std::cout << "observs.k_array[0] = " << observs.k_array[0] << std::endl;
      //std::cout << "observs.pks_array = " << observs.pks_array << std::endl;
      //std::cout << "observs.pkt_array = " << observs.pkt_array << std::endl;

      std::cout << "calc_full_pk = " << calc_full_pk << std::endl;
      if (calc_full_pk == 0)
      {
        std::cout << "observs.As " << observs.As << std::endl;
        std::cout << "observs.A_iso " << observs.A_iso << std::endl;
        std::cout << "observs.A_pnad " << observs.A_pnad << std::endl;
        std::cout << "observs.A_ent " << observs.A_ent << std::endl;
        std::cout << "observs.A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
        //std::cout << "observs.A_bundle " << observs.A_bundle << std::endl;
        std::cout << "observs.ns " << observs.ns << std::endl;
        std::cout << "observs.nt " << observs.nt << std::endl;
        std::cout << "observs.n_iso " << observs.n_iso << std::endl;
        std::cout << "observs.n_pnad " << observs.n_pnad << std::endl;
        std::cout << "observs.n_ent " << observs.n_ent << std::endl;
        std::cout << "observs.r " << observs.r << std::endl;
        std::cout << "observs.alpha_s " << observs.alpha_s << std::endl;
        std::cout << "observs.runofrun " << observs.runofrun << std::endl;
        std::cout << "observs.f_NL " << observs.f_NL << std::endl;
        std::cout << "observs.tau_NL " << observs.tau_NL << std::endl;

        strcpy(cosmo.fc.name[1],"omega_b");
        strcpy(cosmo.fc.name[2],"omega_cdm");
        strcpy(cosmo.fc.name[3],"H0");
        strcpy(cosmo.fc.name[4],"ln10^{10}A_s");
        strcpy(cosmo.fc.name[5],"n_s");
        strcpy(cosmo.fc.name[6],"tau_reio");
        strcpy(cosmo.fc.name[7],"r");

        sprintf(cosmo.fc.value[1],"%e",*Param["omega_b"]);
        sprintf(cosmo.fc.value[2],"%e",*Param["omega_cdm"]);
        sprintf(cosmo.fc.value[3],"%e",*Param["H0"]);
        sprintf(cosmo.fc.value[4],"%e",observs.As);
        sprintf(cosmo.fc.value[5],"%e",observs.ns);
        sprintf(cosmo.fc.value[6],"%e",*Param["tau_reio"]);
        sprintf(cosmo.fc.value[7],"%e",observs.r);
      }
    }

    void class_set_parameter_inf_1mono32Inf_LCDMt(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_inf_1mono32Inf_LCDMt" << std::endl;
      using namespace Pipes::class_set_parameter_inf_1mono32Inf_LCDMt;

      int l_max=cosmo.lmax;

      // What follows is a loop to read all generic additional parameters for CLASS
      // needed to (hopefully) reproduce the LCDM values of Planck, which we pass
      // as options a YAML-node named "class_dict" within the 'Rules' section of the input,
      // and which will be added if there are present.

      int len_of_input = 11;
      Options class_dict;
      std::vector<str> names;

      if (runOptions->hasKey("class_dict"))
      {
        class_dict = Options(runOptions->getValue<YAML::Node>("class_dict"));
        names = class_dict.getNames();
        len_of_input += names.size();
      }

      BEreq::class_parser_initialize(&cosmo.fc,byVal(len_of_input),"",cosmo.class_errmsg);

      strcpy(cosmo.fc.name[0],"output");
      strcpy(cosmo.fc.value[0],"tCl pCl lCl");
      strcpy(cosmo.fc.name[8],"l_max_scalars");
      sprintf(cosmo.fc.value[8],"%d",l_max);
      strcpy(cosmo.fc.name[9],"modes");
      strcpy(cosmo.fc.value[9],"s,t");
      strcpy(cosmo.fc.name[10],"lensing");
      strcpy(cosmo.fc.value[10],"yes");

      int i = 11;

      if (runOptions->hasKey("class_dict"))
      {
        for (std::vector<str>::iterator name_it = names.begin(); name_it != names.end(); name_it++)
        {
          //std::cout << "Key = " << *name_it << " and Value = " << class_dict.getValue<str>(*name_it) << std::endl;
          str key, value;
          key = *name_it;
          value = class_dict.getValue<str>(*name_it);
          sprintf(cosmo.fc.name[i],"%s",key.c_str());
          sprintf(cosmo.fc.value[i],"%s",value.c_str());
          i++;
        }
      }

      // Initialization parameters controlling main characteristics.
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");

      // Control the output of analytic approximations for comparison.
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");

      // Parameters to control how the ICs are sampled.
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");

      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this

      // Parameters to control how the vparams are sampled.
      int param_sampling = runOptions->getValue<int> ("param_sampling");

      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");

      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");

      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");

      // Parameters to be passed to the potential
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

      //-------------------------------------------------
      //----------------- Sampling \lambda ------------------
      //-------------------------------------------------
      vparams[0] = *Param["lambda"];

      //-------------------------------------------------
      //--------------- Sampling N_pivot ----------------
      //-------------------------------------------------
      double N_pivot = *Param["N_pivot"];

      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");

      // Priors on the IC and N_pivot ranges
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");

      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");

      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");

      std::cout << "calc_full_pk = " << calc_full_pk << std::endl;

      calc_full_pk = 0; // lets hack for now as I decide how to implement full spectra to Patrick's version.

      std::cout << "steps = " << steps << std::endl;

      gambit_inflation_observables observs;

      std::cout << "num_inflaton = " << num_inflaton << std::endl;
      std::cout << "save_iso_N = " << save_iso_N << std::endl;
      std::cout << "N_pivot = " << N_pivot << std::endl;
      std::cout << "kmin = " << kmin << std::endl;
      std::cout << "dlnk = " << dlnk << std::endl;
      std::cout << "lambda = " << vparams[0] << std::endl;

      // The function below calls the multimodecode backend for a given choice of inflationary model,
      // which calculates the observables.
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
           byVal(&dphi_init0[0]), // dphi correct at the multimodecode
           byVal(&vparams[0]),
           N_pivot,
           k_pivot,
           dlnk,
           calc_full_pk,
           steps,
           kmin,
           //kmax,
           byVal(&phi0_priors_min[0]),
           byVal(&phi0_priors_max[0]),
           byVal(&dphi0_priors_min[0]),
           byVal(&dphi0_priors_max[0]),
           N_pivot_prior_min,
           N_pivot_prior_max);


      std::cout << " we are out of the multimodecode_gambit_driver function" << std::endl;

      std::cout << "observs.k_array = " << sizeof(observs.k_array) << std::endl;
      std::cout << "observs.k_array[0] = " << observs.k_array[0] << std::endl;
      //std::cout << "observs.pks_array = " << observs.pks_array << std::endl;
      //std::cout << "observs.pkt_array = " << observs.pkt_array << std::endl;

      std::cout << "calc_full_pk = " << calc_full_pk << std::endl;
      if (calc_full_pk == 0)
      {
        std::cout << "observs.As " << observs.As << std::endl;
        std::cout << "observs.A_iso " << observs.A_iso << std::endl;
        std::cout << "observs.A_pnad " << observs.A_pnad << std::endl;
        std::cout << "observs.A_ent " << observs.A_ent << std::endl;
        std::cout << "observs.A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
        //std::cout << "observs.A_bundle " << observs.A_bundle << std::endl;
        std::cout << "observs.ns " << observs.ns << std::endl;
        std::cout << "observs.nt " << observs.nt << std::endl;
        std::cout << "observs.n_iso " << observs.n_iso << std::endl;
        std::cout << "observs.n_pnad " << observs.n_pnad << std::endl;
        std::cout << "observs.n_ent " << observs.n_ent << std::endl;
        std::cout << "observs.r " << observs.r << std::endl;
        std::cout << "observs.alpha_s " << observs.alpha_s << std::endl;
        std::cout << "observs.runofrun " << observs.runofrun << std::endl;
        std::cout << "observs.f_NL " << observs.f_NL << std::endl;
        std::cout << "observs.tau_NL " << observs.tau_NL << std::endl;

        strcpy(cosmo.fc.name[1],"omega_b");
        strcpy(cosmo.fc.name[2],"omega_cdm");
        strcpy(cosmo.fc.name[3],"H0");
        strcpy(cosmo.fc.name[4],"ln10^{10}A_s");
        strcpy(cosmo.fc.name[5],"n_s");
        strcpy(cosmo.fc.name[6],"tau_reio");
        strcpy(cosmo.fc.name[7],"r");

        sprintf(cosmo.fc.value[1],"%e",*Param["omega_b"]);
        sprintf(cosmo.fc.value[2],"%e",*Param["omega_cdm"]);
        sprintf(cosmo.fc.value[3],"%e",*Param["H0"]);
        sprintf(cosmo.fc.value[4],"%e",observs.As);
        sprintf(cosmo.fc.value[5],"%e",observs.ns);
        sprintf(cosmo.fc.value[6],"%e",*Param["tau_reio"]);
        sprintf(cosmo.fc.value[7],"%e",observs.r);
      }
    }

    void class_set_parameter_inf_1linearInf_LCDMt(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_inf_1linearInf_LCDMt" << std::endl;
      using namespace Pipes::class_set_parameter_inf_1linearInf_LCDMt;

      int l_max=cosmo.lmax;

      // What follows is a loop to read all generic additional parameters for CLASS
      // needed to (hopefully) reproduce the LCDM values of Planck, which we pass
      // as options a YAML-node named "class_dict" within the 'Rules' section of the input,
      // and which will be added if there are present.

      int len_of_input = 11;
      Options class_dict;
      std::vector<str> names;

      if (runOptions->hasKey("class_dict"))
      {
        class_dict = Options(runOptions->getValue<YAML::Node>("class_dict"));
        names = class_dict.getNames();
        len_of_input += names.size();
      }

      BEreq::class_parser_initialize(&cosmo.fc,byVal(len_of_input),"",cosmo.class_errmsg);

      strcpy(cosmo.fc.name[0],"output");
      strcpy(cosmo.fc.value[0],"tCl pCl lCl");
      strcpy(cosmo.fc.name[8],"l_max_scalars");
      sprintf(cosmo.fc.value[8],"%d",l_max);
      strcpy(cosmo.fc.name[9],"modes");
      strcpy(cosmo.fc.value[9],"s,t");
      strcpy(cosmo.fc.name[10],"lensing");
      strcpy(cosmo.fc.value[10],"yes");

      int i = 11;

      if (runOptions->hasKey("class_dict"))
      {
        for (std::vector<str>::iterator name_it = names.begin(); name_it != names.end(); name_it++)
        {
          //std::cout << "Key = " << *name_it << " and Value = " << class_dict.getValue<str>(*name_it) << std::endl;
          str key, value;
          key = *name_it;
          value = class_dict.getValue<str>(*name_it);
          sprintf(cosmo.fc.name[i],"%s",key.c_str());
          sprintf(cosmo.fc.value[i],"%s",value.c_str());
          i++;
        }
      }

      // Initialization parameters controlling main characteristics.
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");

      // Control the output of analytic approximations for comparison.
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");

      // Parameters to control how the ICs are sampled.
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");

      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this

      // Parameters to control how the vparams are sampled.
      int param_sampling = runOptions->getValue<int> ("param_sampling");

      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");

      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");

      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");

      // Parameters to be passed to the potential
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

      //-------------------------------------------------
      //----------------- Sampling \lambda ------------------
      //-------------------------------------------------
      vparams[0] = *Param["lambda"];

      //-------------------------------------------------
      //--------------- Sampling N_pivot ----------------
      //-------------------------------------------------
      double N_pivot = *Param["N_pivot"];

      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");

      // Priors on the IC and N_pivot ranges
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");

      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");

      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");

      std::cout << "calc_full_pk = " << calc_full_pk << std::endl;

      calc_full_pk = 0; // lets hack for now as I decide how to implement full spectra to Patrick's version.

      std::cout << "steps = " << steps << std::endl;

      gambit_inflation_observables observs;

      std::cout << "num_inflaton = " << num_inflaton << std::endl;
      std::cout << "save_iso_N = " << save_iso_N << std::endl;
      std::cout << "N_pivot = " << N_pivot << std::endl;
      std::cout << "kmin = " << kmin << std::endl;
      std::cout << "dlnk = " << dlnk << std::endl;
      std::cout << "lambda = " << vparams[0] << std::endl;

      // The function below calls the multimodecode backend for a given choice of inflationary model,
      // which calculates the observables.
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
           byVal(&dphi_init0[0]), // dphi correct at the multimodecode
           byVal(&vparams[0]),
           N_pivot,
           k_pivot,
           dlnk,
           calc_full_pk,
           steps,
           kmin,
           //kmax,
           byVal(&phi0_priors_min[0]),
           byVal(&phi0_priors_max[0]),
           byVal(&dphi0_priors_min[0]),
           byVal(&dphi0_priors_max[0]),
           N_pivot_prior_min,
           N_pivot_prior_max);

      std::cout << " we are out of the multimodecode_gambit_driver function" << std::endl;

      std::cout << "observs.k_array = " << sizeof(observs.k_array) << std::endl;
      std::cout << "observs.k_array[0] = " << observs.k_array[0] << std::endl;
      //std::cout << "observs.pks_array = " << observs.pks_array << std::endl;
      //std::cout << "observs.pkt_array = " << observs.pkt_array << std::endl;

      std::cout << "calc_full_pk = " << calc_full_pk << std::endl;
      if (calc_full_pk == 0)
      {
        std::cout << "observs.As " << observs.As << std::endl;
        std::cout << "observs.A_iso " << observs.A_iso << std::endl;
        std::cout << "observs.A_pnad " << observs.A_pnad << std::endl;
        std::cout << "observs.A_ent " << observs.A_ent << std::endl;
        std::cout << "observs.A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
        //std::cout << "observs.A_bundle " << observs.A_bundle << std::endl;
        std::cout << "observs.ns " << observs.ns << std::endl;
        std::cout << "observs.nt " << observs.nt << std::endl;
        std::cout << "observs.n_iso " << observs.n_iso << std::endl;
        std::cout << "observs.n_pnad " << observs.n_pnad << std::endl;
        std::cout << "observs.n_ent " << observs.n_ent << std::endl;
        std::cout << "observs.r " << observs.r << std::endl;
        std::cout << "observs.alpha_s " << observs.alpha_s << std::endl;
        std::cout << "observs.runofrun " << observs.runofrun << std::endl;
        std::cout << "observs.f_NL " << observs.f_NL << std::endl;
        std::cout << "observs.tau_NL " << observs.tau_NL << std::endl;


        strcpy(cosmo.fc.name[1],"omega_b");
        strcpy(cosmo.fc.name[2],"omega_cdm");
        strcpy(cosmo.fc.name[3],"H0");
        strcpy(cosmo.fc.name[4],"ln10^{10}A_s");
        strcpy(cosmo.fc.name[5],"n_s");
        strcpy(cosmo.fc.name[6],"tau_reio");
        strcpy(cosmo.fc.name[7],"r");

        sprintf(cosmo.fc.value[1],"%e",*Param["omega_b"]);
        sprintf(cosmo.fc.value[2],"%e",*Param["omega_cdm"]);
        sprintf(cosmo.fc.value[3],"%e",*Param["H0"]);
        sprintf(cosmo.fc.value[4],"%e",observs.As);
        sprintf(cosmo.fc.value[5],"%e",observs.ns);
        sprintf(cosmo.fc.value[6],"%e",*Param["tau_reio"]);
        sprintf(cosmo.fc.value[7],"%e",observs.r);
      }
    }

    void class_set_parameter_inf_smashInf_LCDMt(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_set_parameter_inf_smashInf_LCDMt" << std::endl;
      using namespace Pipes::class_set_parameter_inf_smashInf_LCDMt;

      int l_max=cosmo.lmax;

      // What follows is a loop to read all generic additional parameters for CLASS
      // needed to (hopefully) reproduce the LCDM values of Planck, which we pass
      // as options a YAML-node named "class_dict" within the 'Rules' section of the input,
      // and which will be added if there are present.

      int len_of_input = 11;
      Options class_dict;
      std::vector<str> names;

      if (runOptions->hasKey("class_dict"))
      {
        class_dict = Options(runOptions->getValue<YAML::Node>("class_dict"));
        names = class_dict.getNames();
        len_of_input += names.size();
      }

      BEreq::class_parser_initialize(&cosmo.fc,byVal(len_of_input),"",cosmo.class_errmsg);

      strcpy(cosmo.fc.name[0],"output");
      strcpy(cosmo.fc.value[0],"tCl pCl lCl");
      strcpy(cosmo.fc.name[8],"l_max_scalars");
      sprintf(cosmo.fc.value[8],"%d",l_max);
      strcpy(cosmo.fc.name[9],"modes");
      strcpy(cosmo.fc.value[9],"s,t");
      strcpy(cosmo.fc.name[10],"lensing");
      strcpy(cosmo.fc.value[10],"yes");

      int i = 11;

      if (runOptions->hasKey("class_dict"))
      {
        for (std::vector<str>::iterator name_it = names.begin(); name_it != names.end(); name_it++)
        {
          //std::cout << "Key = " << *name_it << " and Value = " << class_dict.getValue<str>(*name_it) << std::endl;
          str key, value;
          key = *name_it;
          value = class_dict.getValue<str>(*name_it);
          sprintf(cosmo.fc.name[i],"%s",key.c_str());
          sprintf(cosmo.fc.value[i],"%s",value.c_str());
          i++;
        }
      }

      /*
        This part needs cleaning up and commenting - will get on to in asap.
      */
      // Parameters to be passed to the potential
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

      //-------------------------------------------------
      //----------- at the moment I override-------------
      //-------------------------------------------------
      vparams[0] = *Param["log10_xi"];
      vparams[1] = *Param["log10_beta"];
      vparams[2] = *Param["log10_lambda"];

      std::cout << "log10[xi] = " << vparams[0] << std::endl;
      std::cout << "log10[beta] = " << vparams[1] << std::endl;
      std::cout << "log10[lambda] = " << vparams[2] << std::endl;

      /* coefficients of P(x) =  -8 + b*\[Phi]**2  */
      /*                            + b*\[Xi]*\[Phi]**4
            + 6*\[Xi]**2*\[Phi]**4 */
      double smashp1[5] = { -8.0, 0, pow(10.0,vparams[1]), 0, (pow(10.0,vparams[1])*
                     pow(10.0,vparams[0])+
                     6.0*pow(10.0,vparams[0])*
                     pow(10.0,vparams[0]))};

      double smashd1[8];

      gsl_poly_complex_workspace * wsmash = gsl_poly_complex_workspace_alloc (5);
      gsl_poly_complex_solve (smashp1, 5, wsmash, smashd1);
      gsl_poly_complex_workspace_free (wsmash);

      for (i = 0; i < 4; i++)
      {
        printf ("z%d = %+.18f %+.18f\n", i, smashd1[2*i], smashd1[2*i+1]);
      }

      double badway[4] = {smashd1[0],smashd1[2],smashd1[4],smashd1[6]};
      double tempbw = 0;

      for(int i=0;i<4;i++)
      {
        if(badway[i]>tempbw) tempbw=badway[i];
      }
      //gsl_sf_log((1.0+pow(10.0,vparams[0])*smashp1[4]*smashp1[4])/
      //       (1.0+pow(10.0,vparams[0])*)

      vparams[3] = tempbw;

      std::cout << "phi_end = " << tempbw << std::endl;

      //-------------------------------------------------
      //--------------- Sampling N_pivot ----------------
      //-------------------------------------------------
      double N_pivot = *Param["N_pivot"];

      ///////////// solving the equality to calculate the
      ////////////  slow roll field amplitute at N_pivot.
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

      printf ("using %s method\n",gsl_min_fminimizer_name (s));

      printf ("%5s [%9s, %9s] %9s %10s %9s\n","iter", "lower", "upper", "min","err", "err(est)");

      printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",iter, a, b, m, m - m_expected, b - a);

      do
      {
        iter++;
        status = gsl_min_fminimizer_iterate (s);

        m = gsl_min_fminimizer_x_minimum (s);
        a = gsl_min_fminimizer_x_lower (s);
        b = gsl_min_fminimizer_x_upper (s);

        status = gsl_min_test_interval (a, b, 0.01, 0.0);

        if (status == GSL_SUCCESS) printf ("Converged:\n");

        printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
                  iter, a, b,
                  m, m - m_expected, b - a);
      }
      while (status == GSL_CONTINUE && iter < max_iter);

      gsl_min_fminimizer_free (s);

      //phi_init0[0] = a*1.01;

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
      r_self = r_SR(smash_eps,smash_eta);

      /* Have calculated the field value at N_pivot
        predicted by slow-roll approximation.
        Choosing a very close slightly higher
        value will be satisfactory for initial conditions.
      */

      strcpy(cosmo.fc.name[1],"omega_b");
      strcpy(cosmo.fc.name[2],"omega_cdm");
      strcpy(cosmo.fc.name[3],"H0");
      strcpy(cosmo.fc.name[4],"ln10^{10}A_s");
      strcpy(cosmo.fc.name[5],"n_s");
      strcpy(cosmo.fc.name[6],"tau_reio");
      strcpy(cosmo.fc.name[7],"r");

      sprintf(cosmo.fc.value[1],"%e",*Param["omega_b"]);
      sprintf(cosmo.fc.value[2],"%e",*Param["omega_cdm"]);
      sprintf(cosmo.fc.value[3],"%e",*Param["H0"]);
      sprintf(cosmo.fc.value[4],"%e",As_self);
      sprintf(cosmo.fc.value[5],"%e",ns_self);
      sprintf(cosmo.fc.value[6],"%e",*Param["tau_reio"]);
      sprintf(cosmo.fc.value[7],"%e",r_self);
    }

    void class_run_func(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_run_func" << std::endl;
      using namespace Pipes::class_run_func;

      char error_printout[1024];
      cosmo = *Dep::class_set_parameter;

      if (BEreq::class_input_initialize(&cosmo.fc,&cosmo.pr,&cosmo.ba,&cosmo.th,&cosmo.pt,&cosmo.tr,&cosmo.pm,&cosmo.sp,&cosmo.nl,&cosmo.le,&cosmo.op,cosmo.class_errmsg) == _FAILURE_)
      {
        sprintf(error_printout,"Error in class_input_initialize\n=>%s\n",cosmo.class_errmsg);
        invalid_point().raise(error_printout);
      }
      if (BEreq::class_background_initialize(&cosmo.pr,&cosmo.ba) == _FAILURE_)
      {
        sprintf(error_printout,"Error in class_background_initialize\n=>%s\n",cosmo.ba.error_message);
        invalid_point().raise(error_printout);
      }
      if (BEreq::class_thermodynamics_initialize(&cosmo.pr,&cosmo.ba,&cosmo.th) == _FAILURE_)
      {
        sprintf(error_printout,"Error in class_thermodynamics_initialize\n=>%s\n",cosmo.th.error_message);
        invalid_point().raise(error_printout);
      }
      if (BEreq::class_perturb_initialize(&cosmo.pr,&cosmo.ba,&cosmo.th,&cosmo.pt) == _FAILURE_)
      {
        sprintf(error_printout,"Error in class_perturb_initialize\n=>%s\n",cosmo.pt.error_message);
        invalid_point().raise(error_printout);
      }
      if (BEreq::class_primordial_initialize(&cosmo.pr,&cosmo.pt,&cosmo.pm) == _FAILURE_)
      {
        sprintf(error_printout,"Error in class_primordial_initialize\n=>%s\n",cosmo.pm.error_message);
        invalid_point().raise(error_printout);
      }
      if (BEreq::class_nonlinear_initialize(&cosmo.pr,&cosmo.ba,&cosmo.th,&cosmo.pt,&cosmo.pm,&cosmo.nl) == _FAILURE_)
      {
        sprintf(error_printout,"Error in class_nonlinear_initialize\n=>%s\n",cosmo.nl.error_message);
        invalid_point().raise(error_printout);
      }
      if (BEreq::class_transfer_initialize(&cosmo.pr,&cosmo.ba,&cosmo.th,&cosmo.pt,&cosmo.nl,&cosmo.tr) == _FAILURE_)
      {
        sprintf(error_printout,"Error in class_transfer_initialize\n=>%s\n",cosmo.tr.error_message);
        invalid_point().raise(error_printout);
      }
      if (BEreq::class_spectra_initialize(&cosmo.pr,&cosmo.ba,&cosmo.pt,&cosmo.pm,&cosmo.nl,&cosmo.tr,&cosmo.sp) == _FAILURE_)
      {
        sprintf(error_printout,"Error in class_spectra_initialize\n=>%s\n",cosmo.sp.error_message);
        invalid_point().raise(error_printout);
      }
      if (BEreq::class_lensing_initialize(&cosmo.pr,&cosmo.pt,&cosmo.sp,&cosmo.nl,&cosmo.le) == _FAILURE_)
      {
        sprintf(error_printout,"Error in class_lensing_initialize\n=>%s\n",cosmo.le.error_message);
        invalid_point().raise(error_printout);
      }
      cosmo.non_free_pointer = true;
    }

    void class_get_spectra_func(Class_container& cosmo)
    {
      //std::cout << "Last seen alive in: class_get_spectra_func" << std::endl;
      using namespace Pipes::class_get_spectra_func;

      cosmo = *Dep::class_run;

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

      // Now that all calculations with CLASS are done, free the pointers which were allocated in the meantime.
      if (cosmo.non_free_pointer)
      {
        BEreq::class_lensing_free(&cosmo.le);
        BEreq::class_spectra_free(&cosmo.sp);
        BEreq::class_transfer_free(&cosmo.tr);
        BEreq::class_nonlinear_free(&cosmo.nl);
        BEreq::class_primordial_free(&cosmo.pm);
        BEreq::class_perturb_free(&cosmo.pt);
        BEreq::class_thermodynamics_free(&cosmo.th);
        BEreq::class_background_free(&cosmo.ba);
        cosmo.non_free_pointer = false;
      }
    }

    double** return_vanilla_cls(double omega_b,double omega_cdm,double H0,double ln10A_s,double n_s,double tau_reio)
    {
      using namespace Pipes::function_vanilla_lowp_TT_loglike;

      struct Class::precision pr;        /* for precision parameters */
      struct Class::background ba;       /* for cosmological background */
      struct Class::thermo th;           /* for thermodynamics */
      struct Class::perturbs pt;         /* for source functions */
      struct Class::transfers tr;        /* for transfer functions */
      struct Class::primordial pm;       /* for primordial spectra */
      struct Class::spectra sp;          /* for output spectra */
      struct Class::nonlinear nl;        /* for non-linear spectra */
      struct Class::lensing le;          /* for lensed spectra */
      struct Class::output op;           /* for output files */
      Class::ErrorMsg class_errmsg;      /* for error messages */

      int l,l_max;
      int num_ct_max=7;

      l_max=3000;

      char *class_null;// defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
       could add another array with P(k), or extract other results from
       the code - here we assume that we are interested in the C_l's
       only */

      double* cl[l_max];
      for(int i = 0; i < l_max; ++i)
        cl[i] = new double[num_ct_max];

      struct Class::file_content fc;

      BEreq::class_parser_initialize(&fc,9,"",class_errmsg);

      strcpy(fc.name[0],"output");
      strcpy(fc.value[0],"tCl,pCl");

      strcpy(fc.name[1],"omega_b");
      strcpy(fc.name[2],"omega_cdm");
      strcpy(fc.name[3],"H0");
      strcpy(fc.name[4],"ln10^{10}A_s");
      strcpy(fc.name[5],"n_s");
      strcpy(fc.name[6],"tau_reio");

      //strcpy(fc.name[7],"l_scalar_max");
      strcpy(fc.name[7],"l_max_scalars");
      sprintf(fc.value[7],"%d",3000);

      sprintf(fc.value[1],"%e",omega_b);
      sprintf(fc.value[2],"%e",omega_cdm);
      sprintf(fc.value[3],"%e",H0);
      sprintf(fc.value[4],"%e",ln10A_s);
      sprintf(fc.value[5],"%e",n_s);
      sprintf(fc.value[6],"%e",tau_reio);

      std::cout << "omega_b = " << omega_b << std::endl;
      std::cout << "omega_cdm = " << omega_cdm << std::endl;
      std::cout << "H0 = " << H0 << std::endl;
      std::cout << "ln10A_s = " << ln10A_s << std::endl;
      std::cout << "n_s = " << n_s << std::endl;
      std::cout << "tau_reio = " << tau_reio << std::endl;

      BEreq::class_input_initialize(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,class_errmsg);

      //DEBUGGING
      //std::cout << pr.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_background_initialize(&pr,&ba);

      //DEBUGGING
      //std::cout << ba.error_message << std::endl;
      //ENDOF DEBUGGING */

      /* for bbn */
      sprintf(pr.sBBN_file,"");
      strcat(pr.sBBN_file,"CosmoBit/src/sBBN.dat");

      BEreq::class_thermodynamics_initialize(&pr,&ba,&th);

      //DEBUGGING
      //std::cout << th.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_perturb_initialize(&pr,&ba,&th,&pt);

      //DEBUGGING
      //std::cout << pt.error_message << std::endl;
      //ENDOF DEBUGGING */

      //DEBUGGING
      //std::cout << tr.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_primordial_initialize(&pr,&pt,&pm);

      //DEBUGGING
      //std::cout << pm.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_nonlinear_initialize(&pr,&ba,&th,&pt,&pm,&nl);

      //DEBUGGING
      //std::cout << nl.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_transfer_initialize(&pr,&ba,&th,&pt,&nl,&tr);

      //DEBUGGING
      //std::cout << sp.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_spectra_initialize(&pr,&ba,&pt,&pm,&nl,&tr,&sp);

      //DEBUGGING
      //std::cout << tr.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_lensing_initialize(&pr,&pt,&sp,&nl,&le);

      //DEBUGGING
      //std::cout << le.error_message << std::endl;
      //ENDOF DEBUGGING */

      /****** write the Cl values in the input array cl[l]  *******/

      for (l=2; l < l_max; l++) {

        int errval = BEreq::class_output_total_cl_at_l(&sp,&le,&op,byVal(l),byVal(cl[l]));

        cl[l][sp.index_ct_tt] = cl[l][sp.index_ct_tt]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_te] = cl[l][sp.index_ct_te]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_ee] = cl[l][sp.index_ct_ee]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_bb] = cl[l][sp.index_ct_bb]*pow(ba.T_cmb*1.e6,2);

      }

      cout << "we are outside the loop for multipoles " << endl;

      BEreq::class_lensing_free(&le);
      BEreq::class_spectra_free(&sp);
      BEreq::class_transfer_free(&tr);
      BEreq::class_nonlinear_free(&nl);
      BEreq::class_primordial_free(&pm);
      BEreq::class_perturb_free(&pt);
      BEreq::class_thermodynamics_free(&th);
      BEreq::class_background_free(&ba);

      return cl;
    }

    double** return_LCDMtensor_cls(double omega_b,double omega_cdm,double H0,double ln10A_s,double n_s,double tau_reio, double r_tensor)
    {
      using namespace Pipes::function_LCDMtensor_lowp_TT_loglike;

      struct Class::precision pr;        /* for precision parameters */
      struct Class::background ba;       /* for cosmological background */
      struct Class::thermo th;           /* for thermodynamics */
      struct Class::perturbs pt;         /* for source functions */
      struct Class::transfers tr;        /* for transfer functions */
      struct Class::primordial pm;       /* for primordial spectra */
      struct Class::spectra sp;          /* for output spectra */
      struct Class::nonlinear nl;        /* for non-linear spectra */
      struct Class::lensing le;          /* for lensed spectra */
      struct Class::output op;           /* for output files */
      Class::ErrorMsg class_errmsg;      /* for error messages */

      int l,l_max;
      int num_ct_max=7;

      l_max=3000;

      char *class_null;// defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
       could add another array with P(k), or extract other results from
       the code - here we assume that we are interested in the C_l's
       only */

      double* cl[l_max];
      for(int i = 0; i < l_max; ++i)
        cl[i] = new double[num_ct_max];

      struct Class::file_content fc;

      BEreq::class_parser_initialize(&fc,11,"",class_errmsg);

      strcpy(fc.name[0],"output");
      strcpy(fc.value[0],"tCl,pCl");

      strcpy(fc.name[1],"omega_b");
      strcpy(fc.name[2],"omega_cdm");
      strcpy(fc.name[3],"H0");
      strcpy(fc.name[4],"ln10^{10}A_s");
      strcpy(fc.name[5],"n_s");
      strcpy(fc.name[6],"tau_reio");

      //strcpy(fc.name[8],"l_scalar_max");
      strcpy(fc.name[8],"l_max_scalars");
      sprintf(fc.value[8],"%d",3000);

      strcpy(fc.name[7],"modes");
      strcpy(fc.value[7],"s,t");

      sprintf(fc.value[1],"%e",omega_b);
      sprintf(fc.value[2],"%e",omega_cdm);
      sprintf(fc.value[3],"%e",H0);
      sprintf(fc.value[4],"%e",ln10A_s);
      sprintf(fc.value[5],"%e",n_s);
      sprintf(fc.value[6],"%e",tau_reio);

      strcpy(fc.name[9],"r");
      sprintf(fc.value[9],"%e",r_tensor);

      std::cout << "omega_b = " << omega_b << std::endl;
      std::cout << "omega_cdm = " << omega_cdm << std::endl;
      std::cout << "H0 = " << H0 << std::endl;
      std::cout << "ln10A_s = " << ln10A_s << std::endl;
      std::cout << "n_s = " << n_s << std::endl;
      std::cout << "tau_reio = " << tau_reio << std::endl;
      std::cout << "r_tensor = " << r_tensor << std::endl;

      BEreq::class_input_initialize(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,class_errmsg);

      //DEBUGGING
      //std::cout << pr.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_background_initialize(&pr,&ba);

      //DEBUGGING
      //std::cout << ba.error_message << std::endl;
      //ENDOF DEBUGGING */

      /* for bbn */
      sprintf(pr.sBBN_file,"");
      strcat(pr.sBBN_file,"CosmoBit/src/sBBN.dat");

      BEreq::class_thermodynamics_initialize(&pr,&ba,&th);

      //DEBUGGING
      //std::cout << th.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_perturb_initialize(&pr,&ba,&th,&pt);

      //DEBUGGING
      //std::cout << pt.error_message << std::endl;
      //ENDOF DEBUGGING */

      //DEBUGGING
      //std::cout << tr.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_primordial_initialize(&pr,&pt,&pm);

      //DEBUGGING
      //std::cout << pm.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_nonlinear_initialize(&pr,&ba,&th,&pt,&pm,&nl);

      //DEBUGGING
      //std::cout << nl.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_transfer_initialize(&pr,&ba,&th,&pt,&nl,&tr);

      //DEBUGGING
      //std::cout << sp.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_spectra_initialize(&pr,&ba,&pt,&pm,&nl,&tr,&sp);

      //DEBUGGING
      //std::cout << tr.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_lensing_initialize(&pr,&pt,&sp,&nl,&le);

      //DEBUGGING
      //std::cout << le.error_message << std::endl;
      //ENDOF DEBUGGING */

      /****** write the Cl values in the input array cl[l]  *******/

      for (l=2; l < l_max; l++)
      {
        int errval = BEreq::class_output_total_cl_at_l(&sp,&le,&op,byVal(l),byVal(cl[l]));

        cl[l][sp.index_ct_tt] = cl[l][sp.index_ct_tt]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_te] = cl[l][sp.index_ct_te]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_ee] = cl[l][sp.index_ct_ee]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_bb] = cl[l][sp.index_ct_bb]*pow(ba.T_cmb*1.e6,2);
      }

      cout << "we are outside the loop for multipoles " << endl;

      BEreq::class_lensing_free(&le);
      BEreq::class_spectra_free(&sp);
      BEreq::class_transfer_free(&tr);
      BEreq::class_nonlinear_free(&nl);
      BEreq::class_primordial_free(&pm);
      BEreq::class_perturb_free(&pt);
      BEreq::class_thermodynamics_free(&th);
      BEreq::class_background_free(&ba);

      return cl;
    }

    double** return_LCDMtensor_inflation_cls(double omega_b,double omega_cdm,double H0,double ln10A_s,double n_s,double tau_reio, double r_tensor)
    {
      using namespace Pipes::function_LCDMtensor_inflation_lowp_TT_loglike;

      struct Class::precision pr;        /* for precision parameters */
      struct Class::background ba;       /* for cosmological background */
      struct Class::thermo th;           /* for thermodynamics */
      struct Class::perturbs pt;         /* for source functions */
      struct Class::transfers tr;        /* for transfer functions */
      struct Class::primordial pm;       /* for primordial spectra */
      struct Class::spectra sp;          /* for output spectra */
      struct Class::nonlinear nl;        /* for non-linear spectra */
      struct Class::lensing le;          /* for lensed spectra */
      struct Class::output op;           /* for output files */
      Class::ErrorMsg class_errmsg;      /* for error messages */

      int l,l_max;
      int num_ct_max=7;

      l_max=3000;

      char *class_null;// defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
       could add another array with P(k), or extract other results from
       the code - here we assume that we are interested in the C_l's
       only */

      double* cl[l_max];
      for(int i = 0; i < l_max; ++i)
        cl[i] = new double[num_ct_max];

      struct Class::file_content fc;

      BEreq::class_parser_initialize(&fc,11,"",class_errmsg);


      strcpy(fc.name[0],"output");
      strcpy(fc.value[0],"tCl,pCl");

      strcpy(fc.name[1],"omega_b");
      strcpy(fc.name[2],"omega_cdm");
      strcpy(fc.name[3],"H0");
      strcpy(fc.name[4],"ln10^{10}A_s");
      strcpy(fc.name[5],"n_s");
      strcpy(fc.name[6],"tau_reio");

      //strcpy(fc.name[8],"l_scalar_max");
      strcpy(fc.name[8],"l_max_scalars");
      sprintf(fc.value[8],"%d",3000);

      strcpy(fc.name[7],"modes");
      sprintf(fc.value[7],"s,t");

      sprintf(fc.value[1],"%e",omega_b);
      sprintf(fc.value[2],"%e",omega_cdm);
      sprintf(fc.value[3],"%e",H0);
      sprintf(fc.value[4],"%e",ln10A_s);
      sprintf(fc.value[5],"%e",n_s);
      sprintf(fc.value[6],"%e",tau_reio);

      sprintf(fc.value[9],"%e",r_tensor);

      std::cout << "omega_b = " << omega_b << std::endl;
      std::cout << "omega_cdm = " << omega_cdm << std::endl;
      std::cout << "H0 = " << H0 << std::endl;
      std::cout << "ln10A_s = " << ln10A_s << std::endl;
      std::cout << "n_s = " << n_s << std::endl;
      std::cout << "tau_reio = " << tau_reio << std::endl;
      std::cout << "r_tensor = " << r_tensor << std::endl;

      BEreq::class_input_initialize(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,class_errmsg);

      //DEBUGGING
      //std::cout << pr.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_background_initialize(&pr,&ba);

      //DEBUGGING
      //std::cout << ba.error_message << std::endl;
      //ENDOF DEBUGGING */

      /* for bbn */
      sprintf(pr.sBBN_file,"");
      strcat(pr.sBBN_file,"CosmoBit/src/sBBN.dat");

      BEreq::class_thermodynamics_initialize(&pr,&ba,&th);

      //DEBUGGING
      //std::cout << th.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_perturb_initialize(&pr,&ba,&th,&pt);

      //DEBUGGING
      //std::cout << pt.error_message << std::endl;
      //ENDOF DEBUGGING */

      //DEBUGGING
      //std::cout << tr.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_primordial_initialize(&pr,&pt,&pm);

      //DEBUGGING
      //std::cout << pm.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_nonlinear_initialize(&pr,&ba,&th,&pt,&pm,&nl);

      //DEBUGGING
      //std::cout << nl.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_transfer_initialize(&pr,&ba,&th,&pt,&nl,&tr);

      //DEBUGGING
      //std::cout << sp.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_spectra_initialize(&pr,&ba,&pt,&pm,&nl,&tr,&sp);

      //DEBUGGING
      //std::cout << tr.error_message << std::endl;
      //ENDOF DEBUGGING */

      BEreq::class_lensing_initialize(&pr,&pt,&sp,&nl,&le);

      //DEBUGGING
      //std::cout << le.error_message << std::endl;
      //ENDOF DEBUGGING */

      /****** write the Cl values in the input array cl[l]  *******/

      for (l=2; l < l_max; l++)
      {
        int errval = BEreq::class_output_total_cl_at_l(&sp,&le,&op,byVal(l),byVal(cl[l]));

        cl[l][sp.index_ct_tt] = cl[l][sp.index_ct_tt]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_te] = cl[l][sp.index_ct_te]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_ee] = cl[l][sp.index_ct_ee]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_bb] = cl[l][sp.index_ct_bb]*pow(ba.T_cmb*1.e6,2);
      }

      cout << "we are outside the loop for multipoles " << endl;

      BEreq::class_lensing_free(&le);
      BEreq::class_spectra_free(&sp);
      BEreq::class_transfer_free(&tr);
      BEreq::class_nonlinear_free(&nl);
      BEreq::class_primordial_free(&pm);
      BEreq::class_perturb_free(&pt);
      BEreq::class_thermodynamics_free(&th);
      BEreq::class_background_free(&ba);

      return cl;
    }

    double** return_LCDMtensor_1quadInf_cls(double omega_b,double omega_cdm,double H0,double A_s,double n_s,double tau_reio, double r_tensor)
    {
      using namespace Pipes::function_1quadInf_rLCDMtensor_lowp_TT_loglike;

      struct Class::precision pr;        /* for precision parameters */
      struct Class::background ba;       /* for cosmological background */
      struct Class::thermo th;           /* for thermodynamics */
      struct Class::perturbs pt;         /* for source functions */
      struct Class::transfers tr;        /* for transfer functions */
      struct Class::primordial pm;       /* for primordial spectra */
      struct Class::spectra sp;          /* for output spectra */
      struct Class::nonlinear nl;        /* for non-linear spectra */
      struct Class::lensing le;          /* for lensed spectra */
      struct Class::output op;           /* for output files */
      Class::ErrorMsg class_errmsg;      /* for error messages */

      int l,l_max;
      int num_ct_max=7;

      l_max=3000;

      char *class_null;// defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
      could add another array with P(k), or extract other results from
      the code - here we assume that we are interested in the C_l's
      only */

      double* cl[l_max];
      for(int i = 0; i < l_max; ++i)
        cl[i] = new double[num_ct_max];

      struct Class::file_content fc;

      BEreq::class_parser_initialize(&fc,11,"",class_errmsg);

      strcpy(fc.name[0],"output");
      strcpy(fc.name[1],"omega_b");
      strcpy(fc.name[2],"omega_cdm");
      strcpy(fc.name[3],"H0");
      strcpy(fc.name[4],"A_s");
      strcpy(fc.name[5],"n_s");
      strcpy(fc.name[6],"tau_reio");
      strcpy(fc.name[7],"modes");
      strcpy(fc.name[8],"l_scalar_max");
      strcpy(fc.name[9],"r");

      strcpy(fc.value[0],"tCl,pCl");
      sprintf(fc.value[1],"%e",omega_b);
      sprintf(fc.value[2],"%e",omega_cdm);
      sprintf(fc.value[3],"%e",H0);
      sprintf(fc.value[4],"%e",A_s);
      sprintf(fc.value[5],"%e",n_s);
      sprintf(fc.value[6],"%e",tau_reio);
      sprintf(fc.value[7],"s,t");
      sprintf(fc.value[8],"%d",3000);
      sprintf(fc.value[9],"%e",r_tensor);

      std::cout << "Here are the theory parameter values I generate the Cl's with" << std::endl;
      std::cout << "omega_b = " << omega_b << std::endl;
      std::cout << "omega_cdm = " << omega_cdm << std::endl;
      std::cout << "H0 = " << H0 << std::endl;
      std::cout << "A_s = " << A_s << std::endl;
      std::cout << "n_s = " << n_s << std::endl;
      std::cout << "tau_reio = " << tau_reio << std::endl;
      std::cout << "r_tensor = " << r_tensor << std::endl;

      BEreq::class_input_initialize(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,class_errmsg);
      BEreq::class_background_initialize(&pr,&ba);

      /* for bbn */
      sprintf(pr.sBBN_file,"");
      strcat(pr.sBBN_file,"/Users/selimhotinli/Dropbox/gambit/Backends/installed/class/2.6.1/bbn/sBBN.dat");

      BEreq::class_thermodynamics_initialize(&pr,&ba,&th);
      BEreq::class_perturb_initialize(&pr,&ba,&th,&pt);
      BEreq::class_primordial_initialize(&pr,&pt,&pm);
      BEreq::class_nonlinear_initialize(&pr,&ba,&th,&pt,&pm,&nl);
      BEreq::class_transfer_initialize(&pr,&ba,&th,&pt,&nl,&tr);
      BEreq::class_spectra_initialize(&pr,&ba,&pt,&pm,&nl,&tr,&sp);
      BEreq::class_lensing_initialize(&pr,&pt,&sp,&nl,&le);

      /****** write the Cl values in the input array cl[l]  *******/

      for (l=2; l < l_max; l++)
      {
        int errval = BEreq::class_output_total_cl_at_l(&sp,&le,&op,byVal(l),byVal(cl[l]));

        cl[l][sp.index_ct_tt] = cl[l][sp.index_ct_tt]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_te] = cl[l][sp.index_ct_te]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_ee] = cl[l][sp.index_ct_ee]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_bb] = cl[l][sp.index_ct_bb]*pow(ba.T_cmb*1.e6,2);

        cout << " cl["<<l<<"]["<<sp.index_ct_tt<<"] ="<< cl[l][sp.index_ct_tt]<< endl;
      }

      cout << "we are outside the loop for multipoles " << endl;

      BEreq::class_lensing_free(&le);
      BEreq::class_spectra_free(&sp);
      BEreq::class_transfer_free(&tr);
      BEreq::class_nonlinear_free(&nl);
      BEreq::class_primordial_free(&pm);
      BEreq::class_perturb_free(&pt);
      BEreq::class_thermodynamics_free(&th);
      BEreq::class_background_free(&ba);

      return cl;
    }

    /* ----------------   SH ---------------- */
    // testing feeding CLASS the full power spectrum
    double** return_LCDMtensor_1quadInf_cls_fullPk(double omega_b,
                                                   double omega_cdm,
                                                   double H0,
                                                   double tau_reio,
                                                   double k_array[],
                                                   double pks_array[],
                                                   double pkt_array[],
                                                   int k_array_size)
    {
      using namespace Pipes::function_1quadInf_rLCDMtensor_lowp_TT_loglike;

      struct Class::precision pr;        /* for precision parameters */
      struct Class::background ba;       /* for cosmological background */
      struct Class::thermo th;           /* for thermodynamics */
      struct Class::perturbs pt;         /* for source functions */
      struct Class::transfers tr;        /* for transfer functions */
      struct Class::primordial pm;       /* for primordial spectra */
      struct Class::spectra sp;          /* for output spectra */
      struct Class::nonlinear nl;        /* for non-linear spectra */
      struct Class::lensing le;          /* for lensed spectra */
      struct Class::output op;           /* for output files */
      Class::ErrorMsg class_errmsg;      /* for error messages */

      int l,l_max;
      int num_ct_max=7;

      l_max=3000;

      char *class_null;// defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
      could add another array with P(k), or extract other results from
      the code - here we assume that we are interested in the C_l's
      only */
      std::cout << "We are inside!" << std::endl;

//      std::cout << "k_array = " << (*k_array)[0] << std::endl;
//      std::cout << "pks_array = " << (*pks_array)[0] << std::endl;
//      std::cout << "pkt_array = " << (*pkt_array)[0] << std::endl;

      double* cl[l_max];
      for(int i = 0; i < l_max; ++i)
        cl[i] = new double[num_ct_max];

      struct Class::file_content fc;

      BEreq::class_parser_initialize(&fc,11,"",class_errmsg);

      strcpy(fc.name[0],"output");
      strcpy(fc.name[1],"omega_b");
      strcpy(fc.name[2],"omega_cdm");
      strcpy(fc.name[3],"H0");

      //strcpy(fc.name[4],"A_s");
      //strcpy(fc.name[5],"n_s");
      strcpy(fc.name[4], "P_k_ini type");
      strcpy(fc.name[5], "primordial_verbose");

      strcpy(fc.name[6],"tau_reio");
      strcpy(fc.name[7],"modes");
      strcpy(fc.name[8],"l_scalar_max");
      //strcpy(fc.name[9],"r");

      strcpy(fc.value[0],"tCl,pCl");
      sprintf(fc.value[1],"%e",omega_b);
      sprintf(fc.value[2],"%e",omega_cdm);
      sprintf(fc.value[3],"%e",H0);

      //sprintf(fc.value[4],"%e",A_s);
      //sprintf(fc.value[5],"%e",n_s);
      strcpy(fc.value[4],"gambit_Pk"); // patched CLASS to have gambit_Pk option.
      sprintf(fc.value[5],"%d",1);

      sprintf(fc.value[6],"%e",tau_reio);
      sprintf(fc.value[7],"s,t");
      sprintf(fc.value[8],"%d",3000);
      //sprintf(fc.value[9],"%e",r_tensor);

      std::cout << "Here are the theory parameter values I generate the Cl's with" << std::endl;
      std::cout << "omega_b = " << omega_b << std::endl;
      std::cout << "omega_cdm = " << omega_cdm << std::endl;
      std::cout << "H0 = " << H0 << std::endl;
      //std::cout << "A_s = " << A_s << std::endl;
      //std::cout << "n_s = " << n_s << std::endl;
      std::cout << "tau_reio = " << tau_reio << std::endl;
      //std::cout << "r_tensor = " << r_tensor << std::endl;

      BEreq::class_input_initialize(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,class_errmsg);
      BEreq::class_background_initialize(&pr,&ba);

      /* for bbn */
      sprintf(pr.sBBN_file,"");
      strcat(pr.sBBN_file,"/Users/selimhotinli/Dropbox/gambit/Backends/installed/class/2.6.1/bbn/sBBN.dat");

      BEreq::class_thermodynamics_initialize(&pr,&ba,&th);
      BEreq::class_perturb_initialize(&pr,&ba,&th,&pt);

      std::cout << "we are fine after initialized and upto perturb-init " << std::endl;

      //BEreq::class_primordial_initialize(&pr,&pt,&pm);
      //similar to the lines 3392-3412 in primordial.c CLASS

      /** - Make room */
      pm.lnk_size = k_array_size;

      pm.lnk = (double *)malloc(k_array_size*sizeof(double));

      std::cout << "we pass pm.lnk\n" << std::endl;

      pm.lnpk = (double **)malloc(pt.md_size*sizeof(double));
      pm.ddlnpk = (double **)malloc(pt.md_size*sizeof(double));
      pm.ic_size = (int *)malloc(pt.md_size*sizeof(int));
      pm.ic_ic_size = (int *)malloc(pt.md_size*sizeof(int));
      pm.is_non_zero = (short **)malloc(pt.md_size*sizeof(short));

      int index_md;
      for (index_md = 0; index_md < pt.md_size; index_md++)
      {
        pm.ic_size[index_md] = pt.ic_size[index_md];
        pm.ic_ic_size[index_md] = (pm.ic_size[index_md]*(pm.ic_size[index_md]+1))/2;

        std::cout << "pm.ic_size["<<index_md<<"] =  " << pt.ic_size[index_md] << std::endl;
        std::cout << "pm.ic_ic_size["<<index_md<<"] = " << pm.ic_ic_size[index_md]<<std::endl;

        pm.lnpk[index_md] = (double *)malloc(pm.lnk_size*pm.ic_ic_size[index_md]*sizeof(double));
        pm.ddlnpk[index_md] = (double *)malloc(pm.lnk_size*pm.ic_ic_size[index_md]*sizeof(double));
        pm.is_non_zero[index_md] = (short *)malloc(pm.lnk_size*pm.ic_ic_size[index_md]*sizeof(short));
      }

//      pm.lnpk[pt.index_md_scalars] = (double *)malloc(k_array_size*sizeof(double));
//      std::cout << "we pass pm.lnpk " << std::endl;
//      pm.ddlnpk[pt.index_md_scalars] = (double *)malloc(k_array_size*sizeof(double));
//      // SH: We shall assume we have tensors --> can correct this as it is post input-initialize.
//      pm.lnpk[pt.index_md_tensors] = (double *)malloc(k_array_size*sizeof(double));
//      pm.ddlnpk[pt.index_md_tensors] = (double *)malloc(k_array_size*sizeof(double));

      std::cout << "we pass all " << std::endl;

      std::cout << "we are fine pass -make room- " << std::endl;

      /** - Store values */
      for (int index_k=0; index_k<pm.lnk_size; index_k++)
      {
        std::cout << "k_array["<<index_k<<"]="<<k_array[index_k]<<std::endl;
        std::cout << "pks_array["<<index_k<<"]="<<pks_array[index_k]<<std::endl;
        std::cout << "pkt_array["<<index_k<<"]="<<pkt_array[index_k]<<std::endl;

        pm.lnk[index_k] = std::log(k_array[index_k]);
        pm.lnpk[pt.index_md_scalars][index_k] = std::log(pks_array[index_k]);
        if (pt.has_tensors == _TRUE_)
          pm.lnpk[pt.index_md_tensors][index_k] = std::log(pkt_array[index_k]);

        std::cout << "pm.lnk["<<index_k<<"]="<<pm.lnk[index_k]<<std::endl;
        std::cout << "pm.lnpk["<<pt.index_md_scalars<<"]["<<index_k<<"]="<<pm.lnpk[pt.index_md_scalars][index_k]<<std::endl;
        std::cout << "pm.lnpk["<<pt.index_md_tensors<<"]["<<index_k<<"]="<<pm.lnpk[pt.index_md_tensors][index_k]<<std::endl;
      }

      std::cout << "we are fine pass -store values- " << std::endl;

      /** - Tell CLASS that there are scalar (and tensor) modes */
      pm.is_non_zero[pt.index_md_scalars][pt.index_ic_ad] = _TRUE_;
      if (pt.has_tensors == _TRUE_)
        pm.is_non_zero[pt.index_md_tensors][pt.index_ic_ten] = _TRUE_;

      BEreq::class_primordial_initialize(&pr,&pt,&pm);

      std::cout << "we are out of class_primordial_initialize " << std::endl;

      //DEBUGGING
//      std::cout << pm.error_message << std::endl;
      //ENDOF DEBUGGING

      BEreq::class_nonlinear_initialize(&pr,&ba,&th,&pt,&pm,&nl);

      std::cout << "we are out of class_nonlinear_initialize " << std::endl;

      //DEBUGGING
//      std::cout << nl.error_message << std::endl;
      //ENDOF DEBUGGING

      BEreq::class_transfer_initialize(&pr,&ba,&th,&pt,&nl,&tr);

      std::cout << "tr.l_size[0] = "<<tr.l_size[0]<<std::endl;
      std::cout << "tr.l_size[1] = "<<tr.l_size[1]<<std::endl;

      std::cout << "we are out of class_transfer_initialize " << std::endl;

      //DEBUGGING
//      std::cout << tr.error_message << std::endl;
      //ENDOF DEBUGGING

      BEreq::class_spectra_initialize(&pr,&ba,&pt,&pm,&nl,&tr,&sp);

      std::cout << "sp.l_size[0] = "<<sp.l_size[0]<<std::endl;
      std::cout << "sp.l_size[1] = "<<sp.l_size[1]<<std::endl;

      std::cout << "we are out of class_spectra_initialize " << std::endl;

      //DEBUGGING
//      std::cout << sp.error_message << std::endl;
      //ENDOF DEBUGGING

      BEreq::class_lensing_initialize(&pr,&pt,&sp,&nl,&le);

      std::cout << "we are out of class_lensing_initialize " << std::endl;

      //DEBUGGING
//      std::cout << le.error_message << std::endl;
      //ENDOF DEBUGGING

      /****** write the Cl values in the input array cl[l]  *******/

      for (l=2; l < l_max; l++)
      {
        int errval = BEreq::class_output_total_cl_at_l(&sp,&le,&op,byVal(l),byVal(cl[l]));
//        cout << "we are outside class_output_total_cl_at_l " << endl;

        cl[l][sp.index_ct_tt] = cl[l][sp.index_ct_tt]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_te] = cl[l][sp.index_ct_te]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_ee] = cl[l][sp.index_ct_ee]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_bb] = cl[l][sp.index_ct_bb]*pow(ba.T_cmb*1.e6,2);

//        std::cout << "cl["<<l<<"]["<<sp.index_ct_tt<<"] ="<< cl[l][sp.index_ct_tt] << std::endl;
//        std::cout << "cl["<<l<<"]["<<sp.index_ct_te<<"] ="<< cl[l][sp.index_ct_te] << std::endl;
//        std::cout << "cl["<<l<<"]["<<sp.index_ct_ee<<"] ="<< cl[l][sp.index_ct_ee] << std::endl;
//        std::cout << "cl["<<l<<"]["<<sp.index_ct_bb<<"] ="<< cl[l][sp.index_ct_bb] << std::endl;
      }

      cout << "we are outside the loop for multipoles " << endl;

      BEreq::class_lensing_free(&le);

      std::cout << "we are out of class_lensing_free " << std::endl;

      BEreq::class_spectra_free(&sp);

      std::cout << "we are out of class_spectra_free " << std::endl;

      BEreq::class_transfer_free(&tr);

      std::cout << "we are out of class_transfer_free " << std::endl;

      BEreq::class_nonlinear_free(&nl);

      std::cout << "we are out of class_nonlinear_free " << std::endl;

      BEreq::class_primordial_free(&pm);

      std::cout << "we are out of class_primordial_free " << std::endl;

      BEreq::class_perturb_free(&pt);

      std::cout << "we are out of class_perturb_free " << std::endl;

      BEreq::class_thermodynamics_free(&th);

      std::cout << "we are out of class_thermodynamics_free " << std::endl;

      BEreq::class_background_free(&ba);

      std::cout << "we are out of class_background_free " << std::endl;

      return cl;
    }

    double** return_LCDMtensor_1quarInf_cls(double omega_b,double omega_cdm,double H0,double A_s,double n_s,double tau_reio, double r_tensor)
    {
      using namespace Pipes::function_1quarInf_rLCDMtensor_lowp_TT_loglike;
      using namespace Class;

      struct Class::precision pr;        /* for precision parameters */
      struct Class::background ba;       /* for cosmological background */
      struct Class::thermo th;           /* for thermodynamics */
      struct Class::perturbs pt;         /* for source functions */
      struct Class::transfers tr;        /* for transfer functions */
      struct Class::primordial pm;       /* for primordial spectra */
      struct Class::spectra sp;          /* for output spectra */
      struct Class::nonlinear nl;        /* for non-linear spectra */
      struct Class::lensing le;          /* for lensed spectra */
      struct Class::output op;           /* for output files */
      Class::ErrorMsg class_errmsg;      /* for error messages */

      int l,l_max;
      int num_ct_max=7;

      l_max=3000;

      char *class_null;// defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
      could add another array with P(k), or extract other results from
      the code - here we assume that we are interested in the C_l's
      only */

      double* cl[l_max];
      for(int i = 0; i < l_max; ++i)
        cl[i] = new double[num_ct_max];

      struct Class::file_content fc;

      BEreq::class_parser_initialize(&fc,11,"",class_errmsg);

      strcpy(fc.name[0],"output");
      strcpy(fc.name[1],"omega_b");
      strcpy(fc.name[2],"omega_cdm");
      strcpy(fc.name[3],"H0");
      strcpy(fc.name[4],"A_s");
      strcpy(fc.name[5],"n_s");
      strcpy(fc.name[6],"tau_reio");
      strcpy(fc.name[7],"modes");
      strcpy(fc.name[8],"l_scalar_max");
      strcpy(fc.name[9],"r");

      strcpy(fc.value[0],"tCl,pCl");
      sprintf(fc.value[1],"%e",omega_b);
      sprintf(fc.value[2],"%e",omega_cdm);
      sprintf(fc.value[3],"%e",H0);
      sprintf(fc.value[4],"%e",A_s);
      sprintf(fc.value[5],"%e",n_s);
      sprintf(fc.value[6],"%e",tau_reio);
      sprintf(fc.value[7],"s,t");
      sprintf(fc.value[8],"%d",3000);
      sprintf(fc.value[9],"%e",r_tensor);

      std::cout << "Here are the theory parameter values I generate the Cl's with" << std::endl;
      std::cout << "omega_b = " << omega_b << std::endl;
      std::cout << "omega_cdm = " << omega_cdm << std::endl;
      std::cout << "H0 = " << H0 << std::endl;
      std::cout << "A_s = " << A_s << std::endl;
      std::cout << "n_s = " << n_s << std::endl;
      std::cout << "tau_reio = " << tau_reio << std::endl;
      std::cout << "r_tensor = " << r_tensor << std::endl;

      BEreq::class_input_initialize(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,class_errmsg);
      BEreq::class_background_initialize(&pr,&ba);

      /* for bbn */
      sprintf(pr.sBBN_file,"");
      strcat(pr.sBBN_file,"/Users/selimhotinli/Dropbox/gambit/Backends/installed/class/2.6.1/bbn/sBBN.dat");

      BEreq::class_thermodynamics_initialize(&pr,&ba,&th);
      BEreq::class_perturb_initialize(&pr,&ba,&th,&pt);
      BEreq::class_primordial_initialize(&pr,&pt,&pm);
      BEreq::class_nonlinear_initialize(&pr,&ba,&th,&pt,&pm,&nl);
      BEreq::class_transfer_initialize(&pr,&ba,&th,&pt,&nl,&tr);
      BEreq::class_spectra_initialize(&pr,&ba,&pt,&pm,&nl,&tr,&sp);
      BEreq::class_lensing_initialize(&pr,&pt,&sp,&nl,&le);

      /****** write the Cl values in the input array cl[l]  *******/

      for (l=2; l < l_max; l++)
      {
        int errval = BEreq::class_output_total_cl_at_l(&sp,&le,&op,byVal(l),byVal(cl[l]));

        cl[l][sp.index_ct_tt] = cl[l][sp.index_ct_tt]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_te] = cl[l][sp.index_ct_te]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_ee] = cl[l][sp.index_ct_ee]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_bb] = cl[l][sp.index_ct_bb]*pow(ba.T_cmb*1.e6,2);
      }

      cout << "we are outside the loop for multipoles " << endl;

      BEreq::class_lensing_free(&le);
      BEreq::class_spectra_free(&sp);
      BEreq::class_transfer_free(&tr);
      BEreq::class_nonlinear_free(&nl);
      BEreq::class_primordial_free(&pm);
      BEreq::class_perturb_free(&pt);
      BEreq::class_thermodynamics_free(&th);
      BEreq::class_background_free(&ba);

      return cl;
    }

    double** return_LCDMtensor_1mono32Inf_cls(double omega_b,double omega_cdm,double H0,double A_s,double n_s,double tau_reio, double r_tensor)
    {
      using namespace Pipes::function_1mono32Inf_rLCDMtensor_lowp_TT_loglike;

      struct Class::precision pr;        /* for precision parameters */
      struct Class::background ba;       /* for cosmological background */
      struct Class::thermo th;           /* for thermodynamics */
      struct Class::perturbs pt;         /* for source functions */
      struct Class::transfers tr;        /* for transfer functions */
      struct Class::primordial pm;       /* for primordial spectra */
      struct Class::spectra sp;          /* for output spectra */
      struct Class::nonlinear nl;        /* for non-linear spectra */
      struct Class::lensing le;          /* for lensed spectra */
      struct Class::output op;           /* for output files */
      Class::ErrorMsg class_errmsg;      /* for error messages */

      int l,l_max;
      int num_ct_max=7;

      l_max=3000;

      char *class_null;// defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
      could add another array with P(k), or extract other results from
      the code - here we assume that we are interested in the C_l's
      only */

      double* cl[l_max];
      for(int i = 0; i < l_max; ++i)
      cl[i] = new double[num_ct_max];

      struct Class::file_content fc;

      BEreq::class_parser_initialize(&fc,11,"",class_errmsg);

      strcpy(fc.name[0],"output");
      strcpy(fc.name[1],"omega_b");
      strcpy(fc.name[2],"omega_cdm");
      strcpy(fc.name[3],"H0");
      strcpy(fc.name[4],"A_s");
      strcpy(fc.name[5],"n_s");
      strcpy(fc.name[6],"tau_reio");
      strcpy(fc.name[7],"modes");
      strcpy(fc.name[8],"l_scalar_max");
      strcpy(fc.name[9],"r");

      strcpy(fc.value[0],"tCl,pCl");
      sprintf(fc.value[1],"%e",omega_b);
      sprintf(fc.value[2],"%e",omega_cdm);
      sprintf(fc.value[3],"%e",H0);
      sprintf(fc.value[4],"%e",A_s);
      sprintf(fc.value[5],"%e",n_s);
      sprintf(fc.value[6],"%e",tau_reio);
      sprintf(fc.value[7],"s,t");
      sprintf(fc.value[8],"%d",3000);
      sprintf(fc.value[9],"%e",r_tensor);

      std::cout << "Here are the theory parameter values I generate the Cl's with" << std::endl;
      std::cout << "omega_b = " << omega_b << std::endl;
      std::cout << "omega_cdm = " << omega_cdm << std::endl;
      std::cout << "H0 = " << H0 << std::endl;
      std::cout << "A_s = " << A_s << std::endl;
      std::cout << "n_s = " << n_s << std::endl;
      std::cout << "tau_reio = " << tau_reio << std::endl;
      std::cout << "r_tensor = " << r_tensor << std::endl;

      BEreq::class_input_initialize(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,class_errmsg);
      BEreq::class_background_initialize(&pr,&ba);

      /* for bbn */
      sprintf(pr.sBBN_file,"");
      strcat(pr.sBBN_file,"/Users/selimhotinli/Dropbox/gambit/Backends/installed/class/2.6.1/bbn/sBBN.dat");

      BEreq::class_thermodynamics_initialize(&pr,&ba,&th);
      BEreq::class_perturb_initialize(&pr,&ba,&th,&pt);
      BEreq::class_primordial_initialize(&pr,&pt,&pm);
      BEreq::class_nonlinear_initialize(&pr,&ba,&th,&pt,&pm,&nl);
      BEreq::class_transfer_initialize(&pr,&ba,&th,&pt,&nl,&tr);
      BEreq::class_spectra_initialize(&pr,&ba,&pt,&pm,&nl,&tr,&sp);
      BEreq::class_lensing_initialize(&pr,&pt,&sp,&nl,&le);

      /****** write the Cl values in the input array cl[l]  *******/

      for (l=2; l < l_max; l++)
      {
        int errval = BEreq::class_output_total_cl_at_l(&sp,&le,&op,byVal(l),byVal(cl[l]));

        cl[l][sp.index_ct_tt] = cl[l][sp.index_ct_tt]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_te] = cl[l][sp.index_ct_te]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_ee] = cl[l][sp.index_ct_ee]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_bb] = cl[l][sp.index_ct_bb]*pow(ba.T_cmb*1.e6,2);
      }

      cout << "we are outside the loop for multipoles " << endl;

      BEreq::class_lensing_free(&le);
      BEreq::class_spectra_free(&sp);
      BEreq::class_transfer_free(&tr);
      BEreq::class_nonlinear_free(&nl);
      BEreq::class_primordial_free(&pm);
      BEreq::class_perturb_free(&pt);
      BEreq::class_thermodynamics_free(&th);
      BEreq::class_background_free(&ba);

      return cl;
    }

    double** return_LCDMtensor_1linearInf_cls(double omega_b,double omega_cdm,double H0,double A_s,double n_s,double tau_reio, double r_tensor)
    {
      using namespace Pipes::function_1linearInf_rLCDMtensor_lowp_TT_loglike;

      struct Class::precision pr;        /* for precision parameters */
      struct Class::background ba;       /* for cosmological background */
      struct Class::thermo th;           /* for thermodynamics */
      struct Class::perturbs pt;         /* for source functions */
      struct Class::transfers tr;        /* for transfer functions */
      struct Class::primordial pm;       /* for primordial spectra */
      struct Class::spectra sp;          /* for output spectra */
      struct Class::nonlinear nl;        /* for non-linear spectra */
      struct Class::lensing le;          /* for lensed spectra */
      struct Class::output op;           /* for output files */
      Class::ErrorMsg class_errmsg;      /* for error messages */

      int l,l_max;
      int num_ct_max=7;

      l_max=3000;

      char *class_null;// defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
      could add another array with P(k), or extract other results from
      the code - here we assume that we are interested in the C_l's
      only */

      double* cl[l_max];
      for(int i = 0; i < l_max; ++i)
        cl[i] = new double[num_ct_max];

      struct Class::file_content fc;

      BEreq::class_parser_initialize(&fc,11,"",class_errmsg);

      strcpy(fc.name[0],"output");
      strcpy(fc.name[1],"omega_b");
      strcpy(fc.name[2],"omega_cdm");
      strcpy(fc.name[3],"H0");
      strcpy(fc.name[4],"A_s");
      strcpy(fc.name[5],"n_s");
      strcpy(fc.name[6],"tau_reio");
      strcpy(fc.name[7],"modes");
      strcpy(fc.name[8],"l_scalar_max");
      strcpy(fc.name[9],"r");

      strcpy(fc.value[0],"tCl,pCl");
      sprintf(fc.value[1],"%e",omega_b);
      sprintf(fc.value[2],"%e",omega_cdm);
      sprintf(fc.value[3],"%e",H0);
      sprintf(fc.value[4],"%e",A_s);
      sprintf(fc.value[5],"%e",n_s);
      sprintf(fc.value[6],"%e",tau_reio);
      sprintf(fc.value[7],"s,t");
      sprintf(fc.value[8],"%d",3000);
      sprintf(fc.value[9],"%e",r_tensor);

      std::cout << "Here are the theory parameter values I generate the Cl's with" << std::endl;
      std::cout << "omega_b = " << omega_b << std::endl;
      std::cout << "omega_cdm = " << omega_cdm << std::endl;
      std::cout << "H0 = " << H0 << std::endl;
      std::cout << "A_s = " << A_s << std::endl;
      std::cout << "n_s = " << n_s << std::endl;
      std::cout << "tau_reio = " << tau_reio << std::endl;
      std::cout << "r_tensor = " << r_tensor << std::endl;

      BEreq::class_input_initialize(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,class_errmsg);
      BEreq::class_background_initialize(&pr,&ba);

      /* for bbn */
      sprintf(pr.sBBN_file,"");
      strcat(pr.sBBN_file,"/Users/selimhotinli/Dropbox/gambit/Backends/installed/class/2.6.1/bbn/sBBN.dat");

      BEreq::class_thermodynamics_initialize(&pr,&ba,&th);
      BEreq::class_perturb_initialize(&pr,&ba,&th,&pt);
      BEreq::class_primordial_initialize(&pr,&pt,&pm);
      BEreq::class_nonlinear_initialize(&pr,&ba,&th,&pt,&pm,&nl);
      BEreq::class_transfer_initialize(&pr,&ba,&th,&pt,&nl,&tr);
      BEreq::class_spectra_initialize(&pr,&ba,&pt,&pm,&nl,&tr,&sp);
      BEreq::class_lensing_initialize(&pr,&pt,&sp,&nl,&le);

      /****** write the Cl values in the input array cl[l]  *******/

      for (l=2; l < l_max; l++)
      {
        int errval = BEreq::class_output_total_cl_at_l(&sp,&le,&op,byVal(l),byVal(cl[l]));

        cl[l][sp.index_ct_tt] = cl[l][sp.index_ct_tt]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_te] = cl[l][sp.index_ct_te]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_ee] = cl[l][sp.index_ct_ee]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_bb] = cl[l][sp.index_ct_bb]*pow(ba.T_cmb*1.e6,2);

      }

      cout << "we are outside the loop for multipoles " << endl;

      BEreq::class_lensing_free(&le);
      BEreq::class_spectra_free(&sp);
      BEreq::class_transfer_free(&tr);
      BEreq::class_nonlinear_free(&nl);
      BEreq::class_primordial_free(&pm);
      BEreq::class_perturb_free(&pt);
      BEreq::class_thermodynamics_free(&th);
      BEreq::class_background_free(&ba);

      return cl;
    }

    double** return_LCDMtensor_smashInf_cls(double omega_b,double omega_cdm,double H0,double A_s,double n_s,double tau_reio, double r_tensor)
    {
      using namespace Pipes::function_smashInf_rLCDMtensor_lowp_TT_loglike;

      struct Class::precision pr;        /* for precision parameters */
      struct Class::background ba;       /* for cosmological background */
      struct Class::thermo th;           /* for thermodynamics */
      struct Class::perturbs pt;         /* for source functions */
      struct Class::transfers tr;        /* for transfer functions */
      struct Class::primordial pm;       /* for primordial spectra */
      struct Class::spectra sp;          /* for output spectra */
      struct Class::nonlinear nl;        /* for non-linear spectra */
      struct Class::lensing le;          /* for lensed spectra */
      struct Class::output op;           /* for output files */
      Class::ErrorMsg class_errmsg;      /* for error messages */

      int l,l_max;
      int num_ct_max=7;

      l_max=3000;

      char *class_null;// defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
      could add another array with P(k), or extract other results from
      the code - here we assume that we are interested in the C_l's
      only */

      double* cl[l_max];
      for(int i = 0; i < l_max; ++i)
        cl[i] = new double[num_ct_max];

      struct Class::file_content fc;

      BEreq::class_parser_initialize(&fc,11,"",class_errmsg);

      strcpy(fc.name[0],"output");
      strcpy(fc.name[1],"omega_b");
      strcpy(fc.name[2],"omega_cdm");
      strcpy(fc.name[3],"H0");
      strcpy(fc.name[4],"A_s");
      strcpy(fc.name[5],"n_s");
      strcpy(fc.name[6],"tau_reio");
      strcpy(fc.name[7],"modes");
      strcpy(fc.name[8],"l_scalar_max");
      strcpy(fc.name[9],"r");

      strcpy(fc.value[0],"tCl,pCl");
      sprintf(fc.value[1],"%e",omega_b);
      sprintf(fc.value[2],"%e",omega_cdm);
      sprintf(fc.value[3],"%e",H0);
      sprintf(fc.value[4],"%e",A_s);
      sprintf(fc.value[5],"%e",n_s);
      sprintf(fc.value[6],"%e",tau_reio);
      sprintf(fc.value[7],"s,t");
      sprintf(fc.value[8],"%d",3000);
      sprintf(fc.value[9],"%e",r_tensor);

      std::cout << "Here are the theory parameter values I generate the Cl's with" << std::endl;
      std::cout << "omega_b = " << omega_b << std::endl;
      std::cout << "omega_cdm = " << omega_cdm << std::endl;
      std::cout << "H0 = " << H0 << std::endl;
      std::cout << "A_s = " << A_s << std::endl;
      std::cout << "n_s = " << n_s << std::endl;
      std::cout << "tau_reio = " << tau_reio << std::endl;
      std::cout << "r_tensor = " << r_tensor << std::endl;

      BEreq::class_input_initialize(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,class_errmsg);
      BEreq::class_background_initialize(&pr,&ba);

      /* for bbn */
      sprintf(pr.sBBN_file,"");
      strcat(pr.sBBN_file,"/Users/selimhotinli/Dropbox/gambit/Backends/installed/class/2.6.1/bbn/sBBN.dat");

      BEreq::class_thermodynamics_initialize(&pr,&ba,&th);
      BEreq::class_perturb_initialize(&pr,&ba,&th,&pt);
      BEreq::class_primordial_initialize(&pr,&pt,&pm);
      BEreq::class_nonlinear_initialize(&pr,&ba,&th,&pt,&pm,&nl);
      BEreq::class_transfer_initialize(&pr,&ba,&th,&pt,&nl,&tr);
      BEreq::class_spectra_initialize(&pr,&ba,&pt,&pm,&nl,&tr,&sp);
      BEreq::class_lensing_initialize(&pr,&pt,&sp,&nl,&le);

      /****** write the Cl values in the input array cl[l]  *******/

      for (l=2; l < l_max; l++)
      {
        int errval = BEreq::class_output_total_cl_at_l(&sp,&le,&op,byVal(l),byVal(cl[l]));

        cl[l][sp.index_ct_tt] = cl[l][sp.index_ct_tt]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_te] = cl[l][sp.index_ct_te]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_ee] = cl[l][sp.index_ct_ee]*pow(ba.T_cmb*1.e6,2);
        cl[l][sp.index_ct_bb] = cl[l][sp.index_ct_bb]*pow(ba.T_cmb*1.e6,2);
      }

      cout << "we are outside the loop for multipoles " << endl;

      BEreq::class_lensing_free(&le);
      BEreq::class_spectra_free(&sp);
      BEreq::class_transfer_free(&tr);
      BEreq::class_nonlinear_free(&nl);
      BEreq::class_primordial_free(&pm);
      BEreq::class_perturb_free(&pt);
      BEreq::class_thermodynamics_free(&th);
      BEreq::class_background_free(&ba);

      return cl;
    }

    void function_vanilla_lowp_TT_loglike(double& result)
    {
      using namespace Pipes::function_vanilla_lowp_TT_loglike;

      int l,l_max;
      l_max=3000;
      //  int nmax;
      int num_ct_max=7;

      double  highl_TT_cl_and_pars[2525];
      double  lowp_cl_and_pars[121];

      char *class_null;  // defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
       could add another array with P(k), or extract other results from
       the code - here we assume that we are interested in the C_l's
       only */

      double** cl;

      cl = return_vanilla_cls(*Param["omega_b"],
                  *Param["omega_cdm"],
                  *Param["H0"],
                  *Param["ln10A_s"],
                  *Param["n_s"],
                  *Param["tau_reio"]);
      //--------------------------------------------------------------------------
      //------high-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------

      highl_TT_cl_and_pars[0] = 0.0;
      highl_TT_cl_and_pars[1] = 0.0;

      for(int ii = 2; ii < 2509 ; ii++)highl_TT_cl_and_pars[ii] = cl[ii][0];

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      highl_TT_cl_and_pars[2509] = *Param["A_cib_217"];
      highl_TT_cl_and_pars[2510] = *Param["cib_index"];
      highl_TT_cl_and_pars[2511] = *Param["xi_sz_cib"];
      highl_TT_cl_and_pars[2512] = *Param["A_sz"];
      highl_TT_cl_and_pars[2513] = *Param["ps_A_100_100"];
      highl_TT_cl_and_pars[2514] = *Param["ps_A_143_143"];
      highl_TT_cl_and_pars[2515] = *Param["ps_A_143_217"];
      highl_TT_cl_and_pars[2516] = *Param["ps_A_217_217"];
      highl_TT_cl_and_pars[2517] = *Param["ksz_norm"];
      highl_TT_cl_and_pars[2518] = *Param["gal545_A_100"];
      highl_TT_cl_and_pars[2519] = *Param["gal545_A_143"];
      highl_TT_cl_and_pars[2520] = *Param["gal545_A_143_217"];
      highl_TT_cl_and_pars[2521] = *Param["gal545_A_217"];
      highl_TT_cl_and_pars[2522] = *Param["calib_100T"];
      highl_TT_cl_and_pars[2523] = *Param["calib_217T"];
      highl_TT_cl_and_pars[2524] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //-------low-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[0] = 0.0;
      lowp_cl_and_pars[1] = 0.0;
      for (l=2;l<=29;l++)  {
        lowp_cl_and_pars[l] = cl[l][0];
      }
      int k;
      lowp_cl_and_pars[30] = 0.0;
      lowp_cl_and_pars[31] = 0.0;
      for (l=32;l<=59;l++)  {
        k = l-30;
        lowp_cl_and_pars[l] = cl[k][1];
      }
      lowp_cl_and_pars[60] = 0.0;
      lowp_cl_and_pars[61] = 0.0;
      for (l=62;l<=89;l++)  {
        k = l-60;
        //lowp_cl_and_pars[l] = cl[k][2];
        lowp_cl_and_pars[l] = cl[k][3];
      }
      lowp_cl_and_pars[90] = 0.0;
      lowp_cl_and_pars[91] = 0.0;
      for (l=92;l<=119;l++)  {
        k = l-90;
        //lowp_cl_and_pars[l] = cl[k][3];
        lowp_cl_and_pars[l] = cl[k][2];
      }
      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[120] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_object* high_clikid;
      clik_object* lowl_clikid;
      clik_error *_err;
      double lowp_log;
      double highl_TT_log;

      lowl_clikid = BEreq::return_lowp_TT();
      high_clikid = BEreq::return_high_TT();

      _err = BEreq::clik_initialize_error();

      highl_TT_log = BEreq::clik_compute_loglike(byVal(high_clikid),
                             byVal(highl_TT_cl_and_pars),
                             &_err);
      lowp_log    = BEreq::clik_compute_loglike(byVal(lowl_clikid),
                             byVal(lowp_cl_and_pars),
                             &_err);


      result = highl_TT_log+lowp_log;

      std::cout << "Log likelihood is : " << result << std::endl;
    }


    void function_Planck_high_TT_loglike(double& result)
    {
      //std::cout << "Last seen alive in: function_Planck_high_TT_loglike" << std::endl;
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

      //std::cout << "Log likelihood (of high_TT) is : " << result << std::endl;
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

      //std::cout << "Log likelihood (of high_TTTEEE) is : " << result << std::endl;
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

      //std::cout << "Log likelihood (of high_TT_lite) is : " << result << std::endl;
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

      //std::cout << "Log likelihood (of lensing) is : " << result << std::endl;
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

      //std::cout << "Log likelihood (of lowp_TT) is : " << result << std::endl;
    }

    void function_LCDMtensor_lowp_TT_loglike(double& result)
    {
      using namespace Pipes::function_LCDMtensor_lowp_TT_loglike;

      int l,l_max;
      l_max=3000;
      //  int nmax;
      int num_ct_max=7;

      double  highl_TT_cl_and_pars[2525];
      double  lowp_cl_and_pars[121];

      char *class_null;  // defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
       could add another array with P(k), or extract other results from
       the code - here we assume that we are interested in the C_l's
       only */

      double** cl;

      cl = return_LCDMtensor_cls(*Param["omega_b"],
                     *Param["omega_cdm"],
                     *Param["H0"],
                     *Param["ln10A_s"],
                     *Param["n_s"],
                     *Param["tau_reio"],
                     *Param["r_tensor"]);
      //--------------------------------------------------------------------------
      //------high-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------

      highl_TT_cl_and_pars[0] = 0.0;
      highl_TT_cl_and_pars[1] = 0.0;

      for(int ii = 2; ii < 2509 ; ii++)highl_TT_cl_and_pars[ii] = cl[ii][0];

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      highl_TT_cl_and_pars[2509] = *Param["A_cib_217"];
      highl_TT_cl_and_pars[2510] = *Param["cib_index"];
      highl_TT_cl_and_pars[2511] = *Param["xi_sz_cib"];
      highl_TT_cl_and_pars[2512] = *Param["A_sz"];
      highl_TT_cl_and_pars[2513] = *Param["ps_A_100_100"];
      highl_TT_cl_and_pars[2514] = *Param["ps_A_143_143"];
      highl_TT_cl_and_pars[2515] = *Param["ps_A_143_217"];
      highl_TT_cl_and_pars[2516] = *Param["ps_A_217_217"];
      highl_TT_cl_and_pars[2517] = *Param["ksz_norm"];
      highl_TT_cl_and_pars[2518] = *Param["gal545_A_100"];
      highl_TT_cl_and_pars[2519] = *Param["gal545_A_143"];
      highl_TT_cl_and_pars[2520] = *Param["gal545_A_143_217"];
      highl_TT_cl_and_pars[2521] = *Param["gal545_A_217"];
      highl_TT_cl_and_pars[2522] = *Param["calib_100T"];
      highl_TT_cl_and_pars[2523] = *Param["calib_217T"];
      highl_TT_cl_and_pars[2524] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //-------low-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[0] = 0.0;
      lowp_cl_and_pars[1] = 0.0;
      for (l=2;l<=29;l++)  {
        lowp_cl_and_pars[l] = cl[l][0];
      }
      int k;
      lowp_cl_and_pars[30] = 0.0;
      lowp_cl_and_pars[31] = 0.0;
      for (l=32;l<=59;l++)  {
        k = l-30;
        lowp_cl_and_pars[l] = cl[k][1];
      }
      lowp_cl_and_pars[60] = 0.0;
      lowp_cl_and_pars[61] = 0.0;
      for (l=62;l<=89;l++)  {
        k = l-60;
        //lowp_cl_and_pars[l] = cl[k][2];
        lowp_cl_and_pars[l] = cl[k][3];
      }
      lowp_cl_and_pars[90] = 0.0;
      lowp_cl_and_pars[91] = 0.0;
      for (l=92;l<=119;l++)  {
        k = l-90;
        //lowp_cl_and_pars[l] = cl[k][3];
        lowp_cl_and_pars[l] = cl[k][2];
      }
      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[120] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_object* high_clikid;
      clik_object* lowl_clikid;
      clik_error *_err;
      double lowp_log;
      double highl_TT_log;

      lowl_clikid = BEreq::return_lowp_TT();
      high_clikid = BEreq::return_high_TT();

      _err = BEreq::clik_initialize_error();

      highl_TT_log = BEreq::clik_compute_loglike(byVal(high_clikid),
             byVal(highl_TT_cl_and_pars),
             &_err);
      lowp_log    = BEreq::clik_compute_loglike(byVal(lowl_clikid),
            byVal(lowp_cl_and_pars),
            &_err);


      result = highl_TT_log+lowp_log;

      std::cout << "Log likelihood is : " << result << std::endl;
    }

    // declare an array and size from the yaml. doesn't have to be const. reading it from an array and make it static - static int - first time created -
    // instead of std vectors .
    void function_LCDMtensor_inflation_lowp_TT_loglike(double& result)
    {
      using namespace Pipes::function_LCDMtensor_inflation_lowp_TT_loglike;

      // Initialization parameters controlling main characteristics.
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");


      // Control the output of analytic approximations for comparison.
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");

      // Parameters to control how the ICs are sampled.
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");

      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this

      // Parameters to control how the vparams are sampled.
      int param_sampling = runOptions->getValue<int> ("param_sampling");


      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");

      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");

      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");

      // Parameters to be passed to the potential
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

      double N_pivot = runOptions->getValue<double> ("N_pivot");
      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");


      // Priors on the IC and N_pivot ranges
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");

      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");

      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");

      gambit_inflation_observables observs;

      std::cout << "num_inflaton = " << num_inflaton << std::endl;

      std::cout << "save_iso_N = " << save_iso_N << std::endl;

      std::cout << "N_pivot = " << N_pivot << std::endl;

      std::cout << "kmin = " << kmin << std::endl;

      std::cout << "dlnk = " << dlnk << std::endl;


      // The function below calls the multimodecode backend for a given choice of inflationary model,
      // which calculates the observables.
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
                         byVal(&dphi_init0[0]), // dphi correct at the multimodecode
                         byVal(&vparams[0]),
                         N_pivot,
                         k_pivot,
                         dlnk,
                         calc_full_pk,
                         steps,
                         kmin,
//                         kmax,
                         byVal(&phi0_priors_min[0]),
                         byVal(&phi0_priors_max[0]),
                         byVal(&dphi0_priors_min[0]),
                         byVal(&dphi0_priors_max[0]),
                         N_pivot_prior_min,
                         N_pivot_prior_max);


      std::cout << "observs.As " << observs.As << std::endl;
      std::cout << "observs.A_iso " << observs.A_iso << std::endl;
      std::cout << "observs.A_pnad " << observs.A_pnad << std::endl;
      std::cout << "observs.A_ent " << observs.A_ent << std::endl;
      std::cout << "observs.A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
//      std::cout << "observs.A_bundle " << observs.A_bundle << std::endl;
      std::cout << "observs.ns " << observs.ns << std::endl;
      std::cout << "observs.nt " << observs.nt << std::endl;
      std::cout << "observs.n_iso " << observs.n_iso << std::endl;
      std::cout << "observs.n_pnad " << observs.n_pnad << std::endl;
      std::cout << "observs.n_ent " << observs.n_ent << std::endl;
      std::cout << "observs.r " << observs.r << std::endl;
      std::cout << "observs.alpha_s " << observs.alpha_s << std::endl;
      std::cout << "observs.runofrun " << observs.runofrun << std::endl;
      std::cout << "observs.f_NL " << observs.f_NL << std::endl;
      std::cout << "observs.tau_NL " << observs.tau_NL << std::endl;

//*/
      int l,l_max;
      l_max=3000;
      //  int nmax;
      int num_ct_max=7;

      double  highl_TT_cl_and_pars[2525];
      double  lowp_cl_and_pars[121];

      char *class_null;  // defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
       could add another array with P(k), or extract other results from
       the code - here we assume that we are interested in the C_l's
       only */

      double** cl;

      cl = return_LCDMtensor_inflation_cls(*Param["omega_b"],
                         *Param["omega_cdm"],
                         *Param["H0"],
                         *Param["ln10A_s"],
//                         *Param["n_s"],
                         observs.ns,
                         *Param["tau_reio"],
//                         *Param["r_tensor"]
                         observs.r);
      //--------------------------------------------------------------------------
      //------high-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------

      highl_TT_cl_and_pars[0] = 0.0;
      highl_TT_cl_and_pars[1] = 0.0;

      for(int ii = 2; ii < 2509 ; ii++)highl_TT_cl_and_pars[ii] = cl[ii][0];

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      highl_TT_cl_and_pars[2509] = *Param["A_cib_217"];
      highl_TT_cl_and_pars[2510] = *Param["cib_index"];
      highl_TT_cl_and_pars[2511] = *Param["xi_sz_cib"];
      highl_TT_cl_and_pars[2512] = *Param["A_sz"];
      highl_TT_cl_and_pars[2513] = *Param["ps_A_100_100"];
      highl_TT_cl_and_pars[2514] = *Param["ps_A_143_143"];
      highl_TT_cl_and_pars[2515] = *Param["ps_A_143_217"];
      highl_TT_cl_and_pars[2516] = *Param["ps_A_217_217"];
      highl_TT_cl_and_pars[2517] = *Param["ksz_norm"];
      highl_TT_cl_and_pars[2518] = *Param["gal545_A_100"];
      highl_TT_cl_and_pars[2519] = *Param["gal545_A_143"];
      highl_TT_cl_and_pars[2520] = *Param["gal545_A_143_217"];
      highl_TT_cl_and_pars[2521] = *Param["gal545_A_217"];
      highl_TT_cl_and_pars[2522] = *Param["calib_100T"];
      highl_TT_cl_and_pars[2523] = *Param["calib_217T"];
      highl_TT_cl_and_pars[2524] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //-------low-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[0] = 0.0;
      lowp_cl_and_pars[1] = 0.0;
      for (l=2;l<=29;l++)  {
        lowp_cl_and_pars[l] = cl[l][0];
      }
      int k;
      lowp_cl_and_pars[30] = 0.0;
      lowp_cl_and_pars[31] = 0.0;
      for (l=32;l<=59;l++)  {
        k = l-30;
        lowp_cl_and_pars[l] = cl[k][1];
      }
      lowp_cl_and_pars[60] = 0.0;
      lowp_cl_and_pars[61] = 0.0;
      for (l=62;l<=89;l++)  {
        k = l-60;
        lowp_cl_and_pars[l] = cl[k][2];
      }
      lowp_cl_and_pars[90] = 0.0;
      lowp_cl_and_pars[91] = 0.0;
      for (l=92;l<=119;l++)  {
        k = l-90;
        lowp_cl_and_pars[l] = cl[k][3];
      }
      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[120] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_object* high_clikid;
      clik_object* lowl_clikid;
      clik_error *_err;
      double lowp_log;
      double highl_TT_log;

      lowl_clikid = BEreq::return_lowp_TT();
      high_clikid = BEreq::return_high_TT();

      _err = BEreq::clik_initialize_error();

      highl_TT_log = BEreq::clik_compute_loglike(byVal(high_clikid),
                             byVal(highl_TT_cl_and_pars),
                             &_err);
      lowp_log    = BEreq::clik_compute_loglike(byVal(lowl_clikid),
                             byVal(lowp_cl_and_pars),
                             &_err);


      result = highl_TT_log+lowp_log;

      std::cout << "Log likelihood is : " << result << std::endl;



    }


  // Simplest quadratic inflationary model with 7 params with Planck log-like.
  void function_1quadInf_rLCDMtensor_lowp_TT_loglike(double& result)
  {
    using namespace Pipes::function_1quadInf_rLCDMtensor_lowp_TT_loglike;

    // Initialization parameters controlling main characteristics.
    int num_inflaton = runOptions->getValue<int> ("num_inflaton");
    int potential_choice = runOptions->getValue<int> ("potential_choice");
    int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
    int instreheat = runOptions->getValue<int> ("instreheat");
    int vparam_rows = runOptions->getValue<int> ("vparam_rows");


    // Control the output of analytic approximations for comparison.
    int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
    int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
    int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
    int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");

    // Parameters to control how the ICs are sampled.
    int ic_sampling = runOptions->getValue<int> ("ic_sampling");
    double energy_scale = runOptions->getValue<double> ("energy_scale");
    int numb_samples = runOptions->getValue<int> ("numb_samples");
    int save_iso_N = runOptions->getValue<int> ("save_iso_N");

    double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this

    // Parameters to control how the vparams are sampled.
    int param_sampling = runOptions->getValue<int> ("param_sampling");


    std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
    std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");

    int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
    int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");

    std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
    std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");

    // Parameters to be passed to the potential
    std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

    //-------------------------------------------------
    //----------------- Sampling m^2 ------------------
    //-------------------------------------------------
    vparams[0] = *Param["m2_inflaton"];

    //-------------------------------------------------
    //--------------- Sampling N_pivot ----------------
    //-------------------------------------------------
    double N_pivot = *Param["N_pivot"];

    double k_pivot = runOptions->getValue<double> ("k_pivot");
    double dlnk = runOptions->getValue<double> ("dlnk");


    // Priors on the IC and N_pivot ranges
    std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
    std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
    std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
    std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");

    double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
    double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");

    // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
    int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
    int steps = runOptions->getValue<int> ("steps");
    double kmin = runOptions->getValue<double> ("kmin");
    double kmax = runOptions->getValue<double> ("kmax");

    gambit_inflation_observables observs;

    std::cout << "steps = " << steps << std::endl;

//    if (calc_full_pk == 1){
//
//      std::vector<double> vec1(steps,0.0);
//
//      observs.k_array   = vec1;
//      observs.pks_array = vec1;
//      observs.pkt_array = vec1;
//    }

    std::cout << "num_inflaton = " << num_inflaton << std::endl;

    std::cout << "save_iso_N = " << save_iso_N << std::endl;

    std::cout << "N_pivot = " << N_pivot << std::endl;

    std::cout << "kmin = " << kmin << std::endl;

    std::cout << "dlnk = " << dlnk << std::endl;

    std::cout << "m2_inflaton = " << vparams[0] << std::endl;

    std::cout << "observs.k_array = " << observs.k_array << std::endl;
    std::cout << "observs.pks_array = " << observs.pks_array << std::endl;
    std::cout << "observs.pkt_array = " << observs.pkt_array << std::endl;

    double** cl;

    /** - Make room */
//    observs.k_array   = (double *) malloc(steps*sizeof(double));
//    observs.pks_array = (double *) malloc(steps*sizeof(double));
//    observs.pkt_array = (double *) malloc(steps*sizeof(double));

    // The function below calls the multimodecode backend for a given choice of inflationary model,
    // which calculates the observables.
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
                       byVal(&dphi_init0[0]), // dphi correct at the multimodecode
                       byVal(&vparams[0]),
                       N_pivot,
                       k_pivot,
                       dlnk,
                       calc_full_pk,
                       steps,
                       kmin,
                       //kmax,
                       byVal(&phi0_priors_min[0]),
                       byVal(&phi0_priors_max[0]),
                       byVal(&dphi0_priors_min[0]),
                       byVal(&dphi0_priors_max[0]),
                       N_pivot_prior_min,
                       N_pivot_prior_max);

    std::cout << " we are out of the multimodecode_gambit_driver function" << std::endl;

    std::cout << "observs.k_array = " << sizeof(observs.k_array) << std::endl;
    std::cout << "observs.k_array[0] = " << observs.k_array[0] << std::endl;
//    std::cout << "observs.pks_array = " << observs.pks_array << std::endl;
//    std::cout << "observs.pkt_array = " << observs.pkt_array << std::endl;


    if (calc_full_pk == 0)
    {
    std::cout << "observs.As " << observs.As << std::endl;
    std::cout << "observs.A_iso " << observs.A_iso << std::endl;
    std::cout << "observs.A_pnad " << observs.A_pnad << std::endl;
    std::cout << "observs.A_ent " << observs.A_ent << std::endl;
    std::cout << "observs.A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
    //std::cout << "observs.A_bundle " << observs.A_bundle << std::endl;
    std::cout << "observs.ns " << observs.ns << std::endl;
    std::cout << "observs.nt " << observs.nt << std::endl;
    std::cout << "observs.n_iso " << observs.n_iso << std::endl;
    std::cout << "observs.n_pnad " << observs.n_pnad << std::endl;
    std::cout << "observs.n_ent " << observs.n_ent << std::endl;
    std::cout << "observs.r " << observs.r << std::endl;
    std::cout << "observs.alpha_s " << observs.alpha_s << std::endl;
    std::cout << "observs.runofrun " << observs.runofrun << std::endl;
    std::cout << "observs.f_NL " << observs.f_NL << std::endl;
    std::cout << "observs.tau_NL " << observs.tau_NL << std::endl;

    cl = return_LCDMtensor_1quadInf_cls(*Param["omega_b"],
                       *Param["omega_cdm"],
                       *Param["H0"],
                       observs.As,
                       observs.ns,
                       *Param["tau_reio"],
                       observs.r);
    }

    else
    {
     cl = return_LCDMtensor_1quadInf_cls_fullPk(*Param["omega_b"],
                          *Param["omega_cdm"],
                          *Param["H0"],
                          *Param["tau_reio"],
                          observs.k_array,
                          observs.pks_array,
                          observs.pkt_array,
                          steps);
    }


    int l,l_max;
    l_max=3000;
    //  int nmax;
    int num_ct_max=7;

    double  highl_TT_cl_and_pars[2525];
    double  lowp_cl_and_pars[121];

    char *class_null;  // defining an empty string for parser initialization.

    /* allocate the array where calculated Cl's will be written (we
     could add another array with P(k), or extract other results from
     the code - here we assume that we are interested in the C_l's
     only */


    //--------------------------------------------------------------------------
    //------high-l likelihood calculation making of Cls-------------------------
    //--------------------------------------------------------------------------

    highl_TT_cl_and_pars[0] = 0.0;
    highl_TT_cl_and_pars[1] = 0.0;

    for(int ii = 2; ii < 2509 ; ii++)highl_TT_cl_and_pars[ii] = cl[ii][0];

    //--------------------------------------------------------------------------
    //------addition of nuisance parameters to Cl array-------------------------
    //--------------------------------------------------------------------------
    highl_TT_cl_and_pars[2509] = *Param["A_cib_217"];
    highl_TT_cl_and_pars[2510] = *Param["cib_index"];
    highl_TT_cl_and_pars[2511] = *Param["xi_sz_cib"];
    highl_TT_cl_and_pars[2512] = *Param["A_sz"];
    highl_TT_cl_and_pars[2513] = *Param["ps_A_100_100"];
    highl_TT_cl_and_pars[2514] = *Param["ps_A_143_143"];
    highl_TT_cl_and_pars[2515] = *Param["ps_A_143_217"];
    highl_TT_cl_and_pars[2516] = *Param["ps_A_217_217"];
    highl_TT_cl_and_pars[2517] = *Param["ksz_norm"];
    highl_TT_cl_and_pars[2518] = *Param["gal545_A_100"];
    highl_TT_cl_and_pars[2519] = *Param["gal545_A_143"];
    highl_TT_cl_and_pars[2520] = *Param["gal545_A_143_217"];
    highl_TT_cl_and_pars[2521] = *Param["gal545_A_217"];
    highl_TT_cl_and_pars[2522] = *Param["calib_100T"];
    highl_TT_cl_and_pars[2523] = *Param["calib_217T"];
    highl_TT_cl_and_pars[2524] = *Param["A_planck"];

    //--------------------------------------------------------------------------
    //-------low-l likelihood calculation making of Cls-------------------------
    //--------------------------------------------------------------------------
    lowp_cl_and_pars[0] = 0.0;
    lowp_cl_and_pars[1] = 0.0;
    for (l=2;l<=29;l++)  {
      lowp_cl_and_pars[l] = cl[l][0];
    }
    int k;
    lowp_cl_and_pars[30] = 0.0;
    lowp_cl_and_pars[31] = 0.0;
    for (l=32;l<=59;l++)  {
      k = l-30;
      lowp_cl_and_pars[l] = cl[k][1];
    }
    lowp_cl_and_pars[60] = 0.0;
    lowp_cl_and_pars[61] = 0.0;
    for (l=62;l<=89;l++)  {
      k = l-60;
      lowp_cl_and_pars[l] = cl[k][2];
    }
    lowp_cl_and_pars[90] = 0.0;
    lowp_cl_and_pars[91] = 0.0;
    for (l=92;l<=119;l++)  {
      k = l-90;
      lowp_cl_and_pars[l] = cl[k][3];
    }
    //--------------------------------------------------------------------------
    //------addition of nuisance parameters to Cl array-------------------------
    //--------------------------------------------------------------------------
    lowp_cl_and_pars[120] = *Param["A_planck"];

    //--------------------------------------------------------------------------
    //------calculation of the planck loglikelihood-----------------------------
    //--------------------------------------------------------------------------
    clik_object* high_clikid;
    clik_object* lowl_clikid;
    clik_error *_err;
    double lowp_log;
    double highl_TT_log;

    lowl_clikid = BEreq::return_lowp_TT();
    high_clikid = BEreq::return_high_TT();

    _err = BEreq::clik_initialize_error();

    highl_TT_log = BEreq::clik_compute_loglike(byVal(high_clikid),
                           byVal(highl_TT_cl_and_pars),
                           &_err);
    lowp_log    = BEreq::clik_compute_loglike(byVal(lowl_clikid),
                           byVal(lowp_cl_and_pars),
                           &_err);


    result = highl_TT_log+lowp_log;

    std::cout << "Log likelihood is : " << result << std::endl;

  }

    // Simplest quartic inflationary model with 7 params with Planck log-like.
    void function_1quarInf_rLCDMtensor_lowp_TT_loglike(double& result)
    {
      using namespace Pipes::function_1quarInf_rLCDMtensor_lowp_TT_loglike;

      // Initialization parameters controlling main characteristics.
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");

      // Control the output of analytic approximations for comparison.
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");

      // Parameters to control how the ICs are sampled.
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");

      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this

      // Parameters to control how the vparams are sampled.
      int param_sampling = runOptions->getValue<int> ("param_sampling");


      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");

      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");

      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");

      // Parameters to be passed to the potential
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

      //-------------------------------------------------
      //----------------- Sampling \lambda ------------------
      //-------------------------------------------------
      vparams[0] = *Param["lambda"];

      //-------------------------------------------------
      //--------------- Sampling N_pivot ----------------
      //-------------------------------------------------
      double N_pivot = *Param["N_pivot"];

      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");


      // Priors on the IC and N_pivot ranges
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");

      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");

      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");

      gambit_inflation_observables observs;

      std::cout << "num_inflaton = " << num_inflaton << std::endl;

      std::cout << "save_iso_N = " << save_iso_N << std::endl;

      std::cout << "N_pivot = " << N_pivot << std::endl;

      std::cout << "kmin = " << kmin << std::endl;

      std::cout << "dlnk = " << dlnk << std::endl;

      std::cout << "lambda = " << vparams[0] << std::endl;

      // The function below calls the multimodecode backend for a given choice of inflationary model,
      // which calculates the observables.
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
                         byVal(&dphi_init0[0]), // dphi correct at the multimodecode
                         byVal(&vparams[0]),
                         N_pivot,
                         k_pivot,
                         dlnk,
                         calc_full_pk,
                         steps,
                         kmin,
                         //                         kmax,
                         byVal(&phi0_priors_min[0]),
                         byVal(&phi0_priors_max[0]),
                         byVal(&dphi0_priors_min[0]),
                         byVal(&dphi0_priors_max[0]),
                         N_pivot_prior_min,
                         N_pivot_prior_max);


      std::cout << "observs.As " << observs.As << std::endl;
      std::cout << "observs.A_iso " << observs.A_iso << std::endl;
      std::cout << "observs.A_pnad " << observs.A_pnad << std::endl;
      std::cout << "observs.A_ent " << observs.A_ent << std::endl;
      std::cout << "observs.A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
      //std::cout << "observs.A_bundle " << observs.A_bundle << std::endl;
      std::cout << "observs.ns " << observs.ns << std::endl;
      std::cout << "observs.nt " << observs.nt << std::endl;
      std::cout << "observs.n_iso " << observs.n_iso << std::endl;
      std::cout << "observs.n_pnad " << observs.n_pnad << std::endl;
      std::cout << "observs.n_ent " << observs.n_ent << std::endl;
      std::cout << "observs.r " << observs.r << std::endl;
      std::cout << "observs.alpha_s " << observs.alpha_s << std::endl;
      std::cout << "observs.runofrun " << observs.runofrun << std::endl;
      std::cout << "observs.f_NL " << observs.f_NL << std::endl;
      std::cout << "observs.tau_NL " << observs.tau_NL << std::endl;

      //*/
      int l,l_max;
      l_max=3000;
      //  int nmax;
      int num_ct_max=7;

      double  highl_TT_cl_and_pars[2525];
      double  lowp_cl_and_pars[121];

      char *class_null;  // defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
       could add another array with P(k), or extract other results from
       the code - here we assume that we are interested in the C_l's
       only */

      double** cl;

      cl = return_LCDMtensor_1quarInf_cls(*Param["omega_b"],
                       *Param["omega_cdm"],
                       *Param["H0"],
                       //*Param["ln10A_s"],
                       observs.As,
                       //*Param["n_s"],
                       observs.ns,
                       *Param["tau_reio"],
                       //*Param["r_tensor"]
                       observs.r);
      //--------------------------------------------------------------------------
      //------high-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------

      highl_TT_cl_and_pars[0] = 0.0;
      highl_TT_cl_and_pars[1] = 0.0;

      for(int ii = 2; ii < 2509 ; ii++)highl_TT_cl_and_pars[ii] = cl[ii][0];

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      highl_TT_cl_and_pars[2509] = *Param["A_cib_217"];
      highl_TT_cl_and_pars[2510] = *Param["cib_index"];
      highl_TT_cl_and_pars[2511] = *Param["xi_sz_cib"];
      highl_TT_cl_and_pars[2512] = *Param["A_sz"];
      highl_TT_cl_and_pars[2513] = *Param["ps_A_100_100"];
      highl_TT_cl_and_pars[2514] = *Param["ps_A_143_143"];
      highl_TT_cl_and_pars[2515] = *Param["ps_A_143_217"];
      highl_TT_cl_and_pars[2516] = *Param["ps_A_217_217"];
      highl_TT_cl_and_pars[2517] = *Param["ksz_norm"];
      highl_TT_cl_and_pars[2518] = *Param["gal545_A_100"];
      highl_TT_cl_and_pars[2519] = *Param["gal545_A_143"];
      highl_TT_cl_and_pars[2520] = *Param["gal545_A_143_217"];
      highl_TT_cl_and_pars[2521] = *Param["gal545_A_217"];
      highl_TT_cl_and_pars[2522] = *Param["calib_100T"];
      highl_TT_cl_and_pars[2523] = *Param["calib_217T"];
      highl_TT_cl_and_pars[2524] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //-------low-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[0] = 0.0;
      lowp_cl_and_pars[1] = 0.0;
      for (l=2;l<=29;l++)  {
        lowp_cl_and_pars[l] = cl[l][0];
      }
      int k;
      lowp_cl_and_pars[30] = 0.0;
      lowp_cl_and_pars[31] = 0.0;
      for (l=32;l<=59;l++)  {
        k = l-30;
        lowp_cl_and_pars[l] = cl[k][1];
      }
      lowp_cl_and_pars[60] = 0.0;
      lowp_cl_and_pars[61] = 0.0;
      for (l=62;l<=89;l++)  {
        k = l-60;
        lowp_cl_and_pars[l] = cl[k][2];
      }
      lowp_cl_and_pars[90] = 0.0;
      lowp_cl_and_pars[91] = 0.0;
      for (l=92;l<=119;l++)  {
        k = l-90;
        lowp_cl_and_pars[l] = cl[k][3];
      }
      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[120] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_object* high_clikid;
      clik_object* lowl_clikid;
      clik_error *_err;
      double lowp_log;
      double highl_TT_log;

      lowl_clikid = BEreq::return_lowp_TT();
      high_clikid = BEreq::return_high_TT();

      _err = BEreq::clik_initialize_error();

      highl_TT_log = BEreq::clik_compute_loglike(byVal(high_clikid),
                             byVal(highl_TT_cl_and_pars),
                             &_err);
      lowp_log    = BEreq::clik_compute_loglike(byVal(lowl_clikid),
                             byVal(lowp_cl_and_pars),
                             &_err);


      result = highl_TT_log+lowp_log;

      std::cout << "Log likelihood is : " << result << std::endl;



    }


    // Simplest mono32 potential inflationary model with 7 params with Planck log-like.
    void function_1mono32Inf_rLCDMtensor_lowp_TT_loglike(double& result)
    {
      using namespace Pipes::function_1mono32Inf_rLCDMtensor_lowp_TT_loglike;

      // Initialization parameters controlling main characteristics.
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");

      // Control the output of analytic approximations for comparison.
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");

      // Parameters to control how the ICs are sampled.
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");

      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this

      // Parameters to control how the vparams are sampled.
      int param_sampling = runOptions->getValue<int> ("param_sampling");


      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");

      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");

      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");

      // Parameters to be passed to the potential
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

      //-------------------------------------------------
      //----------------- Sampling \lambda ------------------
      //-------------------------------------------------
      vparams[0] = *Param["lambda"];

      //-------------------------------------------------
      //--------------- Sampling N_pivot ----------------
      //-------------------------------------------------
      double N_pivot = *Param["N_pivot"];

      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");


      // Priors on the IC and N_pivot ranges
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");

      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");

      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");

      gambit_inflation_observables observs;

      std::cout << "num_inflaton = " << num_inflaton << std::endl;

      std::cout << "save_iso_N = " << save_iso_N << std::endl;

      std::cout << "N_pivot = " << N_pivot << std::endl;

      std::cout << "kmin = " << kmin << std::endl;

      std::cout << "dlnk = " << dlnk << std::endl;

      std::cout << "lambda = " << vparams[0] << std::endl;

      // The function below calls the multimodecode backend for a given choice of inflationary model,
      // which calculates the observables.
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
                         byVal(&dphi_init0[0]), // dphi correct at the multimodecode
                         byVal(&vparams[0]),
                         N_pivot,
                         k_pivot,
                         dlnk,
                         calc_full_pk,
                         steps,
                         kmin,
                         //                         kmax,
                         byVal(&phi0_priors_min[0]),
                         byVal(&phi0_priors_max[0]),
                         byVal(&dphi0_priors_min[0]),
                         byVal(&dphi0_priors_max[0]),
                         N_pivot_prior_min,
                         N_pivot_prior_max);


      std::cout << "observs.As " << observs.As << std::endl;
      std::cout << "observs.A_iso " << observs.A_iso << std::endl;
      std::cout << "observs.A_pnad " << observs.A_pnad << std::endl;
      std::cout << "observs.A_ent " << observs.A_ent << std::endl;
      std::cout << "observs.A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
      //std::cout << "observs.A_bundle " << observs.A_bundle << std::endl;
      std::cout << "observs.ns " << observs.ns << std::endl;
      std::cout << "observs.nt " << observs.nt << std::endl;
      std::cout << "observs.n_iso " << observs.n_iso << std::endl;
      std::cout << "observs.n_pnad " << observs.n_pnad << std::endl;
      std::cout << "observs.n_ent " << observs.n_ent << std::endl;
      std::cout << "observs.r " << observs.r << std::endl;
      std::cout << "observs.alpha_s " << observs.alpha_s << std::endl;
      std::cout << "observs.runofrun " << observs.runofrun << std::endl;
      std::cout << "observs.f_NL " << observs.f_NL << std::endl;
      std::cout << "observs.tau_NL " << observs.tau_NL << std::endl;

      //*/
      int l,l_max;
      l_max=3000;
      //  int nmax;
      int num_ct_max=7;

      double  highl_TT_cl_and_pars[2525];
      double  lowp_cl_and_pars[121];

      char *class_null;  // defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
       could add another array with P(k), or extract other results from
       the code - here we assume that we are interested in the C_l's
       only */

      double** cl;

      cl = return_LCDMtensor_1mono32Inf_cls(*Param["omega_b"],
                       *Param["omega_cdm"],
                       *Param["H0"],
                       //*Param["ln10A_s"],
                       observs.As,
                       //*Param["n_s"],
                       observs.ns,
                       *Param["tau_reio"],
                       //*Param["r_tensor"]
                       observs.r);
      //--------------------------------------------------------------------------
      //------high-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------

      highl_TT_cl_and_pars[0] = 0.0;
      highl_TT_cl_and_pars[1] = 0.0;

      for(int ii = 2; ii < 2509 ; ii++)highl_TT_cl_and_pars[ii] = cl[ii][0];

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      highl_TT_cl_and_pars[2509] = *Param["A_cib_217"];
      highl_TT_cl_and_pars[2510] = *Param["cib_index"];
      highl_TT_cl_and_pars[2511] = *Param["xi_sz_cib"];
      highl_TT_cl_and_pars[2512] = *Param["A_sz"];
      highl_TT_cl_and_pars[2513] = *Param["ps_A_100_100"];
      highl_TT_cl_and_pars[2514] = *Param["ps_A_143_143"];
      highl_TT_cl_and_pars[2515] = *Param["ps_A_143_217"];
      highl_TT_cl_and_pars[2516] = *Param["ps_A_217_217"];
      highl_TT_cl_and_pars[2517] = *Param["ksz_norm"];
      highl_TT_cl_and_pars[2518] = *Param["gal545_A_100"];
      highl_TT_cl_and_pars[2519] = *Param["gal545_A_143"];
      highl_TT_cl_and_pars[2520] = *Param["gal545_A_143_217"];
      highl_TT_cl_and_pars[2521] = *Param["gal545_A_217"];
      highl_TT_cl_and_pars[2522] = *Param["calib_100T"];
      highl_TT_cl_and_pars[2523] = *Param["calib_217T"];
      highl_TT_cl_and_pars[2524] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //-------low-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[0] = 0.0;
      lowp_cl_and_pars[1] = 0.0;
      for (l=2;l<=29;l++)  {
        lowp_cl_and_pars[l] = cl[l][0];
      }
      int k;
      lowp_cl_and_pars[30] = 0.0;
      lowp_cl_and_pars[31] = 0.0;
      for (l=32;l<=59;l++)  {
        k = l-30;
        lowp_cl_and_pars[l] = cl[k][1];
      }
      lowp_cl_and_pars[60] = 0.0;
      lowp_cl_and_pars[61] = 0.0;
      for (l=62;l<=89;l++)  {
        k = l-60;
        lowp_cl_and_pars[l] = cl[k][2];
      }
      lowp_cl_and_pars[90] = 0.0;
      lowp_cl_and_pars[91] = 0.0;
      for (l=92;l<=119;l++)  {
        k = l-90;
        lowp_cl_and_pars[l] = cl[k][3];
      }
      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[120] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_object* high_clikid;
      clik_object* lowl_clikid;
      clik_error *_err;
      double lowp_log;
      double highl_TT_log;

      lowl_clikid = BEreq::return_lowp_TT();
      high_clikid = BEreq::return_high_TT();

      _err = BEreq::clik_initialize_error();

      highl_TT_log = BEreq::clik_compute_loglike(byVal(high_clikid),
                             byVal(highl_TT_cl_and_pars),
                             &_err);
      lowp_log    = BEreq::clik_compute_loglike(byVal(lowl_clikid),
                             byVal(lowp_cl_and_pars),
                             &_err);


      result = highl_TT_log+lowp_log;

      std::cout << "Log likelihood is : " << result << std::endl;



    }


    // Simplest linear potential inflationary model with 7 params with Planck log-like.
    void function_1linearInf_rLCDMtensor_lowp_TT_loglike(double& result)
    {
      using namespace Pipes::function_1linearInf_rLCDMtensor_lowp_TT_loglike;

      // Initialization parameters controlling main characteristics.
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");

      // Control the output of analytic approximations for comparison.
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");

      // Parameters to control how the ICs are sampled.
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");

      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this

      // Parameters to control how the vparams are sampled.
      int param_sampling = runOptions->getValue<int> ("param_sampling");


      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");

      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");

      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");

      // Parameters to be passed to the potential
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

      //-------------------------------------------------
      //----------------- Sampling \lambda ------------------
      //-------------------------------------------------
      vparams[0] = *Param["lambda"];


      //-------------------------------------------------
      //--------------- Sampling N_pivot ----------------
      //-------------------------------------------------
      double N_pivot = *Param["N_pivot"];

      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");


      // Priors on the IC and N_pivot ranges
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");

      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");

      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");

      gambit_inflation_observables observs;

      std::cout << "num_inflaton = " << num_inflaton << std::endl;

      std::cout << "save_iso_N = " << save_iso_N << std::endl;

      std::cout << "N_pivot = " << N_pivot << std::endl;

      std::cout << "kmin = " << kmin << std::endl;

      std::cout << "dlnk = " << dlnk << std::endl;

      std::cout << "lambda = " << vparams[0] << std::endl;

      // The function below calls the multimodecode backend for a given choice of inflationary model,
      // which calculates the observables.
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
                         byVal(&dphi_init0[0]), // dphi correct at the multimodecode
                         byVal(&vparams[0]),
                         N_pivot,
                         k_pivot,
                         dlnk,
                         calc_full_pk,
                         steps,
                         kmin,
                         //                         kmax,
                         byVal(&phi0_priors_min[0]),
                         byVal(&phi0_priors_max[0]),
                         byVal(&dphi0_priors_min[0]),
                         byVal(&dphi0_priors_max[0]),
                         N_pivot_prior_min,
                         N_pivot_prior_max);


      std::cout << "observs.As " << observs.As << std::endl;
      std::cout << "observs.A_iso " << observs.A_iso << std::endl;
      std::cout << "observs.A_pnad " << observs.A_pnad << std::endl;
      std::cout << "observs.A_ent " << observs.A_ent << std::endl;
      std::cout << "observs.A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
      //std::cout << "observs.A_bundle " << observs.A_bundle << std::endl;
      std::cout << "observs.ns " << observs.ns << std::endl;
      std::cout << "observs.nt " << observs.nt << std::endl;
      std::cout << "observs.n_iso " << observs.n_iso << std::endl;
      std::cout << "observs.n_pnad " << observs.n_pnad << std::endl;
      std::cout << "observs.n_ent " << observs.n_ent << std::endl;
      std::cout << "observs.r " << observs.r << std::endl;
      std::cout << "observs.alpha_s " << observs.alpha_s << std::endl;
      std::cout << "observs.runofrun " << observs.runofrun << std::endl;
      std::cout << "observs.f_NL " << observs.f_NL << std::endl;
      std::cout << "observs.tau_NL " << observs.tau_NL << std::endl;

      //*/
      int l,l_max;
      l_max=3000;
      //  int nmax;
      int num_ct_max=7;

      double  highl_TT_cl_and_pars[2525];
      double  lowp_cl_and_pars[121];

      char *class_null;  // defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
       could add another array with P(k), or extract other results from
       the code - here we assume that we are interested in the C_l's
       only */

      double** cl;

      cl = return_LCDMtensor_1linearInf_cls(*Param["omega_b"],
                       *Param["omega_cdm"],
                       *Param["H0"],
                       //*Param["ln10A_s"],
                       observs.As,
                       //*Param["n_s"],
                       observs.ns,
                       *Param["tau_reio"],
                       //*Param["r_tensor"]
                       observs.r);
      //--------------------------------------------------------------------------
      //------high-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------

      highl_TT_cl_and_pars[0] = 0.0;
      highl_TT_cl_and_pars[1] = 0.0;

      for(int ii = 2; ii < 2509 ; ii++)highl_TT_cl_and_pars[ii] = cl[ii][0];

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      highl_TT_cl_and_pars[2509] = *Param["A_cib_217"];
      highl_TT_cl_and_pars[2510] = *Param["cib_index"];
      highl_TT_cl_and_pars[2511] = *Param["xi_sz_cib"];
      highl_TT_cl_and_pars[2512] = *Param["A_sz"];
      highl_TT_cl_and_pars[2513] = *Param["ps_A_100_100"];
      highl_TT_cl_and_pars[2514] = *Param["ps_A_143_143"];
      highl_TT_cl_and_pars[2515] = *Param["ps_A_143_217"];
      highl_TT_cl_and_pars[2516] = *Param["ps_A_217_217"];
      highl_TT_cl_and_pars[2517] = *Param["ksz_norm"];
      highl_TT_cl_and_pars[2518] = *Param["gal545_A_100"];
      highl_TT_cl_and_pars[2519] = *Param["gal545_A_143"];
      highl_TT_cl_and_pars[2520] = *Param["gal545_A_143_217"];
      highl_TT_cl_and_pars[2521] = *Param["gal545_A_217"];
      highl_TT_cl_and_pars[2522] = *Param["calib_100T"];
      highl_TT_cl_and_pars[2523] = *Param["calib_217T"];
      highl_TT_cl_and_pars[2524] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //-------low-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[0] = 0.0;
      lowp_cl_and_pars[1] = 0.0;
      for (l=2;l<=29;l++)  {
        lowp_cl_and_pars[l] = cl[l][0];
      }
      int k;
      lowp_cl_and_pars[30] = 0.0;
      lowp_cl_and_pars[31] = 0.0;
      for (l=32;l<=59;l++)  {
        k = l-30;
        lowp_cl_and_pars[l] = cl[k][1];
      }
      lowp_cl_and_pars[60] = 0.0;
      lowp_cl_and_pars[61] = 0.0;
      for (l=62;l<=89;l++)  {
        k = l-60;
        lowp_cl_and_pars[l] = cl[k][2];
      }
      lowp_cl_and_pars[90] = 0.0;
      lowp_cl_and_pars[91] = 0.0;
      for (l=92;l<=119;l++)  {
        k = l-90;
        lowp_cl_and_pars[l] = cl[k][3];
      }
      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[120] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_object* high_clikid;
      clik_object* lowl_clikid;
      clik_error *_err;
      double lowp_log;
      double highl_TT_log;

      lowl_clikid = BEreq::return_lowp_TT();
      high_clikid = BEreq::return_high_TT();

      _err = BEreq::clik_initialize_error();

      highl_TT_log = BEreq::clik_compute_loglike(byVal(high_clikid),
                             byVal(highl_TT_cl_and_pars),
                             &_err);
      lowp_log    = BEreq::clik_compute_loglike(byVal(lowl_clikid),
                             byVal(lowp_cl_and_pars),
                             &_err);


      result = highl_TT_log+lowp_log;

      std::cout << "Log likelihood is : " << result << std::endl;



    }

    // Simplest smash potential inflationary model with 7 params with Planck log-like.
    void function_smashInf_rLCDMtensor_lowp_TT_loglike(double& result)
    {
      using namespace Pipes::function_smashInf_rLCDMtensor_lowp_TT_loglike;

      // Initialization parameters controlling main characteristics.
      int num_inflaton = runOptions->getValue<int> ("num_inflaton");
      int potential_choice = runOptions->getValue<int> ("potential_choice");
      int slowroll_infl_end = runOptions->getValue<int> ("slowroll_infl_end");
      int instreheat = runOptions->getValue<int> ("instreheat");
      int vparam_rows = runOptions->getValue<int> ("vparam_rows");

      // Control the output of analytic approximations for comparison.
      int use_deltaN_SR = runOptions->getValue<int> ("use_deltaN_SR");
      int evaluate_modes = runOptions->getValue<int> ("evaluate_modes");
      int use_horiz_cross_approx = runOptions->getValue<int> ("use_horiz_cross_approx");
      int get_runningofrunning = runOptions->getValue<int> ("get_runningofrunning");

      // Parameters to control how the ICs are sampled.
      int ic_sampling = runOptions->getValue<int> ("ic_sampling");
      double energy_scale = runOptions->getValue<double> ("energy_scale");
      int numb_samples = runOptions->getValue<int> ("numb_samples");
      int save_iso_N = runOptions->getValue<int> ("save_iso_N");

      double N_iso_ref = runOptions->getValue<int> ("N_iso_ref"); //double check this

      // Parameters to control how the vparams are sampled.
      int param_sampling = runOptions->getValue<int> ("param_sampling");


      std::vector<double> vp_prior_min = runOptions->getValue<std::vector<double> >("vp_prior_min");
      std::vector<double> vp_prior_max = runOptions->getValue<std::vector<double> > ("vp_prior_max");

      int varying_N_pivot = runOptions->getValue<int> ("varying_N_pivot");
      int use_first_priorval = runOptions->getValue<int> ("use_first_priorval");

      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");
      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");

      // Parameters to be passed to the potential
      std::vector<double> vparams = runOptions->getValue<std::vector<double> >("vparams");

      //-------------------------------------------------
      //----------------- Sampling \lambda ------------------
      //-------------------------------------------------
      vparams[0] = *Param["log10_xi"];
      vparams[1] = *Param["log10_beta"];
      vparams[2] = *Param["log10_lambda"];
      // xi = 10**vparams(1,1)
      // beta = 10**vparams(2,1)
      // lambda = 10**vparams(3,1)

      std::cout << "log10[xi] = " << vparams[0] << std::endl;
      std::cout << "log10[beta] = " << vparams[1] << std::endl;
      std::cout << "log10[lambda] = " << vparams[2] << std::endl;

      int i;
      /* coefficients of P(x) =  */
      /*  -8 + b*\[Phi]**2 + b*\[Xi]*\[Phi]**4 + 6*\[Xi]**2*\[Phi]**4 */
      double smashp1[5] = { -8.0, 0, pow(10.0,vparams[1]), 0,
        pow(10.0,vparams[1])*pow(10.0,vparams[0])
        +6.0*pow(10.0,vparams[0])*pow(10.0,vparams[0])};

      double smashd1[8];

      gsl_poly_complex_workspace * wsmash
        = gsl_poly_complex_workspace_alloc (5);

      gsl_poly_complex_solve (smashp1, 5, wsmash, smashd1);

      gsl_poly_complex_workspace_free (wsmash);

      for (i = 0; i < 4; i++)
        {
        printf ("z%d = %+.18f %+.18f\n",
        i, smashd1[2*i], smashd1[2*i+1]);
      }

      double badway[4] = {smashd1[0],smashd1[2],smashd1[4],smashd1[6]};
      double tempbw = 0;

      for(int i=0;i<4;i++)
      {
        if(badway[i]>tempbw)
          tempbw=badway[i];
      }
      //gsl_sf_log((1.0+pow(10.0,vparams[0])*smashp1[4]*smashp1[4])/
      //       (1.0+pow(10.0,vparams[0])*)

      vparams[3] = tempbw;

      std::cout << "phi_end = " << tempbw << std::endl;

      //-------------------------------------------------
      //--------------- Sampling N_pivot ----------------
      //-------------------------------------------------
      double N_pivot = *Param["N_pivot"];


      ///////////// solving the equality to calculate the
      ////////////  slow roll field amplitute at N_pivot.
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

      printf ("using %s method\n",
          gsl_min_fminimizer_name (s));

      printf ("%5s [%9s, %9s] %9s %10s %9s\n",
          "iter", "lower", "upper", "min",
          "err", "err(est)");

      printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
          iter, a, b,
          m, m - m_expected, b - a);

      do
      {
        iter++;
        status = gsl_min_fminimizer_iterate (s);

        m = gsl_min_fminimizer_x_minimum (s);
        a = gsl_min_fminimizer_x_lower (s);
        b = gsl_min_fminimizer_x_upper (s);

        status = gsl_min_test_interval (a, b, 0.01, 0.0);

      if (status == GSL_SUCCESS) printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] "
          "%.7f %+.7f %.7f\n",
          iter, a, b,
          m, m - m_expected, b - a);
      }
      while (status == GSL_CONTINUE && iter < max_iter);

      gsl_min_fminimizer_free (s);

      //phi_init0[0] = a*1.01;


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
      r_self = r_SR(smash_eps,smash_eta);

      /*
      ///////////// have calculated the field value at N_pivot
      ///////////// predicted by slow-roll approximation.
      ///////////// Choosing a very close slightly higher
      ///////////// value will be satisfactory for initial conditions.

      double k_pivot = runOptions->getValue<double> ("k_pivot");
      double dlnk = runOptions->getValue<double> ("dlnk");


      // Priors on the IC and N_pivot ranges
      std::vector<double> phi0_priors_min = runOptions->getValue<std::vector<double> > ("phi0_priors_min");
      std::vector<double> phi0_priors_max = runOptions->getValue<std::vector<double> > ("phi0_priors_max");
      std::vector<double> dphi0_priors_min = runOptions->getValue<std::vector<double> > ("dphi0_priors_min");
      std::vector<double> dphi0_priors_max = runOptions->getValue<std::vector<double> > ("dphi0_priors_max");

      double N_pivot_prior_min = runOptions->getValue<double> ("N_pivot_prior_min");
      double N_pivot_prior_max = runOptions->getValue<double> ("N_pivot_prior_max");

      // For calculating the full power spectrum P(k). Samples in uniform increments in log(k).
      int calc_full_pk = runOptions->getValue<int> ("calc_full_pk");
      int steps = runOptions->getValue<int> ("steps");
      double kmin = runOptions->getValue<double> ("kmin");
      double kmax = runOptions->getValue<double> ("kmax");

      gambit_inflation_observables observs;

      //std::cout << "num_inflaton = " << num_inflaton << std::endl;

      //std::cout << "save_iso_N = " << save_iso_N << std::endl;

      std::cout << "N_pivot = " << N_pivot << std::endl;

      //std::cout << "kmin = " << kmin << std::endl;

      //std::cout << "dlnk = " << dlnk << std::endl;

      std::cout << "phi_end = " << vparams[3] << std::endl;
      std::cout << "phi_piv = " << a << std::endl;
      std::cout << "phi_init= " << phi_init0[0] << std::endl;

      // The function below calls the multimodecode backend for a given choice of inflationary model,
      // which calculates the observables.
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
                         byVal(&dphi_init0[0]), // dphi correct at the multimodecode
                         byVal(&vparams[0]),
                         N_pivot,
                         k_pivot,
                         dlnk,
                         calc_full_pk,
                         steps,
                         kmin,
                         //                         kmax,
                         byVal(&phi0_priors_min[0]),
                         byVal(&phi0_priors_max[0]),
                         byVal(&dphi0_priors_min[0]),
                         byVal(&dphi0_priors_max[0]),
                         N_pivot_prior_min,
                         N_pivot_prior_max);


      std::cout << "observs.As " << observs.As << std::endl;
      std::cout << "observs.A_iso " << observs.A_iso << std::endl;
      std::cout << "observs.A_pnad " << observs.A_pnad << std::endl;
      std::cout << "observs.A_ent " << observs.A_ent << std::endl;
      std::cout << "observs.A_cross_ad_iso " << observs.A_cross_ad_iso << std::endl;
      //std::cout << "observs.A_bundle " << observs.A_bundle << std::endl;
      std::cout << "observs.ns " << observs.ns << std::endl;
      std::cout << "observs.nt " << observs.nt << std::endl;
      std::cout << "observs.n_iso " << observs.n_iso << std::endl;
      std::cout << "observs.n_pnad " << observs.n_pnad << std::endl;
      std::cout << "observs.n_ent " << observs.n_ent << std::endl;
      std::cout << "observs.r " << observs.r << std::endl;
      std::cout << "observs.alpha_s " << observs.alpha_s << std::endl;
      std::cout << "observs.runofrun " << observs.runofrun << std::endl;
      std::cout << "observs.f_NL " << observs.f_NL << std::endl;
      std::cout << "observs.tau_NL " << observs.tau_NL << std::endl;


      */
      //*/
      int l,l_max;
      l_max=3000;
      //  int nmax;
      int num_ct_max=7;

      double  highl_TT_cl_and_pars[2525];
      double  lowp_cl_and_pars[121];

      char *class_null;  // defining an empty string for parser initialization.

      /* allocate the array where calculated Cl's will be written (we
       could add another array with P(k), or extract other results from
       the code - here we assume that we are interested in the C_l's
       only */

      double** cl;

      cl = return_LCDMtensor_smashInf_cls(*Param["omega_b"],
                          *Param["omega_cdm"],
                          *Param["H0"],
                          //*Param["ln10A_s"],
                          //observs.As,
                          As_self,
                          //*Param["n_s"],
                          //observs.ns,
                          ns_self,
                          *Param["tau_reio"],
                          //*Param["r_tensor"]
                          //observs.r
                          r_self);
      //--------------------------------------------------------------------------
      //------high-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------

      highl_TT_cl_and_pars[0] = 0.0;
      highl_TT_cl_and_pars[1] = 0.0;

      for(int ii = 2; ii < 2509 ; ii++)highl_TT_cl_and_pars[ii] = cl[ii][0];

      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      highl_TT_cl_and_pars[2509] = *Param["A_cib_217"];
      highl_TT_cl_and_pars[2510] = *Param["cib_index"];
      highl_TT_cl_and_pars[2511] = *Param["xi_sz_cib"];
      highl_TT_cl_and_pars[2512] = *Param["A_sz"];
      highl_TT_cl_and_pars[2513] = *Param["ps_A_100_100"];
      highl_TT_cl_and_pars[2514] = *Param["ps_A_143_143"];
      highl_TT_cl_and_pars[2515] = *Param["ps_A_143_217"];
      highl_TT_cl_and_pars[2516] = *Param["ps_A_217_217"];
      highl_TT_cl_and_pars[2517] = *Param["ksz_norm"];
      highl_TT_cl_and_pars[2518] = *Param["gal545_A_100"];
      highl_TT_cl_and_pars[2519] = *Param["gal545_A_143"];
      highl_TT_cl_and_pars[2520] = *Param["gal545_A_143_217"];
      highl_TT_cl_and_pars[2521] = *Param["gal545_A_217"];
      highl_TT_cl_and_pars[2522] = *Param["calib_100T"];
      highl_TT_cl_and_pars[2523] = *Param["calib_217T"];
      highl_TT_cl_and_pars[2524] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //-------low-l likelihood calculation making of Cls-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[0] = 0.0;
      lowp_cl_and_pars[1] = 0.0;
      for (l=2;l<=29;l++)  {
        lowp_cl_and_pars[l] = cl[l][0];
      }
      int k;
      lowp_cl_and_pars[30] = 0.0;
      lowp_cl_and_pars[31] = 0.0;
      for (l=32;l<=59;l++)  {
        k = l-30;
        lowp_cl_and_pars[l] = cl[k][1];
      }
      lowp_cl_and_pars[60] = 0.0;
      lowp_cl_and_pars[61] = 0.0;
      for (l=62;l<=89;l++)  {
        k = l-60;
        lowp_cl_and_pars[l] = cl[k][2];
      }
      lowp_cl_and_pars[90] = 0.0;
      lowp_cl_and_pars[91] = 0.0;
      for (l=92;l<=119;l++)  {
        k = l-90;
        lowp_cl_and_pars[l] = cl[k][3];
      }
      //--------------------------------------------------------------------------
      //------addition of nuisance parameters to Cl array-------------------------
      //--------------------------------------------------------------------------
      lowp_cl_and_pars[120] = *Param["A_planck"];

      //--------------------------------------------------------------------------
      //------calculation of the planck loglikelihood-----------------------------
      //--------------------------------------------------------------------------
      clik_object* high_clikid;
      clik_object* lowl_clikid;
      clik_error *_err;
      double lowp_log;
      double highl_TT_log;

      lowl_clikid = BEreq::return_lowp_TT();
      high_clikid = BEreq::return_high_TT();

      _err = BEreq::clik_initialize_error();

      highl_TT_log = BEreq::clik_compute_loglike(byVal(high_clikid),
                             byVal(highl_TT_cl_and_pars),
                             &_err);
      lowp_log    = BEreq::clik_compute_loglike(byVal(lowl_clikid),
                             byVal(lowp_cl_and_pars),
                             &_err);


      result = highl_TT_log+lowp_log;

      std::cout << "Log likelihood is : " << result << std::endl;



    }

  }

}
