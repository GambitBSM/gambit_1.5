//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions of module CosmoBit.
///
///  Put your functions in files like this
///  if you wish to add observables or likelihoods
///  to this module.
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
///  \date 2018 Jan,Feb
///
///  *********************************************
#include <string>
#include <iostream>
#include <cmath>
#include <functional>
#include <omp.h>
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

    void dummy_func(double& result)
    {
      using namespace Pipes::dummy_func;

      std::vector<double> vec;
      BEreq::DA_awesomeness(vec);
      double sum = 0.0;
      for (std::vector<double>::iterator it = vec.begin(); it != vec.end(); ++it)
      {
        sum += *it;
      }
      result = sum;
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

      BEreq::class_parser_initialize(&cosmo.fc,int(len_of_input),"",cosmo.class_errmsg);

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

      BEreq::class_parser_initialize(&cosmo.fc,10,"",cosmo.class_errmsg);

      strcpy(cosmo.fc.name[0],"output");
      strcpy(cosmo.fc.value[0],"tCl pCl lCl");
      strcpy(cosmo.fc.name[8],"l_max_scalars");
      sprintf(cosmo.fc.value[8],"%d",l_max);
      strcpy(cosmo.fc.name[9],"lensing");
      strcpy(cosmo.fc.value[9],"yes");

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

      BEreq::class_parser_initialize(&cosmo.fc,11,"",cosmo.class_errmsg);

      strcpy(cosmo.fc.name[0],"output");
      strcpy(cosmo.fc.value[0],"tCl pCl lCl");
      strcpy(cosmo.fc.name[8],"l_max_scalars");
      sprintf(cosmo.fc.value[8],"%d",l_max);
      strcpy(cosmo.fc.name[9],"modes");
      strcpy(cosmo.fc.value[9],"s,t");
      strcpy(cosmo.fc.name[10],"lensing");
      strcpy(cosmo.fc.value[10],"yes");

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

      //std::cout << "vp_prior_max : " << vp_prior_max << std::endl;
      // Parameters to be passed to the potential
      std::vector<double> phi_init0 = runOptions->getValue<std::vector<double> >("phi_init0");

      std::vector<double> dphi_init0 = runOptions->getValue<std::vector<double> >("dphi_init0");

      //Farray<double,1,10> dphi_init0 = runOptions->getValue<Farray<double,1,10>> ("dphi_init0");
      // ideally - this should be a double vector. This is hard.

      //std::vector<vector<double> > vparams =  runOptions->getValue<std::vector<vector<double> > >("vparams");
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

      //for (std::vector<double>::const_iterator i = vparams.begin(); i != vparams.end(); ++i) std::cout << *i << ' ';

      //std::cout << "Calling the function: multimodecode_gambit_driver" << std::end;


      gambit_inflation_observables observs;

      std::cout << "num_inflaton = " << num_inflaton << std::endl;

      std::cout << "save_iso_N = " << save_iso_N << std::endl;


/*      BEreq::multimodecode_gambit_driver_test(&observs,
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
                          num_inflaton,
                          save_iso_N,
                          N_iso_ref,
                          param_sampling,
                          varying_N_pivot,
                          use_first_priorval,
                          N_pivot,
                          k_pivot,
                          dlnk,
                          calc_full_pk,
                          steps,
                          kmin,
                          kmax,
                          N_pivot_prior_min); // last piece wrong
*/
      std::cout << "N_pivot = " << N_pivot << std::endl;

      std::cout << "kmin = " << kmin << std::endl;

      std::cout << "dlnk = " << dlnk << std::endl;


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

  }

}
