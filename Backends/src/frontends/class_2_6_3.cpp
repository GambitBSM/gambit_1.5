//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for the class 2.6.3 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2018 Apr, May
///  \date 2019 Feb
///
///  \author Selim Hotinli
///  \date 2018 May-June
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 Oct
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/class_2_6_3.hpp"


// Convenience functions (definitions)
BE_NAMESPACE
{
  CosmoBit::Class_container cosmo;

  Class::file_content fc;
  Class::precision pr;
  Class::background ba;
  Class::thermo th;
  Class::perturbs pt;
  Class::transfers tr;
  Class::primordial pm;
  Class::spectra sp;
  Class::nonlinear nl;
  Class::lensing le;
  Class::output op;
  Class::ErrorMsg errmsg;

  static int last_success = 0;
  static bool check_for_unused = true;
  static str maxlevel;

  void class_2_6_3_free()
  {
    if (last_success > 8) class_lensing_free(&le);
    if (last_success > 7) class_spectra_free(&sp);
    if (last_success > 6) class_transfer_free(&tr);
    if (last_success > 5) class_nonlinear_free(&nl);
    if (last_success > 4) class_primordial_free(&pm);
    if (last_success > 3) class_perturb_free(&pt);
    if (last_success > 2) class_thermodynamics_free(&th);
    if (last_success > 1) class_background_free(&ba);
    last_success = 0;
  }

  void class_2_6_3_run()
  {
    //std::cout << "Last seen alive in: class_2_6_3_run" << std::endl;

    std::string error_printout;

    if (class_input_initialize(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_)
    {
      error_printout = "Error in class_input_initialize:\n\n ";
      error_printout += errmsg;
      invalid_point().raise(error_printout.c_str());
    }
    last_success++;
    /** - Check for unused parameters */

    /* (If this function is run for the first time, we check if all parameters were read.
       This is done by a loop through all entries of fc. When the input
       was successfully recognized, fc.read[i] is set to be true.
       If at least one parameter is not used, we throw an error, to prevent the user from
       calculating from things.) */

    if (check_for_unused)
    {
      std::vector< std::string > unused;
      for (int i=0; i<fc.size; i++)
      {
	if (fc.read[i] != _TRUE_)
	{
	  unused.push_back(std::string(fc.name[i]));
	}
      }
      if (!(unused.empty()))
      {
	error_printout = "Some of your parameters are not read by class.\n\nThese are:\n\n";
	for (auto par = unused.begin(); par != unused.end(); par++)
	{
	  error_printout +=  "--> ";
	  error_printout +=  *par;
	  error_printout +=  " <--\n";
	}
	error_printout += "\nPlease fix this. Most probably these unknown parameters are given in the Rules section.";
	backend_error().raise(LOCAL_INFO,error_printout.c_str());
      }
      check_for_unused = false;
    }

    /** - Run CLASS (Go through all class_xxx_initialize) */

    if (class_background_initialize(&pr,&ba) == _FAILURE_)
    {
      error_printout = "Error in class_background_initialize:\n\n ";
      error_printout += ba.error_message;
      invalid_point().raise(error_printout.c_str());
    }
    last_success++;
    if (maxlevel.compare("background") == 0) return;
    if (class_thermodynamics_initialize(&pr,&ba,&th) == _FAILURE_)
    {
      error_printout = "Error in class_thermodynamics_initialize:\n\n ";
      error_printout += th.error_message;
      invalid_point().raise(error_printout.c_str());
    }
    last_success++;
    if (maxlevel.compare("thermodynamics") == 0) return;
    if (class_perturb_initialize(&pr,&ba,&th,&pt) == _FAILURE_)
    {
      error_printout = "Error in class_perturb_initialize:\n\n ";
      error_printout += pt.error_message;
      invalid_point().raise(error_printout.c_str());
    }
    last_success++;
    if (maxlevel.compare("perturb") == 0) return;
    /* If full pk is asked on the set parameter function; this sets ppm->primordial_spec_type
       to be equal to 'gambit_Pk'. Using this flag, we externally fill in the primordial power spectrum.
       (note alterative -without patching class at all- is to use text files and 'cat' command which we
       want to aviod. The existing CLASS on the branch is already patched to have necessary flags and
       functions to allow for hacking and external setting the primordial power spectrum. */
    if(pm.primordial_spec_type == 6) // is 6 when Gambit_pk is asked for.
    {
      /* The bin number of the primordial power spectrum.
       This is not to be set here (and also should be given by the user).
       Will change it ASAP.

       Furthermore, the arrays:
       cosmo.k_ar, cosmo.Pk_S, and cosmo.Pk_T are to be filled
       in the early step of setting the parameters. */

      pm.lnk_size = 30;
      pm.md_size = pt.md_size;

      /** - Make room */
      pm.lnk = (double *)malloc(pm.lnk_size*sizeof(double));
      // std::cout << "DEBUG: we pass pm.lnk \n" << std::endl;

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

	pm.lnpk[index_md] = (double *)malloc(pm.lnk_size*pm.ic_ic_size[index_md]*sizeof(double));
	pm.ddlnpk[index_md] = (double *)malloc(pm.lnk_size*pm.ic_ic_size[index_md]*sizeof(double));
	pm.is_non_zero[index_md] = (short *)malloc(pm.lnk_size*pm.ic_ic_size[index_md]*sizeof(short));
      }

      /** - Store values */
      for (int index_k=0; index_k<pm.lnk_size; index_k++)
      {
	std::cout <<  cosmo.k_ar.at(index_k) << "," << cosmo.Pk_S.at(index_k)<<std::endl;

	 pm.lnk[index_k] = std::log( cosmo.k_ar.at(index_k) );
	 pm.lnpk[pt.index_md_scalars][index_k] = std::log( cosmo.Pk_S.at(index_k) );
	 if (pt.has_tensors == _TRUE_)
	   pm.lnpk[pt.index_md_tensors][index_k] = std::log( cosmo.Pk_T.at(index_k) );
      }

      /** - Tell CLASS that there are scalar (and tensor) modes */
      pm.is_non_zero[pt.index_md_scalars][pt.index_ic_ad] = _TRUE_;
      if (pt.has_tensors == _TRUE_)
	pm.is_non_zero[pt.index_md_tensors][pt.index_ic_ten] = _TRUE_;
    }

    if (class_primordial_initialize(&pr,&pt,&pm) == _FAILURE_)
    {
      error_printout = "Error in class_primordial_initialize:\n\n ";
      error_printout += pm.error_message;
      invalid_point().raise(error_printout.c_str());
    }
    last_success++;
    if (maxlevel.compare("primordial") == 0) return;
    if (class_nonlinear_initialize(&pr,&ba,&th,&pt,&pm,&nl) == _FAILURE_)
    {
      error_printout = "Error in class_nonlinear_initialize:\n\n ";
      error_printout += nl.error_message;
      invalid_point().raise(error_printout.c_str());
    }
    last_success++;
    if (maxlevel.compare("nonlinear") == 0) return;
    if (class_transfer_initialize(&pr,&ba,&th,&pt,&nl,&tr) == _FAILURE_)
    {
      error_printout = "Error in class_transfer_initialize:\n\n ";
      error_printout += tr.error_message;
      invalid_point().raise(error_printout.c_str());
    }
    last_success++;
    if (maxlevel.compare("transfer") == 0) return;
    if (class_spectra_initialize(&pr,&ba,&pt,&pm,&nl,&tr,&sp) == _FAILURE_)
    {
      error_printout = "Error in class_spectra_initialize:\n\n ";
      error_printout += sp.error_message;
      invalid_point().raise(error_printout.c_str());
    }
    last_success++;
    if (maxlevel.compare("spectra") == 0) return;
    if (class_lensing_initialize(&pr,&pt,&sp,&nl,&le) == _FAILURE_)
    {
      error_printout = "Error in class_lensing_initialize:\n\n ";
      error_printout += le.error_message;
      invalid_point().raise(error_printout.c_str());
    }
    last_success++;
  }

  CosmoBit::Class_container get_class_2_6_3()
  {
    return cosmo;
  }

  void class_2_6_3_set_parameter(CosmoBit::ClassInput in_dict)
  {
    std::map<std::string,std::string> in_map = in_dict.get_map();
    int len_of_input = in_map.size();
    class_parser_initialize(&fc,byVal(len_of_input),(char *)"",errmsg);
    for(auto iter=in_map.begin(); iter != in_map.end(); iter++)
    {
      auto pos = std::distance(in_map.begin(), iter );
      std::strncpy(fc.name[pos], iter->first.c_str(), sizeof(Class::FileArg));
      std::strncpy(fc.value[pos], iter->second.c_str(), sizeof(Class::FileArg));
    }
  }

  double class_get_Da(double z)
  {
    double tau;
    int index;
    double *pvecback;
    //transform redshift in conformal time
    background_tau_of_z(&ba,z,&tau);

    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    //call to fill pvecback
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);

    //double H_z=pvecback[ba.index_bg_H];
    double Da=pvecback[ba.index_bg_ang_distance];
    free(pvecback);
    return Da;
  }

  double class_get_Dl(double z)
  {
    double tau=0;
    int index;
    double *pvecback;

    background_tau_of_z(&ba,z,&tau);

    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    //call to fill pvecback
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);

    double Dl = pvecback[ba.index_bg_lum_distance];
    free(pvecback);
    return Dl;
  }

  double class_get_scale_independent_growth_factor(double z)
  {
    double tau;
    int index;
    double *pvecback;
    //transform redshift in conformal time
    background_tau_of_z(&ba,z,&tau);

    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    //call to fill pvecback
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);

    //double H_z=pvecback[ba.index_bg_H];
    double D =pvecback[ba.index_bg_D];
    free(pvecback);
    return D;
  }

  double class_get_scale_independent_growth_factor_f(double z)
  {
    double tau;
    int index;
    double *pvecback;
    //transform redshift in conformal time
    background_tau_of_z(&ba,z,&tau);

    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    //call to fill pvecback
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);

    //double H_z=pvecback[ba.index_bg_H];
    double f =pvecback[ba.index_bg_f];
    free(pvecback);
    return f;
  }

  double class_get_Hz(double z)
  {
    double tau;
    int index;
    double *pvecback;
    //transform redshift in conformal time
    background_tau_of_z(&ba,z,&tau);

    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    //call to fill pvecback
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);

    double H_z=pvecback[ba.index_bg_H];
    free(pvecback);
    return(H_z);
  }

  double class_get_Omega_m()
  {
    return ba.Omega0_cdm + ba.Omega0_dcdmdr + ba.Omega0_b + ba.Omega0_ncdm_tot;
  }

  double class_get_rs()
  {
    return th.rs_d;
  }
  
  std::vector<std::vector<double>> class_get_z_of_r(std::vector<double> z_array)
  {
    double tau;
    int index;
    double *pvecback;
    std::vector<double> r;
    std::vector<double> dzdr;
    std::vector<std::vector<double>> result;

    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    for(int ii = 0; ii < z_array.size(); ii++)
    {
      //transform redshift in conformal time
      background_tau_of_z(&ba,z_array[ii],&tau);

      //call to fill pvecback
      background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);

      r.push_back(pvecback[ba.index_bg_conf_distance]);
      dzdr.push_back(pvecback[ba.index_bg_H]);
    }

    result.push_back(r);
    result.push_back(dzdr);
    free(pvecback);
    
    return result;
  }

  double class_get_sigma8(double z)
  {
    double tau;
    int index;
    double *pvecback;
    double sigma8 = 0.;
    //transform redshift in conformal time
    background_tau_of_z(&ba,z,&tau);

    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    //call to fill pvecback
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);
    spectra_sigma(&ba,&pm,&sp,8./ba.h,z,&sigma8);
    free(pvecback);
    return sigma8;

  }

  std::vector<double> class_get_cl(std::string spectype)
  {
    // Maximal value of l.
    int l_max = cosmo.lmax;

    std::vector<double> cl(l_max+1,0.);

    // Number of Cl-spectra (columns of the Cl-table).
    // The order of the spectra is [TT, EE, TE, BB, PhiPhi, TPhi, EPhi]
    int num_ct_max=7;

    // Define an array which takes the values of Cl of the different
    // spectra at a given value of l
    double* cltemp = new double[num_ct_max];

    // Loop through all l from 0 to l_max (including) and ask for the cl-spectra.
    for (int l=0; l <= l_max; l++)
    {
      if (l < 2)
      {
	// The entries for l=0 and l=1 are zero per defintion
	cl.at(l) = 0.;
      }
      else
      {
	// For l >= 2 ask for the cl-spectra.
	if (class_output_total_cl_at_l(&sp,&le,&op,byVal(l),byVal(cltemp)) == _SUCCESS_)
	{
	  if (spectype.compare("tt") == 0)
	    cl.at(l) = cltemp[sp.index_ct_tt]*pow(ba.T_cmb*1.e6,2);
	  else if (spectype.compare("te") == 0)
	    cl.at(l) = cltemp[sp.index_ct_te]*pow(ba.T_cmb*1.e6,2);
	  else if (spectype.compare("ee") == 0)
	    cl.at(l) = cltemp[sp.index_ct_ee]*pow(ba.T_cmb*1.e6,2);
	  else if (spectype.compare("bb") == 0)
	    cl.at(l) = cltemp[sp.index_ct_bb]*pow(ba.T_cmb*1.e6,2);
	  else if (spectype.compare("pp") == 0)
	    cl.at(l) = cltemp[sp.index_ct_pp];
	  else
	  {
	    std::string error_out;
	    error_out = "The \'spectype\' -> ";
	    error_out +=  spectype;
	    error_out += " <- is not recongnized and cannot be calculated by this Backend.";
	    backend_error().raise(LOCAL_INFO,error_out);
	  }
	}
	else
	{
	  // Failsafe for unexpected behaviour of "class_outpout_at_cl"
	  cl.at(l) = 0.;
	}
      }
    }

    // We do not need "cl" anymore
    delete cltemp;

    return cl;
  }



}
END_BE_NAMESPACE

// Initialisation function (definition)
BE_INI_FUNCTION
{
  static bool scan_level = true;
  if (!scan_level)
  {
    class_2_6_3_free();
  }
  else
  {
    if (*InUse::class_get_Da) maxlevel = "background" ;
    if (*InUse::class_get_Dl) maxlevel = "background" ;
    if (*InUse::class_get_Hz) maxlevel= "background";
    if (*InUse::class_get_Omega_m) maxlevel = "background" ;
    if (*InUse::class_get_rs) maxlevel = "thermodynamics";
    if (*InUse::class_get_sigma8) maxlevel = "lensing";
    if (*InUse::class_get_cl) maxlevel = "lensing";
  }
  scan_level = false;
  cosmo = *Dep::class_set_parameter;
  CosmoBit::ClassInput in_dict = cosmo.input;
  class_2_6_3_set_parameter(in_dict);
  in_dict.clear();
  class_2_6_3_run();
}
END_BE_INI_FUNCTION
