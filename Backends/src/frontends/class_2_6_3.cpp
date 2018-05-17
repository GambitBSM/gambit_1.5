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
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/class_2_6_3.hpp"

// Convenience functions (definitions)
BE_NAMESPACE
{
  CosmoBit::Class_container cosmo;
  static int last_success=0;
  static bool check_for_unused = true;

  void class_2_6_3_free()
  {
    if (last_success > 8) class_lensing_free(&cosmo.le);
    if (last_success > 7) class_spectra_free(&cosmo.sp);
    if (last_success > 6) class_transfer_free(&cosmo.tr);
    if (last_success > 5) class_nonlinear_free(&cosmo.nl);
    if (last_success > 4) class_primordial_free(&cosmo.pm);
    if (last_success > 3) class_perturb_free(&cosmo.pt);
    if (last_success > 2) class_thermodynamics_free(&cosmo.th);
    if (last_success > 1) class_background_free(&cosmo.ba);
    last_success = 0;
  }

  void class_2_6_3_run()
  {
    //std::cout << "Last seen alive in: class_2_6_3_run" << std::endl;

    char error_printout[1024];

    if (class_input_initialize(&cosmo.fc,&cosmo.pr,&cosmo.ba,&cosmo.th,&cosmo.pt,&cosmo.tr,&cosmo.pm,&cosmo.sp,&cosmo.nl,&cosmo.le,&cosmo.op,cosmo.class_errmsg) == _FAILURE_)
    {
      sprintf(error_printout,"Error in class_input_initialize\n=>%s\n",cosmo.class_errmsg);
      invalid_point().raise(error_printout);
    }
    last_success++;
    /** - Check for unused parameters */

    /* (If this function is run for the first time, we check if all parameters were read.
       This is done by a loop through all entries of cosmo.fc. When the input
       was successfully recognized, cosmo.fc.read[i] is set to be true.
       If at least one parameter is not used, we throw an error, to prevent the user from
       calculating from things.) */
    if (check_for_unused)
    {
      std::vector< std::string > unused;
      for (int i=0; i<cosmo.fc.size; i++)
      {
        if (cosmo.fc.read[i] != _TRUE_)
        {
          unused.push_back(std::string(cosmo.fc.name[i]));
        }
      }
      if (!(unused.empty()))
      {
        std::string error_out;
        error_out = "Some of your parameters are not read by class.\n\nThese are:\n\n";
        for (auto par = unused.begin(); par != unused.end(); par++)
        {
          error_out +=  "--> ";
          error_out +=  *par;
          error_out +=  " <--\n";
        }
        error_out += "\nPlease fix this. Most probably these unknown parameters are given in the Rules section.";
        backend_error().raise(LOCAL_INFO,error_out);
      }
      check_for_unused = false;
    }

    /** - Run CLASS (Go through all class_xxx_initialize) */

    if (class_background_initialize(&cosmo.pr,&cosmo.ba) == _FAILURE_)
    {
      sprintf(error_printout,"Error in class_background_initialize\n=>%s\n",cosmo.ba.error_message);
      invalid_point().raise(error_printout);
    }
    last_success++;
    if (class_thermodynamics_initialize(&cosmo.pr,&cosmo.ba,&cosmo.th) == _FAILURE_)
    {
      sprintf(error_printout,"Error in class_thermodynamics_initialize\n=>%s\n",cosmo.th.error_message);
      invalid_point().raise(error_printout);
    }
    last_success++;
    if (class_perturb_initialize(&cosmo.pr,&cosmo.ba,&cosmo.th,&cosmo.pt) == _FAILURE_)
    {
      sprintf(error_printout,"Error in class_perturb_initialize\n=>%s\n",cosmo.pt.error_message);
      invalid_point().raise(error_printout);
    }
    last_success++;

    /* If full pk is asked on the set parameter function; this sets ppm->primordial_spec_type
       to be equal to 'gambit_Pk'. Using this flag, we externally fill in the primordial power spectrum.
       (note alterative -without patching class at all- is to use text files and 'cat' command which we
       want to aviod. The existing CLASS on the branch is already patched to have necessary flags and
       functions to allow for hacking and external setting the primordial power spectrum. */
    if(cosmo.pm.primordial_spec_type == 1) // is 1 when Gambit_pk is asked for.
    {
      /* The bin number of the primordial power spectrum.
       This is not to be set here (and also should be given by the user).
       Will change it ASAP.

       Furthermore, the arrays:
       cosmo.k_ar, cosmo.Pk_S, and cosmo.Pk_T are to be filled
       in the early step of setting the parameters. */

      cosmo.pm.lnk_size = 100;

      /** - Make room */
      cosmo.pm.lnk = (double *)malloc(cosmo.pm.lnk_size*sizeof(double));

      std::cout << "we pass pm.lnk \n" << std::endl;

      cosmo.pm.lnpk = (double **)malloc(cosmo.pt.md_size*sizeof(double));
      cosmo.pm.ddlnpk = (double **)malloc(cosmo.pt.md_size*sizeof(double));
      cosmo.pm.ic_size = (int *)malloc(cosmo.pt.md_size*sizeof(int));
      cosmo.pm.ic_ic_size = (int *)malloc(cosmo.pt.md_size*sizeof(int));
      cosmo.pm.is_non_zero = (short **)malloc(cosmo.pt.md_size*sizeof(short));

      int index_md;
      for (index_md = 0; index_md < cosmo.pt.md_size; index_md++)
      {
        cosmo.pm.ic_size[index_md] = cosmo.pt.ic_size[index_md];
        cosmo.pm.ic_ic_size[index_md] = (cosmo.pm.ic_size[index_md]*(cosmo.pm.ic_size[index_md]+1))/2;

        std::cout << "pm.ic_size["<<index_md<<"] =  " << cosmo.pt.ic_size[index_md] << std::endl;
        std::cout << "pm.ic_ic_size["<<index_md<<"] = " << cosmo.pm.ic_ic_size[index_md]<<std::endl;

        cosmo.pm.lnpk[index_md] = (double *)malloc(cosmo.pm.lnk_size*cosmo.pm.ic_ic_size[index_md]*sizeof(double));
        cosmo.pm.ddlnpk[index_md] = (double *)malloc(cosmo.pm.lnk_size*cosmo.pm.ic_ic_size[index_md]*sizeof(double));
        cosmo.pm.is_non_zero[index_md] = (short *)malloc(cosmo.pm.lnk_size*cosmo.pm.ic_ic_size[index_md]*sizeof(short));
      }

      /** - Store values */
      for (int index_k=0; index_k<cosmo.pm.lnk_size; index_k++)
      {
        std::cout << "k_array["<<index_k<<"]="<< cosmo.k_ar.at(index_k) <<std::endl;
        std::cout << "pks_array["<<index_k<<"]="<< cosmo.Pk_S.at(index_k)<<std::endl;
        std::cout << "pkt_array["<<index_k<<"]="<< cosmo.Pk_T.at(index_k)<<std::endl;

        cosmo.pm.lnk[index_k] = std::log( cosmo.k_ar.at(index_k) );
        cosmo.pm.lnpk[cosmo.pt.index_md_scalars][index_k] = std::log( cosmo.Pk_S.at(index_k) );
        if (cosmo.pt.has_tensors == _TRUE_)
          cosmo.pm.lnpk[cosmo.pt.index_md_tensors][index_k] = std::log( cosmo.Pk_T.at(index_k) );

        std::cout << "pm.lnk["<<index_k<<"]="<<cosmo.pm.lnk[index_k]<<std::endl;
        std::cout << "pm.lnpk["<<cosmo.pt.index_md_scalars<<"]["<<index_k<<"]="<<cosmo.pm.lnpk[cosmo.pt.index_md_scalars][index_k]<<std::endl;
        std::cout << "pm.lnpk["<<cosmo.pt.index_md_tensors<<"]["<<index_k<<"]="<<cosmo.pm.lnpk[cosmo.pt.index_md_tensors][index_k]<<std::endl;
      }

      /** - Tell CLASS that there are scalar (and tensor) modes */
      cosmo.pm.is_non_zero[cosmo.pt.index_md_scalars][cosmo.pt.index_ic_ad] = _TRUE_;
      if (cosmo.pt.has_tensors == _TRUE_)
      cosmo.pm.is_non_zero[cosmo.pt.index_md_tensors][cosmo.pt.index_ic_ten] = _TRUE_;
    }

    if (class_primordial_initialize(&cosmo.pr,&cosmo.pt,&cosmo.pm) == _FAILURE_)
    {
      sprintf(error_printout,"Error in class_primordial_initialize\n=>%s\n",cosmo.pm.error_message);
      invalid_point().raise(error_printout);
    }
    last_success++;
    if (class_nonlinear_initialize(&cosmo.pr,&cosmo.ba,&cosmo.th,&cosmo.pt,&cosmo.pm,&cosmo.nl) == _FAILURE_)
    {
      sprintf(error_printout,"Error in class_nonlinear_initialize\n=>%s\n",cosmo.nl.error_message);
      invalid_point().raise(error_printout);
    }
    last_success++;
    if (class_transfer_initialize(&cosmo.pr,&cosmo.ba,&cosmo.th,&cosmo.pt,&cosmo.nl,&cosmo.tr) == _FAILURE_)
    {
      sprintf(error_printout,"Error in class_transfer_initialize\n=>%s\n",cosmo.tr.error_message);
      invalid_point().raise(error_printout);
    }
    last_success++;
    if (class_spectra_initialize(&cosmo.pr,&cosmo.ba,&cosmo.pt,&cosmo.pm,&cosmo.nl,&cosmo.tr,&cosmo.sp) == _FAILURE_)
    {
      sprintf(error_printout,"Error in class_spectra_initialize\n=>%s\n",cosmo.sp.error_message);
      invalid_point().raise(error_printout);
    }
    last_success++;
    if (class_lensing_initialize(&cosmo.pr,&cosmo.pt,&cosmo.sp,&cosmo.nl,&cosmo.le) == _FAILURE_)
    {
      sprintf(error_printout,"Error in class_lensing_initialize\n=>%s\n",cosmo.le.error_message);
      invalid_point().raise(error_printout);
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
    class_parser_initialize(&cosmo.fc,byVal(len_of_input),(char *)"",cosmo.class_errmsg);
    for(auto iter=in_map.begin(); iter != in_map.end(); iter++)
    {
      auto pos = std::distance(in_map.begin(), iter );
      std::strncpy(cosmo.fc.name[pos], iter->first.c_str(), sizeof(Class::FileArg));
      std::strncpy(cosmo.fc.value[pos], iter->second.c_str(), sizeof(Class::FileArg));
    }
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
  scan_level = false;

  cosmo = *Dep::class_set_parameter;
  CosmoBit::ClassInput in_dict = cosmo.input;
  class_2_6_3_set_parameter(in_dict);
  in_dict.clear();
  class_2_6_3_run();
}
END_BE_INI_FUNCTION
