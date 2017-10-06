//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend header for the plc backend.
///
///  Compile-time registration of available
///  functions and variables from this backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Selim C. Hotinli
///          (selim.hotinli14@pimperial.ac.uk)
///  \date 2016 Sep
///
///  *********************************************

#define BACKENDNAME plc
#define VERSION 2.0
#define SAFE_VERSION 2_0

// Load it
LOAD_LIBRARY

BE_FUNCTION(clik_init, clik_object*, (char*,clik_error**),"clik_init","clik_initialize")
BE_FUNCTION(initError, clik_error* , () ,"initError","clik_initialize_error")
BE_FUNCTION(clik_compute, double, (clik_object*,double*,clik_error**), "clik_compute","clik_compute_loglike")
BE_FUNCTION(clik_get_has_cl, void, (clik_object*,int,clik_error**), "clik_get_has_cl","clik_get_list_cls")
BE_FUNCTION(clik_get_extra_parameter_names, int ,(clik_object*,parname**,clik_error**),"clik_get_extra_parameter_names","clik_get_extra_parameter_names")
BE_FUNCTION(clik_get_lmax, void , (clik_object*,int,clik_error**), "clik_get_lmax", "clik_get_lmax")
BE_FUNCTION(clik_cleanup, void , (clik_object**) , "clik_cleanup", "clik_cleanup")
BE_FUNCTION(clik_get_version, char*, (clik_object*,clik_error**), "clik_get_version" , "clik_get_version")
BE_CONV_FUNCTION(return_clikid_plik_dx11dr2_HM_v18_TT, clik_object*, (), "return_high_TT", (LCDM,LCDMtensor))
BE_CONV_FUNCTION(return_clikid_lowl_SMW_70_dx11d_2014, clik_object*, (), "return_lowp_TT", (LCDM,LCDMtensor))
BE_CONV_FUNCTION(data_initialize,void,(),"data_initialize",(LCDM,LCDMtensor))
BE_CONV_FUNCTION(data_cleanup,void,(),"data_cleanup",(LCDM,LCDMtensor))

BE_INI_FUNCTION
{
	static bool scan_level = true;
	if (scan_level)
	{
		data_initialize();
	}
	scan_level = false;
}
END_BE_INI_FUNCTION


BE_NAMESPACE
{
	
	clik_object * clikid_plik_dx11dr2_HM_v18_TT;
	clik_object * clikid_lowl_SMW_70_dx11d_2014;
	clik_error *_err;
	
	void  data_initialize()
	{
 	  char clik_hpath[] = "/path/to/plik_dx11dr2_HM_v18_TT.clik";
 	  char clik_lpath[] = "/path/to/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik";
 
	  _err = initError();
	  clikid_plik_dx11dr2_HM_v18_TT = clik_init(*&clik_hpath,&_err);
	  clikid_lowl_SMW_70_dx11d_2014 = clik_init(*&clik_lpath,&_err);
 
 	  // cout << "initialized the data files" << endl;
 
	}
	 
	void  data_cleanup()
	{
	  clik_cleanup(byVal(&clikid_plik_dx11dr2_HM_v18_TT));
	  clik_cleanup(byVal(&clikid_lowl_SMW_70_dx11d_2014));
		
      // cout << "cleaned-up the data files" << endl;
	}
	 
	clik_object * return_clikid_plik_dx11dr2_HM_v18_TT()
	{
 	  return clikid_plik_dx11dr2_HM_v18_TT;
	}
	
	clik_object * return_clikid_lowl_SMW_70_dx11d_2014()
	{
 	  return clikid_lowl_SMW_70_dx11d_2014;
	}
}
END_BE_NAMESPACE

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
