//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions of module CosmoBitDEV.
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
///
///  *********************************************
#include <string>
#include <iostream>
#include <cmath>
#include <functional>
#include <omp.h>
#include <stdlib.h>     /* malloc, free, rand */
#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/CosmoBitDEV/CosmoBitDEV_rollcall.hpp"
#include "gambit/Backends/frontend_macros.hpp"

namespace Gambit
{

	namespace CosmoBitDEV
	{
		using namespace LogTags;
		//safe_ptr(Pipes::function_vanilla_lowp_TT_loglike::			BEreq::class_parser_initialize(&fc,9,"",class_errmsg);
		//blah = return_vanilla_cls(ob, cdm, ..., *BEreq::class_parser_initializ)
		
		double** return_vanilla_cls(double omega_b,double omega_cdm,double H0,double ln10A_s,double n_s,double tau_reio)
		{
			using namespace Pipes::function_vanilla_lowp_TT_loglike;
			
			struct precision pr;        /* for precision parameters */
			struct background ba;       /* for cosmological background */
			struct thermo th;           /* for thermodynamics */
			struct perturbs pt;         /* for source functions */
			struct transfers tr;        /* for transfer functions */
			struct primordial pm;       /* for primordial spectra */
			struct spectra sp;          /* for output spectra */
			struct nonlinear nl;        /* for non-linear spectra */
			struct lensing le;          /* for lensed spectra */
			struct output op;           /* for output files */
			ErrorMsg class_errmsg;      /* for error messages */
			
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
			
			struct file_content fc;
			
			BEreq::class_parser_initialize(&fc,9,"",class_errmsg);
			
			
			strcpy(fc.name[0],"output");
			strcpy(fc.value[0],"tCl,pCl");
			
			strcpy(fc.name[1],"omega_b");
			strcpy(fc.name[2],"omega_cdm");
			strcpy(fc.name[3],"H0");
			strcpy(fc.name[4],"ln10^{10}A_s");
			strcpy(fc.name[5],"n_s");
			strcpy(fc.name[6],"tau_reio");
			
			strcpy(fc.name[7],"l_scalar_max");
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
			strcat(pr.sBBN_file,"Backends/installed/class/2.6.1/bbn/sBBN.dat");
			
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
			
			struct precision pr;        /* for precision parameters */
			struct background ba;       /* for cosmological background */
			struct thermo th;           /* for thermodynamics */
			struct perturbs pt;         /* for source functions */
			struct transfers tr;        /* for transfer functions */
			struct primordial pm;       /* for primordial spectra */
			struct spectra sp;          /* for output spectra */
			struct nonlinear nl;        /* for non-linear spectra */
			struct lensing le;          /* for lensed spectra */
			struct output op;           /* for output files */
			ErrorMsg class_errmsg;      /* for error messages */
			
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
			
			struct file_content fc;
			
			BEreq::class_parser_initialize(&fc,11,"",class_errmsg);
			
			
			strcpy(fc.name[0],"output");
			strcpy(fc.value[0],"tCl,pCl");
			
			strcpy(fc.name[1],"omega_b");
			strcpy(fc.name[2],"omega_cdm");
			strcpy(fc.name[3],"H0");
			strcpy(fc.name[4],"ln10^{10}A_s");
			strcpy(fc.name[5],"n_s");
			strcpy(fc.name[6],"tau_reio");
			
			strcpy(fc.name[8],"l_scalar_max");
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
			strcat(pr.sBBN_file,"Backends/installed/class/2.6.1/bbn/sBBN.dat");
			
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
			
			struct precision pr;        /* for precision parameters */
			struct background ba;       /* for cosmological background */
			struct thermo th;           /* for thermodynamics */
			struct perturbs pt;         /* for source functions */
			struct transfers tr;        /* for transfer functions */
			struct primordial pm;       /* for primordial spectra */
			struct spectra sp;          /* for output spectra */
			struct nonlinear nl;        /* for non-linear spectra */
			struct lensing le;          /* for lensed spectra */
			struct output op;           /* for output files */
			ErrorMsg class_errmsg;      /* for error messages */
			
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
			
			struct file_content fc;
			
			BEreq::class_parser_initialize(&fc,11,"",class_errmsg);
			
			
			strcpy(fc.name[0],"output");
			strcpy(fc.value[0],"tCl,pCl");
			
			strcpy(fc.name[1],"omega_b");
			strcpy(fc.name[2],"omega_cdm");
			strcpy(fc.name[3],"H0");
			strcpy(fc.name[4],"ln10^{10}A_s");
			strcpy(fc.name[5],"n_s");
			strcpy(fc.name[6],"tau_reio");
			
			strcpy(fc.name[8],"l_scalar_max");
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
			strcat(pr.sBBN_file,"Backends/installed/class/2.6.1/bbn/sBBN.dat");
			
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
			//	int nmax;
			int num_ct_max=7;
			
			double	highl_TT_cl_and_pars[2525];
			double	lowp_cl_and_pars[121];
		
			char *class_null;	// defining an empty string for parser initialization.
			
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
			for (l=2;l<=29;l++)	{
				lowp_cl_and_pars[l] = cl[l][0];
			}
			int k;
			lowp_cl_and_pars[30] = 0.0;
			lowp_cl_and_pars[31] = 0.0;
			for (l=32;l<=59;l++)	{
				k = l-30;
				lowp_cl_and_pars[l] = cl[k][1];
			}
			lowp_cl_and_pars[60] = 0.0;
			lowp_cl_and_pars[61] = 0.0;
			for (l=62;l<=89;l++)	{
				k = l-60;
				lowp_cl_and_pars[l] = cl[k][2];
			}
			lowp_cl_and_pars[90] = 0.0;
			lowp_cl_and_pars[91] = 0.0;
			for (l=92;l<=119;l++)	{
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
			lowp_log 	 = BEreq::clik_compute_loglike(byVal(lowl_clikid),
													   byVal(lowp_cl_and_pars),
													   &_err);


			result = highl_TT_log+lowp_log;
			
			std::cout << "Log likelihood is : " << result << std::endl;
		}
		
		
		void function_LCDMtensor_lowp_TT_loglike(double& result)
		{
			using namespace Pipes::function_LCDMtensor_lowp_TT_loglike;
			
			int l,l_max;
			l_max=3000;
			//	int nmax;
			int num_ct_max=7;
			
			double	highl_TT_cl_and_pars[2525];
			double	lowp_cl_and_pars[121];
			
			char *class_null;	// defining an empty string for parser initialization.
			
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
			for (l=2;l<=29;l++)	{
				lowp_cl_and_pars[l] = cl[l][0];
			}
			int k;
			lowp_cl_and_pars[30] = 0.0;
			lowp_cl_and_pars[31] = 0.0;
			for (l=32;l<=59;l++)	{
				k = l-30;
				lowp_cl_and_pars[l] = cl[k][1];
			}
			lowp_cl_and_pars[60] = 0.0;
			lowp_cl_and_pars[61] = 0.0;
			for (l=62;l<=89;l++)	{
				k = l-60;
				lowp_cl_and_pars[l] = cl[k][2];
			}
			lowp_cl_and_pars[90] = 0.0;
			lowp_cl_and_pars[91] = 0.0;
			for (l=92;l<=119;l++)	{
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
			lowp_log 	 = BEreq::clik_compute_loglike(byVal(lowl_clikid),
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
//											   kmax,
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
//			std::cout << "observs.A_bundle " << observs.A_bundle << std::endl;
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
			//	int nmax;
			int num_ct_max=7;
			
			double	highl_TT_cl_and_pars[2525];
			double	lowp_cl_and_pars[121];
			
			char *class_null;	// defining an empty string for parser initialization.
			
			/* allocate the array where calculated Cl's will be written (we
			 could add another array with P(k), or extract other results from
			 the code - here we assume that we are interested in the C_l's
			 only */
			
			double** cl;
			
			cl = return_LCDMtensor_inflation_cls(*Param["omega_b"],
												 *Param["omega_cdm"],
												 *Param["H0"],
												 *Param["ln10A_s"],
//												 *Param["n_s"],
												 observs.ns,
												 *Param["tau_reio"],
//												 *Param["r_tensor"]
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
			for (l=2;l<=29;l++)	{
				lowp_cl_and_pars[l] = cl[l][0];
			}
			int k;
			lowp_cl_and_pars[30] = 0.0;
			lowp_cl_and_pars[31] = 0.0;
			for (l=32;l<=59;l++)	{
				k = l-30;
				lowp_cl_and_pars[l] = cl[k][1];
			}
			lowp_cl_and_pars[60] = 0.0;
			lowp_cl_and_pars[61] = 0.0;
			for (l=62;l<=89;l++)	{
				k = l-60;
				lowp_cl_and_pars[l] = cl[k][2];
			}
			lowp_cl_and_pars[90] = 0.0;
			lowp_cl_and_pars[91] = 0.0;
			for (l=92;l<=119;l++)	{
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
			lowp_log 	 = BEreq::clik_compute_loglike(byVal(lowl_clikid),
													   byVal(lowp_cl_and_pars),
													   &_err);
			
			
			result = highl_TT_log+lowp_log;
			
			std::cout << "Log likelihood is : " << result << std::endl;

			
			
		}
		
	}

}

