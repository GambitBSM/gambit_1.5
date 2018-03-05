//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Helper types for plc backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Selim C. Hotinli
///          selim.hotinli14@imperial.ac.uk
///  \date 2017 Jul
///
///  *************************

#include "gambit/Utils/util_types.hpp"

#ifndef __plc_types_hpp__
#define __plc_types_hpp__

namespace Gambit
{
	// Types needed for using the Planck likelihoods
	typedef char parname[256];
	
	typedef void clik_object;
	
	typedef struct _err {
		char errWhere[2048];
		char errText[4192];
		int errValue;
		struct _err* next;
	} clik_error;
	
	typedef char parnam[256];
	
	typedef struct {
		void *plens_payload;
		int lmax[7];
		int type;
		int renorm;
		int ren1;
		double check;
		int has_check;
		double *cl_fid;
	} clik_lensing_object;

	typedef int int_6[6];
}

#endif // defined __plc_types_hpp__
