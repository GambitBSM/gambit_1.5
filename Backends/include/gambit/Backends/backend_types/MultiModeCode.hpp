//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Helper types for MultiModeCode backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Selim C. Hotinli
///          selim.hotinli14@imperial.ac.uk
///  \date 2017 June
///  \date 2017 Oct
///  *************************

#include "gambit/Utils/util_types.hpp"

#ifndef __MultiModeCode_hpp__
#define __MultiModeCode_hpp__

namespace Gambit
{
	// type definition for the multimodecode output.
	typedef struct
	{
		bool check_ic_ok;
		double As;
		double A_iso;
		double A_pnad;
		double A_ent;
		double A_cross_ad_iso;
//		double A_bundle;
		double ns;
		double nt;
		double n_iso;
		double n_pnad;
		double n_ent;
		double r;
		double alpha_s;
		double runofrun;
		double f_NL;
		double tau_NL;
//		std::vector<double> k_array;  //<- Added for the FULL POW SPEC
//		std::vector<double> pks_array;  //<- Added for the FULL POW SPEC
//		std::vector<double> pkt_array;  //<- Added for the FULL POW SPEC
		double k_array[100];  //<- Added for the FULL POW SPEC
		double pks_array[100];  //<- Added for the FULL POW SPEC
		double pkt_array[100];  //<- Added for the FULL POW SPEC
		int k_size;
	} gambit_inflation_observables;

}

#endif // defined __MultiModeCode_hpp__
