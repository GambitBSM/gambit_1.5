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
///
///  *************************

#include "gambit/Utils/util_types.hpp"

#ifndef __MultiModeCode_types_hpp__
#define __MultiModeCode_types_hpp__

namespace Gambit
{
	// type definition for the multimodecode output.
	typedef struct
	{
		double As;
		double A_iso;
		double A_pnad;
		double A_ent;
		double A_cross_ad_iso;
		double A_bundle;
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
	} gambit_inflation_observables;

}

#endif // defined __MultiModeCode_types_hpp__
