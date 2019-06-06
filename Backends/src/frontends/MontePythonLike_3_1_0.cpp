//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend source for the DirectDM backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2019 June
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/MontePythonLike_3_1_0.hpp"

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
//#include <pybind11/eigen.h>

BE_INI_FUNCTION
{
	using namespace pybind11::literals; // to bring in the `_a` literal
	std::cout << "   (MontePythonLike init): begining ini function "<< std::endl;

	// S.B. TODO: Hmmm - somehow need to add a dependency on ExoClass/Class/etc, if the 
	// directory needs to live in the path_dict object. 

	pybind11::dict path_dict = pybind11::dict("MontePython"_a=backendDir+"/montepython/",
											  "data"_a=backendDir+"/data/",
											  "cosmo"_a=backendDir+"/../../exoclass/2.7.0/", 
											  "root"_a=backendDir+"/../../");

	// S.B: So this mcmc_parameters variable will always be empty, right?
	// Since we never use MP itself to scan...
	pybind11::dict mcmc_parameters;
	pybind11::str command_line = ""; // And this?

	// Root likelihood path.
	std::string like_path = backendDir+"/montepython/likelihoods/";
	
	// In the future, this will be an array of experiments used in the scan.
	// For now, just do BAOs
	pybind11::tuple experiments = pybind11::make_tuple("bao");
	// Location of the BAO data
	pybind11::str bao_path = like_path +"/bao/bao.data";

	// Maybe make a dictionary?
	// dict = {
	//			"bao":      "/bao/bao.data",
	//			"Pantheon": "/Pantheon/Pantheon.data"
	//			...
	//		   }
	// In fact we can make this a macro/function, probably, since
	// it looks like the experiment + likelihood share a name in MP. Best day ever!
	
	
	// Import Data object from MontePython
	pybind11::object data = 
		MontePythonLike.attr("Data")("",path_dict,experiments,mcmc_parameters);
	
	// // Pass the Data object to the Likelihood object internal to MontePython.
	pybind11::object Likelihood = 
		MontePythonLike.attr("Likelihood")(bao_path, data, command_line);

	pybind11::module sys = pybind11::module::import("sys");
	sys.attr("path").attr("append")(like_path);

	pybind11::module bao = pybind11::module::import("bao");

	pybind11::object BAO = 
		bao.attr("bao")(bao_path, data, command_line);

	std::cout << "num_points = " << BAO.attr("num_points").cast<int>() << std::endl;

	//std::cout<< "   (MontePythonLike init): before loglkl "<< std::endl;
	//pybind11::dbl chi;
	//std::cout << BAO_like.attr("loglkl")("","") << std::endl;

}
END_BE_INI_FUNCTION

BE_NAMESPACE
{
  
  // std::vector<pybind11::object> Likelihoods init_MP_likes()
  // {
  // 	...

  // }

  void test_MontePythonLike()
  {
  	//pybind11::str like_name = ""
	//auto chi = Likelihood.attr("loglkl")("","");
	//LikeObj.attr("loglkl")("","");
	std::cout << "   		(MontePythonLike test) "  << std::endl;
  }

}
END_BE_NAMESPACE
