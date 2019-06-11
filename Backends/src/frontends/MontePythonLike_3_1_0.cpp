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



BE_NAMESPACE
{
  using namespace pybind11::literals; // to bring in the `_a` literal

    
  double get_MP_loglike(const CosmoBit::MPLike_data_container & mplike, pybind11::object & cosmo)
  {
	//std::string like = "bao"; // add like_string as argument
  	std::string like = "bao";
    std::cout << "   		(MontePythonLike) cosmo h  "<< cosmo.attr("h")().cast<double>() << std::endl;
    std::cout << "   		(MontePythonLike) before calling loglkl " << std::endl;

    // need to use likelihood.at() since it is a const map -> [] can create entry & can't be used on const object
  	double result = mplike.likelihoods.at(like);//.attr("loglkl")(cosmo, ccc.data).cast<double>();
    
    std::cout << "   		(MontePythonLike) computed BAO loglike to be " << result << std::endl;
    
    return result;
  }

  pybind11::object create_data_object()
  {
  	pybind11::dict path_dict = pybind11::dict("MontePython"_a=backendDir+"/montepython/",
  											  "data"_a=backendDir+"/data/",
  											  "cosmo"_a=backendDir+"/../../classy/2.6.3/", 
  											  "root"_a=backendDir+"/../../");
  	pybind11::dict mcmc_parameters;
  	pybind11::str command_line = ""; 

  	// Root likelihood path.
  	std::string like_path = backendDir+"/montepython/likelihoods/";

  	// Maybe make a dictionary?
  	// dict = {
  	//			"bao":      "/bao/bao.data",
  	//			"Pantheon": "/Pantheon/Pantheon.data"
  	//			...
  	//		   }
  	// In fact we can make this a macro/function, probably, since
  	// it looks like the experiment + likelihood share a name in MP. Best day ever!
  	// (JR) xD hahaha, that is a great idea!! I'll be on it as soon as I've dealt with 
  	// the other stuff (and after popping the bottle ;))

  	// In the future, this will be an array of experiments used in the scan.
  	// For now, just do BAOs
  	pybind11::tuple experiments = pybind11::make_tuple("bao");
  	// Location of the BAO data
  	pybind11::str bao_path = like_path +"/bao/bao.data";

  	// Import Data object from MontePython
  	std::cout << "   		(MPLike init_MPLike_Likelihoods) About to init data object" << std::endl;
  	pybind11::object data = MontePythonLike.attr("Data")("",path_dict,experiments,mcmc_parameters);

  	return data;
  }

  //std::map<std::string, pybind11::object> create_likelihood_objects(pybind11::object  & data)
  map_str_dbl create_likelihood_objects(pybind11::object  & data)
  {

  	 pybind11::dict mcmc_parameters;
  	 pybind11::str command_line = ""; 

  	 // Root likelihood path.
  	 std::string like_path = backendDir+"/montepython/likelihoods/";
  	 pybind11::str bao_path = like_path +"/bao/bao.data";

  	pybind11::module sys = pybind11::module::import("sys");
  	sys.attr("path").attr("append")(like_path);

  	// Now we can import the bao likelihood module & create BAO likelihood object
  	pybind11::module bao = pybind11::module::import("bao");

  	std::cout << "   		(MPLike init_MPLike_Likelihoods) About to create BAO" << std::endl;
  	pybind11::object BAO = bao.attr("bao")(bao_path, data, command_line);
  	
  	//std::map<std::string, pybind11::object> likelihoods;
  	map_str_dbl likelihoods;
  	likelihoods["bao"] = 420;

  	std::cout << "   		(MPLike init_MPLike_Likelihoods) finished Likelihood init" << std::endl;
  	return likelihoods;
  		
  }
}
END_BE_NAMESPACE


BE_INI_FUNCTION
{	
	
	
}
END_BE_INI_FUNCTION