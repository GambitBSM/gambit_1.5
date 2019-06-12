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
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date 2019 June
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/MontePythonLike_3_1_0.hpp"

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
//#include <pybind11/eigen.h>

// Maybe make a dictionary?
// dict = {
//			"bao":      "/bao/bao.data",
//			"Pantheon": "/Pantheon/Pantheon.data"
//			...
//		   }
// In fact we can make this a macro/function, probably, since
// it looks like the experiment + likelihood share a name in MP. Best day ever!
// (JR) xD hahaha, that is a great idea!! I'll be on it as soon as I've dealt with 
// the other stuff (and after popping the bottle ;)

BE_NAMESPACE
{
  using namespace pybind11::literals; // to bring in the `_a` literal
    
  double get_MP_loglike(const CosmoBit::MPLike_data_container& mplike, pybind11::object& cosmo)
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

  /// Creates a MontePython 'Data' object. 
  ///This is initialised with a list of the relevant experimental limits to import. 
  pybind11::object create_data_object(std::vector<std::string>& experiments)
  {

  	pybind11::dict path_dict = pybind11::dict("MontePython"_a=backendDir+"/montepython/",
  											  "data"_a=backendDir+"/data/",
  											  "cosmo"_a=backendDir+"/../../classy/2.6.3/", 
  											  "root"_a=backendDir+"/../../");

  	pybind11::dict mcmc_parameters;  // Empty - we do our own sampling, cheers.

  	/*
  	// In the future, this will be an array of experiments used in the scan.
  	// For now, just do BAOs
  	pybind11::tuple experiments = pybind11::make_tuple("bao");
  	// Location of the BAO data
  	pybind11::str bao_path = like_path +"/bao/bao.data";
  	*/

  	// Cast the list of experiments to a tuple, for MP to fire up...
  	pybind11::tuple MP_experiments = pybind11::make_tuple(experiments);

  	// Let's just check, shall we...?
  	pybind11::print(MP_experiments);

  	// Import Data object from MontePython
  	std::cout << "   		(MPLike init_MPLike_Likelihoods) About to init data object" << std::endl;
  	pybind11::object data = MontePythonLike.attr("Data")("", path_dict, MP_experiments, mcmc_parameters);
  	std::cout << "   		(MPLike init_MPLike_Likelihoods) Data object initialised!" << std::endl;

  	return data;
  }

/*  //std::map<std::string, pybind11::object> create_likelihood_objects(pybind11::object  & data)
  map_str_dbl create_likelihood_objects(pybind11::object& data)
  {

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
  	likelihoods["bao"] = 420.;

  	std::cout << "   		(MPLike init_MPLike_Likelihoods) finished Likelihood init" << std::endl;
  	return likelihoods;
  		
  }  */
  // v2
  
  map_str_dbl create_likelihood_objects(pybind11::object& data, std::vector<std::string>& experiments)
  {

  	pybind11::str command_line = ""; 

  	// Root likelihood path.
  	std::string like_path = backendDir+"/montepython/likelihoods/";

	// Add the Likelihood path to sys so we can import it in Python
  	pybind11::module sys = pybind11::module::import("sys");
  	sys.attr("path").attr("append")(like_path);

  	//std::map<std::string, pybind11::object> likelihoods;
  	map_str_dbl likelihoods;

  	// Now go through each experiment one by one, and initialise the Likelihood containers in
  	// MontePython, then add them to a dictionary to pass back to CosmoBit.
  	double i = 0.;
  	for (std::vector<std::string>::const_iterator it = experiments.begin(); it != experiments.end(); ++it)
  	{
  		std::string exp_name = *it;
  		//std::transform(exp_name.begin(), exp_name.end(), exp_name.begin(), ::tolower);
		pybind11::str     exp_path = like_path + "/" + exp_name + "/" + exp_name + ".data";
		pybind11::module  exp_module = pybind11::module::import(exp_name.c_str());
		pybind11::object  EXP_MODULE = exp_module.attr(exp_name.c_str())(exp_path, data, command_line);
		likelihoods[exp_name] = 420. + i;
		i++;
  	}

  	std::cout << "   		(MPLike init_MPLike_Likelihoods) finished Likelihood init" << std::endl;
  	return likelihoods;
  		
  }
}
END_BE_NAMESPACE


BE_INI_FUNCTION
{	
	
	
}
END_BE_INI_FUNCTION