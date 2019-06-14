//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend source for the MontePython backend.
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

BE_NAMESPACE
{
  using namespace pybind11::literals; // to bring in the `_a` literal
    
  /// Convenience function to compute the loglike from a given experiment, given a MontePython likelihood-data container
  /// mplike, using the CLASS Python object cosmo.
  double get_MP_loglike(const CosmoBit::MPLike_data_container& mplike, pybind11::object& cosmo, std::string& experiment)
  {

    std::cout << "   		(MontePythonLike) before calling loglkl " << std::endl;

    // need to use likelihood.at() since it is a const map -> [] can create entry & can't be used on const object
  	double result = mplike.likelihoods.at(experiment).attr("loglkl")(cosmo, mplike.data).cast<double>();
    
    std::cout << "   		(MontePythonLike) computed "<< experiment <<" loglike to be " << result << std::endl;
    
    return result;
  }

  /// Creates a MontePython 'Data' object. 
  /// This is initialised with a list of the relevant experimental limits to import. 
  /// Also needs to know where CLASSY lives. 
  pybind11::object create_data_object(std::vector<std::string>& experiments, std::string& classyDir)
  {
  	pybind11::dict path_dict = pybind11::dict("MontePython"_a=backendDir,
  											  "data"_a=backendDir+"/../data/",
  											  "cosmo"_a=backendDir+"/../../../"+classyDir, 
  											  "root"_a=backendDir+"/../../../");

  	// TODO nuisance parameters and other Cosmology to go here...
  	pybind11::dict mcmc_parameters;

  	// Cast the list of experiments to a tuple, for MP to fire up...
  	pybind11::tuple MP_experiments = pybind11::make_tuple(experiments);

  	// Import Data object from MontePython
  	std::cout << "   		(MPLike init_MPLike_Likelihoods) About to init data object" << std::endl;
  	pybind11::object data = MontePythonLike.attr("Data")("", path_dict, MP_experiments, mcmc_parameters);
  	std::cout << "   		(MPLike init_MPLike_Likelihoods) Data object initialised!" << std::endl;

  	return data;
  }
  
  /// Returns a map of string to Python objects. These Python objects from MontePython are the initialised
  /// "Likelihood" objects used internal to MontePython. GAMBIT interfaces to these by requesting the 
  /// loglike methods of these Python objects.
  map_str_pyobj create_likelihood_objects(pybind11::object& data, std::vector<std::string>& experiments)
  {

  	pybind11::str command_line = ""; 

  	// Root likelihood path.

  	std::string like_path = backendDir+"/likelihoods/";

	// Add the Likelihood path to sys so we can import it in Python
  	pybind11::module sys = pybind11::module::import("sys");
  	sys.attr("path").attr("append")(like_path);

  	map_str_pyobj likelihoods;

  	// Now go through each experiment one by one, and initialise the Likelihood containers in
  	// MontePython, then add them to a dictionary to pass back to CosmoBit.
  	for (std::vector<std::string>::const_iterator it = experiments.begin(); it != experiments.end(); ++it)
  	{
  		std::string exp_name = *it;
		pybind11::str     exp_path = like_path + "/" + exp_name + "/" + exp_name + ".data";
		pybind11::module  exp_module = pybind11::module::import(exp_name.c_str());
		pybind11::object  EXP_MODULE = exp_module.attr(exp_name.c_str())(exp_path, data, command_line);
		likelihoods[exp_name] = EXP_MODULE;
  	}

  	std::cout << "   		(MPLike init_MPLike_Likelihoods) finished Likelihood init" << std::endl;
  	return likelihoods;
  		
  }
}
END_BE_NAMESPACE


BE_INI_FUNCTION
{}
END_BE_INI_FUNCTION