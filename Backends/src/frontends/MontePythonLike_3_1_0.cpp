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
  // Empty ini function.
	using namespace pybind11::literals; // to bring in the `_a` literal
	// py::dict kwargs = py::dict("number"_a=1234, "say"_a="hello", "to"_a=some_instance);
	std::cout << "   (MontePythonLike init): begining ini function "<< std::endl;

	// These things all need to be passed from CosmoBit
	pybind11::str root = "/home/janina/code/gambit_devel/Backends/installed/";
	pybind11::dict path_dict = pybind11::dict("MontePython"_a="/home/janina/code/gambit_devel/Backends/installed/montepythonlike/3.1.0/montepython/", "data"_a="/home/janina/code/gambit_devel/Backends/installed/montepythonlike/3.1.0/data/",
		"cosmo"_a="/home/janina/code/gambit_devel/Backends/installed/exoclass/2.7.0/", "root"_a=root);
	//pybind11::isinstance<py::list>(obj)
	auto experiments = pybind11::make_tuple("bao");

	pybind11::dict cosmo_arguments;
	pybind11::dict mcmc_parameters;
	std::cout << "   (MontePythonLike init): before fill data dict "<< std::endl;
	//pybind11::dict data_dict = pybind11::dict("path"_a=path_dict, "cosmo_arguments"_a=cosmo_arguments);
	pybind11::str command_line = "" ;
	

	pybind11::str path_dat = "/home/janina/code/montepython_public/montepython/likelihoods/bao/bao.data";
	
	std::cout << "   (MontePythonLike init): before Likelihood object "<< std::endl;

	pybind11::object data = MontePythonLike.attr("Data")("",path_dat,experiments,mcmc_parameters);
	
	//pybind11::object bao = py::eval("my_variable + 10", scope).cast<int>();

	pybind11::object Likelihood = MontePythonLike.attr("Likelihood")(path_dict, data, command_line);

	//std::cout<< "   (MontePythonLike init): before loglkl "<< std::endl;
	//pybind11::dbl chi;
	Likelihood.attr("loglkl")("","");

	//std::cout << "   (MontePythonLike init): chi is " << chi << std::endl;
}
END_BE_INI_FUNCTION

/// TODO -
/// Set nuisance parameters for scans within directdm (BE_INI_FUNCTION)

BE_NAMESPACE
{
  
  void test_MontePythonLike()
  {
  	//pybind11::str like_name = ""
	std::cout << "   		(MontePythonLike test) "  << std::endl;
  }

}
END_BE_NAMESPACE
