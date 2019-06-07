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

	// S.B. TODO: Hmmm - somehow need to add a dependency on ExoClass/Class/etc, if the 
	// directory needs to live in the path_dict object.
	
}
END_BE_INI_FUNCTION

BE_NAMESPACE
{
  using namespace pybind11::literals; // to bring in the `_a` literal

  /// Initialise MontePython.
  std::vector<double> get_MP_loglikes(pybind11::object& cosmo)
  {
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

	// Maybe make a dictionary?
	// dict = {
	//			"bao":      "/bao/bao.data",
	//			"Pantheon": "/Pantheon/Pantheon.data"
	//			...
	//		   }
	// In fact we can make this a macro/function, probably, since
	// it looks like the experiment + likelihood share a name in MP. Best day ever!

	// In the future, this will be an array of experiments used in the scan.
	// For now, just do BAOs
	pybind11::tuple experiments = pybind11::make_tuple("bao");
	// Location of the BAO data
	pybind11::str bao_path = like_path +"/bao/bao.data";

	// Import Data object from MontePython
	pybind11::object data = 
		MontePythonLike.attr("Data")("",path_dict,experiments,mcmc_parameters);

	// Pass the Data object to the Likelihood object internal to MontePython.
	// do we need this?
	//pybind11::object Likelihood = 
	//	MontePythonLike.attr("Likelihood")(bao_path, data, command_line);

	// Add the path to the likelihoods in MP to sys.path.
	pybind11::module sys = pybind11::module::import("sys");
	sys.attr("path").attr("append")(like_path);

	// Now we can import the bao likelihood module...
	pybind11::module bao = pybind11::module::import("bao");

	// Create a BAO object
	pybind11::object BAO = 
	    bao.attr("bao")(bao_path, data, command_line);

	// Just check things are working as expected: should return 3
    std::cout << "num_points = " << BAO.attr("num_points").cast<int>() << std::endl;

    // Now use the cosmo python object, to compute the likelihood
    // todo: python does not know how to cast the cosmo object yet.
    // this doesn't work.


    // temporary
    sys.attr("path").attr("append")("/Backends/installed/exoclass/2.7.0/python/build/lib.linux-x86_64-2.7/classy.so");
    pybind11::module classy = pybind11::module::import("classy");

    pybind11::object Class = classy.attr("Class")();
    double age = Class.attr("age")().cast<double>();

    std::cout << age << std::endl;

    pybind11::object baologlike = BAO.attr("loglkl")(Class, data);

    /*	1. create instance of Class() object in CLASS frontend
    *   2. pass this back to CosmoBit 
    *   3. pass from CosmoBit to the MPL frontend
    *   4. 
    */

    // class CosmoTest
    // {
    // 	std::string name_;
    // 	int num_;
    // };
    // PYBIND11_MODULE(cc11binds,m) {
    // 	pybind11::class(m,"CosmoTest")
    // 	.def(pybind11::init<str,)
    // }


    //std::cout << "baologlike = " << baologlike << std::endl;
    std::vector<double> loglikes;
    //std::vector<pybind11::object> loglikesobj;
	//loglikesobj.push_back(baologlike);

	return loglikes;
  }

  void test_MontePythonLike(pybind11::object& cosmo)
  {
  	std::vector<double> loglikes = get_MP_loglikes(cosmo);
	std::cout << "   		(MontePythonLike test) "  << std::endl;
  }

}
END_BE_NAMESPACE
