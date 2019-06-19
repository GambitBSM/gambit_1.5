//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend source for the exoclassy backend.
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
#include "gambit/Backends/frontends/classy_exo_2_7_0.hpp"

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <pybind11/eval.h>


BE_NAMESPACE
{

  // Returns a string of the path to the CLASSY object with respect to backendDir.
  std::string path_to_classy()
  {
    std::string path = "classy/exo_2.7.0/";
    return path;
  }

  void classy_create_class_instance(pybind11::object& result)
  {
		result = classy.attr("Class")();
  }
  
  // clean & empty cosmo object (instance of classy class Class), set the input parameters & run CLASS
  void classy_compute(CosmoBit::Classy_cosmo_container& ccc)
  {
  	
    // Clean CLASS (the equivalent of the struct_free() in the `main` of CLASS -- don't want a memory leak, do we
    ccc.cosmo.attr("struct_cleanup")();

    // Actually only strictly necessary when cosmology is changed completely between two different runs
    // but just to make sure nothing's going wrong do it anyways..
    ccc.cosmo.attr("empty")();

  	// set cosmological parameters
  	ccc.cosmo.attr("set")(ccc.cosmo_input_dict);
  	
    std::cout << "		(classy frontend) after set parameters "<< std::endl;

    // run class
  	ccc.cosmo.attr("compute")();
  	
  	// for testing -- keep it for now just in case.. 
  	//double age = ccc.cosmo.attr("age")().cast<double>();
  	//double h = ccc.cosmo.attr("h")().cast<double>();
    //std::cout << "		(classy frontend) computed age to be "<< age << std::endl;
    //std::cout << "		(classy frontend) computed h "<< h << std::endl;
  }

  
}
END_BE_NAMESPACE


BE_INI_FUNCTION
{ 
  // Don't really need to do anything here -- there is a conditional BE_INI dependency on the 'energy_injection_efficiency' 
  // if a decaying DM model is considered. However, we only need this to be executed before exoclass is called and the 
  // dependency resolver will take care of that. 

}
END_BE_INI_FUNCTION