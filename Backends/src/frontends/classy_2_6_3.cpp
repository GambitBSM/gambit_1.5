//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend source for the classy backend.
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
#include "gambit/Backends/frontends/classy_2_6_3.hpp"

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <pybind11/eval.h>

BE_INI_FUNCTION
{

}
END_BE_INI_FUNCTION

BE_NAMESPACE
{
  using namespace pybind11::literals; // to bring in the `_a` literal

  // Returns a string of the path to the CLASSY object with respect to backendDir.
  std::string path_to_classy()
  {
    std::string path = "classy/2.6.3/";
    return path;
  }

  void classy_create_class_instance(pybind11::object& result)
  {
		result = classy.attr("Class")();
  }
  
  //void classy_compute(pybind11::object& cosmo, pybind11::dict& cosmo_input_dict)
  void classy_compute(CosmoBit::Classy_cosmo_container& ccc)
  {
  	// (JR) Should MP init before that such that output for class gets filled with all necessary entries
  	
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