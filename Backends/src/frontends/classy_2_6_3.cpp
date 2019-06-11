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
#include "gambit/Backends/frontends/classy_2_6_3.hpp"

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <pybind11/eval.h>

BE_INI_FUNCTION
{
	/*static bool first_init = true;
	if (first_init)
	{
	  // create instance of class Class() from classy
		std::cout << "		(classy frontend) Created cosmo struct in frontend" << std::endl;
		CosmoBit::classy_cosmo_container.cosmo = classy.attr("Class")();
	}
	first_init = false;

	//cosmo = classy.attr("Class")();
	//scan_level = false;
	// cleanup cosmo structure every time backend is called to make sure there is no
	// contamination with computations from previous point (it is okay to be called on empty structure)
	std::cout << "		(classy frontend) About to clean up csomo sturct" << std::endl;
	
	CosmoBit::classy_cosmo_container.cosmo.attr("struct_cleanup")(); 

	std::cout << "		(classy frontend) cleaned up csomo sturct" << std::endl;*/
}
END_BE_INI_FUNCTION

BE_NAMESPACE
{
	void classy_2_6_3_create_python_obj(pybind11::object & result)
	{
		result = classy.attr("Class")();
	}
  
  using namespace pybind11::literals; // to bring in the `_a` literal
  
  //pybind11::object cosmo;
  pybind11::dict combine_py_dicts(pybind11::dict a, pybind11::dict b)
  {
  	pybind11::dict combined_dict = pybind11::eval("dict(a.items() + b.items() +[(k, op(a[k], b[k])) for k in set(b) & set(a)])");
  	return combined_dict;
  }

  void classy_2_6_3_set_parameter(pybind11::object& cosmo ,pybind11::dict& cosmo_input_dict)
  {
  	// (JR) Should MP init before that such that output for class gets filled with all necessary entries
  	
  	// (JR) combinding dicts does not work atm, fix later
  	/*pybind11::dict combined_in_dict = combine_py_dicts(ccc.cosmo_input_dict, ccc.cosmo_prec_dict);
  	//ccc.cosmo.attr("set")(combined_in_dict); */
  	
  	// set cosmological parameters
  	cosmo.attr("set")(cosmo_input_dict);
  	
    std::cout << "		(classy frontend) after set parameters "<< std::endl;

    // run class
  	cosmo.attr("compute")();
  	double age = cosmo.attr("age")().cast<double>();
  	double h = cosmo.attr("h")().cast<double>();

    std::cout << "		(classy frontend) computed age to be "<< age << std::endl;
    std::cout << "		(classy frontend) computed h "<< h << std::endl;
  }


}
END_BE_NAMESPACE



/*void init_class_cosmo(pybind11::object & cosmo)
  {
	cosmo = classy.attr("Class")();
	// cleanup cosmo structure every time backend is called to make sure there is no
	// contamination with computations from previous point (it is okay to be called on empty structure)
	cosmo.attr("struct_cleanup")(); 
  }*/
