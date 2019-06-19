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
  DarkAges::fz_table fzt;
  using namespace pybind11::literals; // to bring in the `_a` literal

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

  // This routine is to create a map from the new input parameters (defined for using 
  // exoclass with gambit) to the according pointers to the arrays of the thermodynamics members that need 
  // to be overwritten by the result from DarkAges
  map_str_dblptr classy_set_energy_injection_efficiency_input()
  {
    map_str_dblptr result;
    
    static bool first_run = true;

    int num_lines = fzt.redshift.size();
    std::vector<double> z = fzt.redshift;
    std::vector<double> fh = fzt.f_heat;
    std::vector<double> fly = fzt.f_lya;
    std::vector<double> fhi = fzt.f_hion;
    std::vector<double> fhei = fzt.f_heion;
    std::vector<double> flo = fzt.f_lowe;

    double *annihil_coef_xe,*annihil_coef_heat,*annihil_coef_lya,*annihil_coef_ionH,*annihil_coef_ionHe,*annihil_coef_lowE;
    
    // allocating memory of these pointers here. Later they are passed to CLASS and the pointers of the 
    // respective thermodynamic structure members will be assigned to point to these addresses. 
    // e.g. pth->annihil_coef_xe = annihil_coef_xe
    // CLASS cleans up all structures itself after a run -> the memory is not (and should not) be freed up 
    // within GAMBIT. So if there is a memory leak it should not be caused by this *if* this function is only run 
    // when CLASS is involved. In the original implementation this is made sure by 
    annihil_coef_xe = (double *) malloc(num_lines*sizeof(double));
    annihil_coef_heat = (double *) malloc(num_lines*sizeof(double));
    annihil_coef_lya = (double *) malloc(num_lines*sizeof(double));
    annihil_coef_ionH = (double *) malloc(num_lines*sizeof(double));
    annihil_coef_ionHe = (double *) malloc(num_lines*sizeof(double));
    annihil_coef_lowE = (double *) malloc(num_lines*sizeof(double));


    /* Done in class input.c now
    th.annihil_coef_dd_heat = (double *) malloc(num_lines*sizeof(double));
    th.annihil_coef_dd_lya = (double *) malloc(num_lines*sizeof(double));
    th.annihil_coef_dd_ionH = (double *) malloc(num_lines*sizeof(double));
    th.annihil_coef_dd_ionHe = (double *) malloc(num_lines*sizeof(double));
    th.annihil_coef_dd_lowE = (double *) malloc(num_lines*sizeof(double));
    th.annihil_coef_num_lines = num_lines;
    */

    for (int it = 0; it < num_lines; it++)
    {
      annihil_coef_xe[it] = z.at(it);
      std::cout << "              Filled annihil_coef_xe at "<< it<< " with " << annihil_coef_xe[it]<< std::endl;
      annihil_coef_heat[it] = fh.at(it);
      annihil_coef_lya[it] = fly.at(it);
      annihil_coef_ionH[it] = fhi.at(it);
      annihil_coef_ionHe[it] = fhei.at(it);
      annihil_coef_lowE[it] = flo.at(it);
    }

    //std::cout << "Filled annihil_coef_xe with "<< &annihil_coef_xe[0]<< " and " << &annihil_coef_xe[1]<< std::endl;


    // fill map that 
    //result["annihil_coef_num_lines"] = num_lines;
    result["annihil_coef_xe"] = annihil_coef_xe;
    result["annihil_coef_heat"] = annihil_coef_heat;
    result["annihil_coef_lya"] = annihil_coef_lya;
    result["annihil_coef_ionH"] = annihil_coef_ionH;
    result["annihil_coef_ionHe"] = annihil_coef_ionHe;
    result["annihil_coef_lowE"] = annihil_coef_lowE;

    //pybind11::print("");
    //std::cout << "Trying to print address "<< annihil_coef_xe << result["annihil_coef_xe"]<< std::endl;
    //pybind11::print(" Address of  annihil_coef_xe is ", annihil_coef_xe);
    //pybind11::print("");
    return result;
  }

}
END_BE_NAMESPACE


BE_INI_FUNCTION
{ 
  // check if extra modifications related to energy injections are needed (if no model with decaying DM is scanned
  // the class run will proceed as usual)
  bool has_energy_injection = ModelInUse("TestDecayingDM");

  // get results from DarkAges if needed
  if(has_energy_injection)fzt = *Dep::energy_injection_efficiency;

}
END_BE_INI_FUNCTION