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


BE_NAMESPACE
{
  pybind11::object cosmo;

  // Returns a string of the path to the CLASSY object with respect to backendDir.
  std::string path_to_classy()
  {
    std::string path = "classy/2.6.3/";
    return path;
  }

  // return cosmo object. Need to pass this to MontePython for Likelihoods calculations
  pybind11::object get_classy_cosmo_object()
  {
    return cosmo;
  }

  // getter functions to return a bunch of CLASS outputs. This is here in the frontend
  // to make the capabilities inside CosmoBit independent of types that depend on the 
  // Boltzmann solver in use

  // returns angular diameter distance for given redshift
  double class_get_Da(double z)
  {
    double Da = cosmo.attr("angular_distance")(z).cast<double>();
    // check if units are the same as from class??
    return Da;
  }

  // returns luminosity diameter distance for given redshift
  double class_get_Dl(double z)
  {
    double Dl = cosmo.attr("luminosity_distance")(z).cast<double>();
    return Dl;
  }

  // returns scale_independent_growth_factor for given redshift
  double class_get_scale_independent_growth_factor(double z)
  {
    double growth_fact = cosmo.attr("scale_independent_growth_factor")(z).cast<double>();
    return growth_fact;
  }

  // returns scale_independent_growth_factor for given redshift  TODO: what is different with and without f? think it 
  // is connected to powerspectra with baryons vs. baryons + cdm 
  double class_get_scale_independent_growth_factor_f(double z)
  {
    double growth_fact_f = cosmo.attr("scale_independent_growth_factor_f")(z).cast<double>();
    return growth_fact_f;
  }

  // returns Hubble parameter for given redshift
  double class_get_Hz(double z)
  {
    double H_z = cosmo.attr("Hubble")(z).cast<double>();
    return H_z;
  }

  // returns Omega matter today
  double class_get_Omega0_m()
  {
    double Omega0_m = cosmo.attr("Omega0_m")().cast<double>();
    return Omega0_m;
  }

  // returns Omega nu today
  double class_get_Omega0_nu()
  {
    double Omega0_nu = cosmo.attr("Omega0_nu")().cast<double>();
    return Omega0_nu;
  }

  // returns Omega nu today
  double class_get_Omega0_Lambda()
  {
    double Omega0_Lambda = cosmo.attr("Omega0_Lambda")().cast<double>();
    return Omega0_Lambda;
  }

  // returns sound horizon at drag
  double class_get_rs()
  {
    double rs_d = cosmo.attr("rs_drag")().cast<double>();
    return rs_d;
  }

  // returns sigma8
  double class_get_sigma8()
  {
    // in CosmoBit.cpp test if ClassInput contains mPk -> otherwise SegFault when trying to compute sigma9
    double sigma8 = cosmo.attr("sigma8")().cast<double>();
    return sigma8;
  }

  // returns Neff
  double class_get_Neff()
  {
    // in CosmoBit.cpp test if ClassInput contains mPk -> otherwise SegFault when trying to compute sigma9
    double Neff = cosmo.attr("Neff")().cast<double>();
    return Neff;
  }

}
END_BE_NAMESPACE
  

BE_INI_FUNCTION
{ 

  pybind11::dict cosmo_input_dict = *Dep::get_Classy_cosmo_container;
  
  static bool first_run = true;
  if(first_run)
  {
    cosmo = = classy.attr("Class")();
    first_run = false;
  }

  // Clean CLASS (the equivalent of the struct_free() in the `main` of CLASS -- don't want a memory leak, do we
  cosmo.attr("struct_cleanup")();

  // Actually only strictly necessary when cosmology is changed completely between two different runs
  // but just to make sure nothing's going wrong do it anyways..
  cosmo.attr("empty")();

  // set cosmological parameters
  cosmo.attr("set")(cosmo_input_dict);
  
  std::cout << "    (classy frontend) after set parameters "<< std::endl;

  // run class
  cosmo.attr("compute")();
}
END_BE_INI_FUNCTION
