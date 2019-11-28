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
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2019 July
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

  pybind11::object static cosmo;

  // return cosmo object. Need to pass this to MontePython for Likelihoods calculations
  pybind11::object get_classy_cosmo_object()
  {
		return cosmo;
  }

  /// routine to check class input for consistency. If a case is not treated here class will 
  ///  throw an error saying that one input parameter was not read. Checking for some specific cases
  /// here allows us to give more informative error messages and fix suggestions 
  void class_input_consistency_checks(pybind11::dict classy_input)
  {

    std::ostringstream errMssg;
    // one thing that can go wrong is that the primordial power spectrum is requested but 
    // not output of requiring the perturbations to be solved is asked for => 
    // check if "modes" input is set while "output" is not set
    if(classy_input.contains("modes") and not classy_input.contains("output"))
    {
      errMssg << "You are calling class asking for the following modes to be computed : "<< pybind11::repr(classy_input["modes"]);
      errMssg << "\nHowever, you did not request any output that requires solving the perturbations.\nHence CLASS";
      errMssg << " will not read the input 'modes' and won't run. Add the CLASS input parameter 'output' requesting";
      errMssg << " a spectrum to be computed to the yaml file as run option, e.g. \n  - capability: baseline_classy_input\n";
      errMssg << "    options:\n      classy_dict:\n        output: tCl";
      backend_error().raise(LOCAL_INFO,errMssg.str());
    }
  }
  
  // getter functions to return a bunch of CLASS outputs. This is here in the frontend
  // to make the capabilities inside CosmoBit independent of types that depend on the 
  // Boltzmann solver in use

  // get the CLs as they are needed for the Planck likelihood.
  std::vector<double> class_get_cl(std::string spectype)
  {
    // Get dictionary containing all (lensed) Cl spectra
    pybind11::dict cl_dict = cosmo.attr("lensed_cl")();

    // Get only the relevant Cl as np array and steal the pointer to its data.
    pybind11::object cl_array_obj = cl_dict[pybind11::cast<str>(spectype)];
    pybind11::array_t<double> cl_array = pybind11::cast<pybind11::array_t<double>>(cl_array_obj);
    double* cltemp = (double*) cl_array.request().ptr;
    int len = pybind11::cast<int>(cl_array.attr("__len__")());

    // Class calculates dimensionless spectra.
    // To compare with observation (e.g. Planck),
    // TT EE TE BB need to be multiplied by T_CMB^2.
    double factor = 1.;
    if (spectype.compare("pp") != 0)
      factor = pow( 1.e6*(cosmo.attr("T_cmb")()).cast<double>(), 2);

    std::vector<double> result(len,0.);
    // Loop through all l from 0 to len
    for (int l=0; l < len; l++)
    {
      // The entries for l=0 and l=1 are zero per defintion
      if (l < 2){result[l] = 0;}
      else      {result[l] = factor*cltemp[l];}
    }

    return result;
  }

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

  // returns Omega radiation today
  double class_get_Omega0_r()
  {
    double Omega0_r = cosmo.attr("Omega_r")().cast<double>();
    return Omega0_r;
  }

  // returns Omega of ultra-relativistic species today
  double class_get_Omega0_ur()
  {
    double Omega0_ur = cosmo.attr("Omega_ur")().cast<double>();
    return Omega0_ur;
  }

  // returns Omega matter today
  double class_get_Omega0_m()
  {
    double Omega0_m = cosmo.attr("Omega_m")().cast<double>();
    return Omega0_m;
  }

  // returns Omega ncdm today (contains contributions of all ncdm species)
  double class_get_Omega0_ncdm_tot()
  {
    double Omega0_ncdm = cosmo.attr("Omega_ncdm_tot")().cast<double>();
    return Omega0_ncdm;
  }

/* you *could* also have this function in principle, however since the containts
  the contribution from ALL ncdm components it is not in general true, that Omega_ncdm = Omega_nu
  so I would not reccomend doing it. 
  // returns Omega nu today
  double class_get_Omega0_nu()
  {
    double Omega0_nu = cosmo.attr("Omega_nu")().cast<double>();
    return Omega0_nu;
  }
*/
  // returns Omega_Lambda
  double class_get_Omega0_Lambda()
  {
    double Omega0_Lambda = cosmo.attr("Omega_Lambda")().cast<double>();
    return Omega0_Lambda;
  }


  // returns sound horizon at drag
  double class_get_rs()
  {
    double rs_d = cosmo.attr("rs_drag")().cast<double>();
    return rs_d;
  }

  // returns sigma8 at z = 0
  double class_get_sigma8()
  {
    // in CosmoBit.cpp test if ClassInput contains mPk -> otherwise SegFault when trying to compute sigma8
    double sigma8 = cosmo.attr("sigma8")().cast<double>();
    return sigma8;
  }
  
  // returns Neff
  double class_get_Neff()
  {
    // in CosmoBit.cpp test if ClassInput contains mPk -> otherwise SegFault when trying to compute sigma8
    double Neff = cosmo.attr("Neff")().cast<double>();
    return Neff;
  }

  // print primordial power spectrum for consistency check & debug purposes
  void print_pps()
  {
    std::cout<< "Primordial spectrum from classy: "<< std::string(pybind11::str(cosmo.attr("get_primordial")())) << std::endl;
  }

}
END_BE_NAMESPACE

BE_INI_FUNCTION
{ 
  CosmoBit::ClassyInput input_container= *Dep::get_classy_cosmo_container;
  pybind11::dict cosmo_input_dict = input_container.get_input_dict();

  static bool first_run = true;
  if(first_run)
  {
    cosmo = classy.attr("Class")();
    first_run = false;
    // check input for consistency
    class_input_consistency_checks(cosmo_input_dict);
  }

  // Clean CLASS (the equivalent of the struct_free() in the `main` of CLASS -- don't want a memory leak, do we
  cosmo.attr("struct_cleanup")();

  // Actually only strictly necessary when cosmology is changed completely between two different runs
  // but just to make sure nothing's going wrong do it anyways..
  cosmo.attr("empty")();

  // set cosmological parameters
  logger() << LogTags::debug << "[classy_"<< STRINGIFY(VERSION) <<"] These are the inputs:\n\n";
  logger() << pybind11::repr(cosmo_input_dict) << EOM;
  cosmo.attr("set")(cosmo_input_dict);

  // Try to run class and catch potential errors
  logger() << LogTags::info << "[classy_"<< STRINGIFY(VERSION) <<"] Start to run \"cosmo.compute\"" << EOM;
  try
  {
    cosmo.attr("compute")();
  }
  catch (std::exception &e)
  {
    std::ostringstream errMssg;
    errMssg << "Could not successfully execute cosmo.compute() in classy_"<< STRINGIFY(VERSION)<<"\n";
    std::string rawErrMessage(e.what());
    // If the error is a CosmoSevereError raise an backend_error ...
    if (rawErrMessage.find("CosmoSevereError") != std::string::npos)
    {
      errMssg << "Caught a \'CosmoSevereError\':\n\n";
      errMssg << rawErrMessage;
      backend_error().raise(LOCAL_INFO,errMssg.str());
    }
    // .. but if it is 'only' a CosmoComputationError, just invalidate thew parameter point
    else if (rawErrMessage.find("CosmoComputationError") != std::string::npos)
    {
      errMssg << "Caught a \'CosmoComputationError\':\n\n";
      errMssg << rawErrMessage;
      invalid_point().raise(errMssg.str());
    }
    // any other error (which shouldn't occur) gets also caught as invalid point.
    else
    {
      errMssg << "Caught an unspecified error:\n\n";
      errMssg << rawErrMessage;
      cout << "An unspecified error occurred during compute() in classy_"<< STRINGIFY(VERSION) <<":\n";
      cout << rawErrMessage;
      cout << "\n(This point gets invalidated) " << endl;
      invalid_point().raise(errMssg.str());
    }
  }
  logger() << LogTags::info << "[classy_"<< STRINGIFY(VERSION) <<"] \"cosmo.compute\" was successful" << EOM;

}
END_BE_INI_FUNCTION
