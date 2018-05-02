//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for DarkAges 1.0.0 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2018 Mar
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/DarkAges_1_0_0.hpp"

// Convenience functions (definitions)
BE_NAMESPACE
{
  DarkAges::fz_table gather_results()
  {
    DarkAges::fz_table result;
    pybind11::list tmp_py_list;
    logger().send("Message from 'gather_results' backend convenience function in DarkAges v1.0.0 wrapper",LogTags::info);
    tmp_py_list =  get_result("redshift");
    result.redshift = pybind11::cast< std::vector<double> >(tmp_py_list);
    tmp_py_list =  get_result("Heat");
    result.f_heat = pybind11::cast< std::vector<double> >(tmp_py_list);
    tmp_py_list =  get_result("Ly-A");
    result.f_lya = pybind11::cast< std::vector<double> >(tmp_py_list);
    tmp_py_list =  get_result("H-Ion");
    result.f_hion = pybind11::cast< std::vector<double> >(tmp_py_list);
    tmp_py_list =  get_result("He-Ion");
    result.f_heion = pybind11::cast< std::vector<double> >(tmp_py_list);
    tmp_py_list =  get_result("LowE");
    result.f_lowe = pybind11::cast< std::vector<double> >(tmp_py_list);
    
    return result;
  }
}
END_BE_NAMESPACE

// Initialisation function (definition)
BE_INI_FUNCTION
{
  static bool scan_level = true;
  if (scan_level)
  {
    initialize();
  }
  scan_level = false;
  // Reset the 'alreadyCalculated' flag
  //*already_calculated = false;
  release_flag();
  DarkAges::injectionSpectrum spec = *Dep::injection_spectrum;
  pybind11::list Energies = pybind11::cast(spec.E);
  pybind11::list dNdE_e = pybind11::cast(spec.spec_el);
  pybind11::list dNdE_ph = pybind11::cast(spec.spec_ph);
  double mass = *Dep::DM_mass;
  double lifetime = *Dep::lifetime;
  
  calc_f_decay(Energies,dNdE_e,dNdE_ph,mass,lifetime);
}
END_BE_INI_FUNCTION
