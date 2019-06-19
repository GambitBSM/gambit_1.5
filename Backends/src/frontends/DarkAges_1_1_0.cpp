//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for DarkAges 1.1.0 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2019 Apr
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/DarkAges_1_1_0.hpp"

// Convenience functions (definitions)
BE_NAMESPACE
{
  // Helper functions (translations std::vector <-> numpy array)
  pybind11::array_t<double> cast_std_to_np(const std::vector<double> input)
  {
    // Get size of input
    int size = input.size();

    // Create pybind11::array_t<double> and get the pointer to its data
    pybind11::array_t<double> output(size);
    auto output_buffer = output.request();
    double *output_ptr    = (double *) output_buffer.ptr;

    // copy std::vector -> pybind11::array_t<double>
    std::memcpy(output_ptr,input.data(),size*sizeof(double));

    return output;
  }

  std::vector<double> cast_np_to_std(const pybind11::array_t<double> input)
  {
    // Get size of input
    int size = input.size();

    // Create std::vector
    std::vector<double> output(size);

    // copy py::array -> std::vector
    std::memcpy(output.data(),input.data(),size*sizeof(double));

    return output;
  }

  DarkAges::fz_table gather_results()
  {

    DarkAges::fz_table result;
    logger().send("Message from 'gather_results' backend convenience function in DarkAges v1.0.0 wrapper (Start)",LogTags::info);

    // Using result.<vector>.at(0) instead of .front() since it's safer -> .front() on empty vector will cause a SegFault, .at(0)
    // gives proper error message

    result.redshift = cast_np_to_std(get_result("redshift"));
    result.ptrs_to_member_vecs["annihil_coef_xe"] = &result.redshift.at(0);

    result.f_heat = cast_np_to_std(get_result("Heat"));
    result.ptrs_to_member_vecs["annihil_coef_heat"] = &result.f_heat.at(0);

    result.f_lya = cast_np_to_std(get_result("Ly-A"));
    result.ptrs_to_member_vecs["annihil_coef_lya"] = &result.f_lya.at(0);

    result.f_hion = cast_np_to_std(get_result("H-Ion"));
    result.ptrs_to_member_vecs["annihil_coef_ionH"] = &result.f_hion.at(0);

    result.f_heion = cast_np_to_std(get_result("He-Ion"));
    result.ptrs_to_member_vecs["annihil_coef_ionHe"] = &result.f_heion.at(0);

    result.f_lowe = cast_np_to_std(get_result("LowE"));
    result.ptrs_to_member_vecs["annihil_coef_lowE"] = &result.f_lowe.at(0);

    // get lengths of the tables (should be the same for all). needs to be passed to classy
    result.num_lines = result.redshift.size();

    logger().send("Message from 'gather_results' backend convenience function in DarkAges v1.0.0 wrapper (Done casting NumPy-arrays into std vectors)",LogTags::info);


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
  *already_calculated = false;
  logger().send("Message from \'DarkAges_1.1.0_ini\'. Retrieving the injection spectrum",LogTags::info);
  DarkAges::injectionSpectrum spec = *Dep::injection_spectrum;
  pybind11::array_t<double> Energies = cast_std_to_np(spec.E);
  pybind11::array_t<double> dNdE_e = cast_std_to_np(spec.spec_el);
  pybind11::array_t<double> dNdE_ph = cast_std_to_np(spec.spec_ph);
  double mass = *Dep::DM_mass;
  double lifetime = *Dep::lifetime;

  logger().send("Message from \'DarkAges_1.1.0_ini\'. Starting now to execute \'calc_f_decay\'.",LogTags::info);
  calc_f_decay(Energies,dNdE_e,dNdE_ph,mass,lifetime);
  logger().send("Message from \'DarkAges_1.1.0_ini\'. Done executing \'calc_f_decay\'.",LogTags::info);
}
END_BE_INI_FUNCTION
