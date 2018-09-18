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
    pybind11::array_t<double> tmp;
    logger().send("Message from 'gather_results' backend convenience function in DarkAges v1.0.0 wrapper",LogTags::info);
    tmp =  get_result("redshift");
    result.redshift = cast_np_to_std(tmp);
    tmp =  get_result("Heat");
    result.f_heat = cast_np_to_std(tmp);
    tmp =  get_result("Ly-A");
    result.f_lya = cast_np_to_std(tmp);
    tmp =  get_result("H-Ion");
    result.f_hion = cast_np_to_std(tmp);
    tmp =  get_result("He-Ion");
    result.f_heion = cast_np_to_std(tmp);
    tmp =  get_result("LowE");
    result.f_lowe = cast_np_to_std(tmp);

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
  pybind11::array_t<double> Energies = cast_std_to_np(spec.E);
  pybind11::array_t<double> dNdE_e = cast_std_to_np(spec.spec_el);
  pybind11::array_t<double> dNdE_ph = cast_std_to_np(spec.spec_ph);
  double mass = *Dep::DM_mass;
  double lifetime = *Dep::lifetime;

  calc_f_decay(Energies,dNdE_e,dNdE_ph,mass,lifetime);
}
END_BE_INI_FUNCTION
