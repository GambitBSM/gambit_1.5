//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for DarkAges 1.2.0 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2019 Oct, Nov
///  \date 2020 Jan
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/DarkAges_1_2_0.hpp"

#include "gambit/Utils/numerical_constants.hpp"

#include <algorithm>

// Convenience functions (definitions)
BE_NAMESPACE
{
  using Backends::cast_std_to_np;
  using Backends::cast_np_to_std;

  namespace py = pybind11;

  template<typename T = double>
  using pyArray = typename py::array_t<T>;
  using pyArray_dbl = pyArray<>;

  static bool alreadyCalculated = false;
  static bool hasAnnihilation = false;
  static bool hasDecay = false;
  static bool f_eff_mode;

  static size_t z_size;
  static double zmax;

  static py::module DA;
  static py::module DA_model;

  static pyArray_dbl redshift;
  static std::map<std::string, py::object> transfer_functions;
  static std::map<std::string, pyArray_dbl> result_map;

  inline pyArray_dbl repeat_front_and_end(pyArray_dbl input, double front, double back)
  {
    return Backends::merge(front, Backends::merge(input, back));
  }

  inline pyArray_dbl repeat_front_and_end(pyArray_dbl input)
  {
    size_t len = input.size();
    double front = (input.attr("__getitem__")(0)).cast<double>();
    double back = (input.attr("__getitem__")(len-1)).cast<double>();

    return repeat_front_and_end(input,front,back);
  }

  void calc_f(DarkAges::Energy_injection_spectrum spec, double mass, double lifetime)
  {
    if (alreadyCalculated)
      return;

    static py::module np = py::module::import("numpy");

    static auto transform_spectrum = [](double& x){x *= 1e-9;};
    static auto transform_energy = [](double& x){x = log10(1e9*x);};
    static auto isNonZero = [](double& x){return abs(x) > std::numeric_limits<double>::epsilon();};

    std::for_each(spec.E_el.begin(), spec.E_el.end(), transform_energy);
    std::for_each(spec.spec_el.begin(), spec.spec_el.end(), transform_spectrum);
    bool hasElectrons = std::any_of(spec.spec_el.begin(), spec.spec_el.end(), isNonZero);
    pyArray_dbl logE_el = cast_std_to_np(spec.E_el);
    pyArray_dbl spec_el = cast_std_to_np(spec.spec_el);
    pyArray_dbl null_el = np.attr("zeros_like")(spec_el);

    std::for_each(spec.E_ph.begin(), spec.E_ph.end(), transform_energy);
    std::for_each(spec.spec_ph.begin(), spec.spec_ph.end(), transform_spectrum);
    bool hasPhotons = std::any_of(spec.spec_ph.begin(), spec.spec_ph.end(), isNonZero);
    pyArray_dbl logE_ph = cast_std_to_np(spec.E_ph);
    pyArray_dbl spec_ph = cast_std_to_np(spec.spec_ph);
    pyArray_dbl null_ph = np.attr("zeros_like")(spec_ph);

    mass *= 1.e9;
    const double me_eV = 1.e9 * m_electron;

    static py::object model_el;
    static py::object model_ph;
    if (hasAnnihilation)
    {
      py::object mod = DA_model.attr("annihilating_model");
      model_el = mod(spec_el, null_el, null_el,
                     mass-me_eV, logE_el, redshift,
                     py::arg("normalize_spectrum_by") = "mass");
      model_ph = mod(null_ph, spec_ph, null_ph,
                     mass, logE_ph, redshift,
                     py::arg("normalize_spectrum_by") = "mass");
    }
    else if (hasDecay)
    {
      py::object mod = DA_model.attr("decaying_model");
      model_el = mod(spec_el, null_el, null_el,
                     mass-2*me_eV, lifetime, logE_el, redshift,
                     py::arg("normalize_spectrum_by") = "mass");
      model_ph = mod(null_ph, spec_ph, null_ph,
                     mass, lifetime, logE_ph, redshift,
                     py::arg("normalize_spectrum_by") = "mass");
    }

    for (auto const& it : transfer_functions)
    {
      std::string channel = it.first;
      py::object tf = it.second;
      pyArray_dbl f = np.attr("zeros_like")(redshift);
      if (hasElectrons)
        f = f.attr("__add__")((model_el.attr("calc_f")(tf)).attr("__getitem__")(-1));
      if (hasPhotons)
        f = f.attr("__add__")((model_ph.attr("calc_f")(tf)).attr("__getitem__")(-1));
      //f = f.attr("__getslice__")(0,z_size-1);
      result_map[channel] = repeat_front_and_end(f);
    }

    // Get an independent copy of "redshift" to modify
    // (Shift entries by one to get z and not z+1)
    pyArray_dbl z(z_size,redshift.data());
    //z = z.attr("__getslice__")(0,z_size-1);
    std::for_each(z.mutable_data(), z.mutable_data()+z.size(), [](double& x){x -= 1;});
    result_map["redshift"] =  repeat_front_and_end(z,0,zmax);

    alreadyCalculated = true;

    return;
  }

  DarkAges::Energy_injection_efficiency_table get_energy_injection_efficiency_table()
  {
    DarkAges::Energy_injection_efficiency_table result{};
    logger().send("Message from 'gather_results' backend convenience function in DarkAges v1.0.0 wrapper (Start)",LogTags::info);

    auto safe_retrieve = [](const str& key) -> std::vector<double>
    {
       if (result_map.count(key) != 0)
       {
         return cast_np_to_std(result_map[key]);
       }
       else
       {
         std::ostringstream errMssg;
         errMssg << "Could not successfully gather the results ";
         errMssg << "of DarkAges_v"<< STRINGIFY(VERSION)<<".\n";
         errMssg << "The key -> " << key << " <- is not in \'result map\'.";
         backend_error().raise(LOCAL_INFO,errMssg.str());

         // This will never be returned, but the compiler keeps quiet
         return std::vector<double>{};
       }
    };

    result.f_eff_mode = f_eff_mode;

    result.redshift = safe_retrieve("redshift");

    if (f_eff_mode)
    {
      result.f_eff = safe_retrieve("effective");
    }
    else
    {
      result.f_heat = safe_retrieve("Heat");
      result.f_lya = safe_retrieve("Ly-A");
      result.f_hion = safe_retrieve("H-Ion");
      result.f_heion = safe_retrieve("He-Ion");
      result.f_lowe = safe_retrieve("LowE");
    }

    logger().send("Message from 'gather_results' backend convenience function in DarkAges v1.0.0 wrapper (Done casting NumPy-arrays into std vectors)",LogTags::info);

    return std::move(result);
  }
}
END_BE_NAMESPACE

// Initialisation function (definition)
BE_INI_FUNCTION
{
  result_map.clear();
  static bool scan_level = true;
  static bool print_table = false;
  if (scan_level)
  {
    // Check if annihilating DM or decaying DM is considered.
    // If both is considered, raise an error.
    // The BE_ALLOW_MODELS-macro in the frontend header ensures that
    // hasAnnihilation and hasDecay cannot be false at the same time.
    hasAnnihilation = ModelInUse("AnnihilatingDM_general");
    hasDecay = ModelInUse("DecayingDM_general");
    if (hasAnnihilation && hasDecay)
    {
      std::ostringstream errMssg;
      errMssg << "You asked for a combined scenario of decaying DM and annihilating DM. ";
      errMssg << "The current version of DarkAges (v"<< STRINGIFY(VERSION)<<") cannot handle this.\n";
      errMssg << "Please consider only one scenario exclusively.";
      backend_error().raise(LOCAL_INFO,errMssg.str());
    }

    DA = Gambit::Backends::DarkAges_1_2_0::DarkAges;
    std::string module_name = DA.attr("__name__").cast<std::string>();
    DA_model = py::module::import( (module_name + ".model").c_str() );

    redshift = DA.attr("get_redshift")();
    z_size = redshift.size();

    f_eff_mode = runOptions->getValueOrDef<bool>(false,"f_eff_mode");
    if (f_eff_mode)
    {
      py::object tf_vec = DA.attr("transfer_functions");
      py::object tf_corr = DA.attr("transfer_functions_corr");

      py::object tf_eff = tf_vec.attr("sum")();
      tf_eff = tf_eff.attr("__sub__")(tf_corr);
      transfer_functions["effective"] = tf_eff;
    }
    else
    {
      map_str_int channel_map = DA.attr("channel_dict").cast<map_str_int>();
      for (auto const& it: channel_map)
      {
        std::string channel = it.first;
        int idx = it.second;
        py::object tf = DA.attr("transfer_functions").attr("__getitem__")(idx);
        transfer_functions[channel] = tf;
      }
    }

    zmax = runOptions->getValueOrDef<double>(1e7,"z_max");

    // Should the f(z) be printed (for debugging)
    print_table = runOptions->getValueOrDef<bool>(false,"print_table");
  }
  scan_level = false;

  // Reset the 'alreadyCalculated' flag
  alreadyCalculated = false;
  logger().send("Message from \'DarkAges_1.2.0_ini\'. Retrieving the injection spectrum",LogTags::info);
  DarkAges::Energy_injection_spectrum spec = *Dep::energy_injection_spectrum;
  double mass{};
  double lifetime{};
  if (hasDecay)
  {
    const ModelParameters& NP_params = *Dep::DecayingDM_general_parameters;
    mass = NP_params.at("mass");
    lifetime = NP_params.at("lifetime");
  }
  else if (hasAnnihilation)
  {
    const ModelParameters& NP_params = *Dep::AnnihilatingDM_general_parameters;
    mass = NP_params.at("mass");
    lifetime  = -1.0; // This input will be ignored anyway if hasAnnihilation is true
  }

  logger().send("Message from \'DarkAges_1.2.0_ini\'. Starting now to execute \'calc_f\'.",LogTags::info);
  calc_f(spec,mass,lifetime);
  logger().send("Message from \'DarkAges_1.2.0_ini\'. Done executing \'calc_f\'.",LogTags::info);

  if (print_table)
  {
    std::ostringstream buff;
    buff << "---------------" << "\n";
    if (f_eff_mode)
      buff << "z\tf_eff" << "\n";
    else
      buff << "z\tf_heat\tf_lya\tf_hion\tf_heion\tf_lowe" << "\n";
    for (unsigned int i = 0; i < result_map["redshift"].size(); i++)
    {
      buff << result_map["redshift"].at(i) << "\t";
      if (f_eff_mode)
      {
        buff << result_map["effective"].at(i) << "\n";
      }
      else
      {
        buff << result_map["Heat"].at(i) << "\t";
        buff << result_map["Ly-A"].at(i) << "\t";
        buff << result_map["H-Ion"].at(i) << "\t";
        buff << result_map["He-Ion"].at(i) << "\t";
        buff << result_map["LowE"].at(i)  << "\n";
      }
    }
    buff << "---------------" << "\n";
    std::cout << buff.str();
  }
}
END_BE_INI_FUNCTION
