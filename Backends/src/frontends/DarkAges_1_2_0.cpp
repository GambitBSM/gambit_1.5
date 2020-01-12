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

  void calc_f(DarkAges::injectionSpectrum spec, double mass, double lifetime)
  {
    if (alreadyCalculated)
      return;

    static py::module np = py::module::import("numpy");

    static auto transform_spectrum = [](double& x){x *= 1e-9;};
    static auto transform_energy = [](double& x){x = log10(1e9*x);};

    std::for_each(spec.spec_el.begin(), spec.spec_el.end(), transform_spectrum);
    std::for_each(spec.spec_ph.begin(), spec.spec_ph.end(), transform_spectrum);
    std::for_each(spec.E.begin(), spec.E.end(), transform_energy);

    pyArray_dbl logE = cast_std_to_np(spec.E);
    pyArray_dbl spec_el = cast_std_to_np(spec.spec_el);
    pyArray_dbl spec_ph = cast_std_to_np(spec.spec_ph);
    pyArray_dbl spec_oth = np.attr("zeros_like")(spec_el);

    mass *= 1.e9;

    static py::object model;
    if (hasAnnihilation)
    {
      py::object mod = DA_model.attr("annihilating_model");
      model = mod(spec_el, spec_ph, spec_oth,
                  mass, logE, redshift);
    }
    else if (hasDecay)
    {
      py::object mod = DA_model.attr("decaying_model");
      model = mod(spec_el, spec_ph, spec_oth,
                  mass, lifetime, logE, redshift);
    }

    for (auto const& it : transfer_functions)
    {
      std::string channel = it.first;
      py::object tf = it.second;
      pyArray_dbl f = ((model.attr("calc_f")(tf)).attr("__getitem__")(-1));
      f = f.attr("__getslice__")(0,z_size-1);
      result_map[channel] = repeat_front_and_end(f);
    }

    // Get an independent copy of "redshift" to modify
    // (Stripping the last element and shift by one to get z and not z+1)
    pyArray_dbl z(z_size,redshift.data());
    z = z.attr("__getslice__")(0,z_size-1);
    std::for_each(z.mutable_data(), z.mutable_data()+z.size(), [](double& x){x -= 1;});
    result_map["redshift"] =  repeat_front_and_end(z,0,zmax);

    alreadyCalculated = true;

    return;
  }

  DarkAges::fz_table gather_results()
  {
    DarkAges::fz_table result{};
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
      transfer_functions["effective"] = DA.attr("transfer_functions_corr");
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
  }
  scan_level = false;

  // Reset the 'alreadyCalculated' flag
  alreadyCalculated = false;
  logger().send("Message from \'DarkAges_1.2.0_ini\'. Retrieving the injection spectrum",LogTags::info);
  DarkAges::injectionSpectrum spec = *Dep::injection_spectrum;
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
}
END_BE_INI_FUNCTION
