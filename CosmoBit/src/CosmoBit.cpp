//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Central module file of CosmoBit.  
///  Calculates cosmology-related observables.
///
///  Additionally, contains main routines for 
///  interfacing to CLASS and MontePython.
///
///  Most of the model- or observable-specific code is
///  stored in separate source files.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Selim C. Hotinli
///          (selim.hotinli14@pimperial.ac.uk)
///  \date 2017 Jul
///  \date 2018 May
///  \date 2018 Aug - Sep
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2017 Nov
///  \date 2018 Jan - May
///  \date 2019 Jan, Feb, June, Nov
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 June
///  \date 2019 Mar,June
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date 2019 June, Nov
///
///  \author Sebastian Hoof
///          (hoof@uni-goettingen.de)
///  \date 2020 Mar
///
///  \author Pat Scott
///          (pat.scott@uq.edu.au)
///  \date 2018 Mar
///  \date 2019 Jul
///  \date 2020 Apr
///
///  *********************************************

#include <stdint.h> // save memory addresses as int
#include <boost/algorithm/string/trim.hpp>

#include "gambit/Utils/ascii_table_reader.hpp"
#include "gambit/Utils/ascii_dict_reader.hpp"
#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/CosmoBit/CosmoBit_rollcall.hpp"
#include "gambit/CosmoBit/CosmoBit_types.hpp"
#include "gambit/CosmoBit/CosmoBit_utils.hpp"

namespace Gambit
{

  namespace CosmoBit
  {
    using namespace LogTags;

    /***********************************/
    /* General cosmological quantities */
    /***********************************/

    /// Function for setting k_pivot in Mpc^-1 for consistent use within CosmoBit
    /// (i.e. ensuring a consistent value is used by both CLASS and MultiModeCode)
    void set_k_pivot(double &result)
    {
      result = Pipes::set_k_pivot::runOptions->getValueOrDef<double>(0.05, "k_pivot");
    }

    void get_mNu_tot(double& result)
    {
      using namespace Pipes::get_mNu_tot;

      // The units of StandardModel_SLHA2 are GeV; here we are using eV.
      result = 1e9 * ( *Param["mNu1"] + *Param["mNu2"] + *Param["mNu3"] );
    }

    // Returns the effective number of ultrarelativistic species today
    void get_N_ur(double& result)
    {
      using namespace Pipes::get_N_ur;

      // The units of StandardModel_SLHA2 are GeV; here we are using eV.
      std::vector<double> nuMasses{
        1e9*(*Param["mNu1"]), 1e9*(*Param["mNu2"]), 1e9*(*Param["mNu3"])
      };

      // Count the nonzero entries
      auto isNonZero = [](double i) {return i > 0.;};
      int N_ncdm = std::count_if(nuMasses.begin(), nuMasses.end(), isNonZero);

      // Assing the result to the standard value of N_ur depending on the number of massive neutrinos.
      switch (N_ncdm)
      {
        case 1:
          result = 2.0328;  // N_ur (today) = 2.0328 for 1 massive neutrino at CMB release
          break;
        case 2:
          result = 1.0196;  // N_ur (today) = 1.0196 for 2 massive neutrino at CMB release
          break;
        case 3:
          result = 0.00641;  // N_ur (today) = 0.00641 for 3 massive neutrinos at CMB release
          break;
        case 0:
          result = 3.046;
          break;
        default:
          {
            std::ostringstream err;
            err << "You are asking for more than three massive neutrino species.\n";
            err << "Such a case is not implemented in CosmoBit. ";
            err << "If you want to consider this you can add it to the function ";
            err << "'get_N_ur' of the capability 'N_ur'.";
            CosmoBit_error().raise(LOCAL_INFO, err.str());
          }
      }

      // If "etaBBN_rBBN_rCMB_dNurBBN_dNurCMB" is in use, the result will
      // be scaled and gets extra contributions
      if (ModelInUse("etaBBN_rBBN_rCMB_dNurBBN_dNurCMB"))
      {
        // Check if the input for dNeff is negative (unphysical)
        // NOTE: CosmoBit performs no sanity checks if you allow negative dNEff; you're on your own.
        static bool allow_negative_delta_N_ur = runOptions->getValueOrDef<bool>(false,"allow_negative_delta_N_ur");

        // Get values of the temperature ratio and any ultrarelativistic contribution.
        const ModelParameters& NP_params = *Dep::etaBBN_rBBN_rCMB_dNurBBN_dNurCMB_parameters;
        double dNurCMB =  NP_params.at("dNur_CMB");
        double rCMB =  NP_params.at("r_CMB");

        // Only let the user have negative contributions to dNEff if they've signed off on it.
        if ( (!allow_negative_delta_N_ur) && (dNurCMB < 0.0) )
        {
          std::ostringstream err;
          err << "A negative value for \"dNur_CMB\" is unphysical and is not allowed in CosmoBit by default!\n\n";
          err << "If you want to proceed with negative values, please add\n\n";
          err << "  - module: CosmoBit\n";
          err << "    options:\n";
          err << "      allow_negative_delta_N_ur: true\n\n";
          err << "to the Rules section of the YAML file.";
          CosmoBit_error().raise(LOCAL_INFO,err.str());
        }

        // If the check is passed, set the result.
        result = pow(rCMB,4)*(result) + dNurCMB;
      }
      logger() << "N_ur calculated to be " << result << EOM;
    }

    /// Temperature of non-CDM in the (cosmological) SM.
    void T_ncdm_SM(double &result)
    {
      using namespace Pipes::T_ncdm_SM;

      // Set to 0.71611 in units of photon temperature, above the instantaneous decoupling value (4/11)^(1/3)
      // to recover Sum_i mNu_i/omega = 93.14 eV resulting from studies of active neutrino decoupling (arXiv:hep-ph/0506164)
      result = 0.71611;
      // This standard value enters in many assumptions entering CLASS. Therefore changing this value in
      // the YAML file is disabled at the moment. If you still want to modify it, uncomment the line below and
      // you can set is as a runOption of this capability.
      // result = runOptions->getValueOrDef<double>(0.71611,"T_ncdm");
    }

    /// Temperature of non-CDM in non-standard theories.
    void T_ncdm(double &result)
    {
      using namespace Pipes::T_ncdm;

      double rCMB = 1.0; // Default value if no energy injection is assumed.

      // If the "etaBBN_rBBN_rCMB_dNurBBN_dNurCMB" model is included in the scan,
      // we use rCMB of this model.
      if (ModelInUse("etaBBN_rBBN_rCMB_dNurBBN_dNurCMB"))
      {
        const ModelParameters& NP_params = *Dep::etaBBN_rBBN_rCMB_dNurBBN_dNurCMB_parameters;
        rCMB = NP_params.at("r_CMB");
      }

      // Take the SM value of T_ncdm (T_nu) and multiply it with the value of rCMB
      result = rCMB*(*Dep::T_ncdm_SM);
    }

    /// Baryon-to-photon ratio in LCDM
    void eta0_LCDM(double &result)
    {
      using namespace Pipes::eta0_LCDM;

      double ngamma, nb;
      ngamma = 16*pi*zeta3*pow(*Param["T_cmb"]*_kB_eV_over_K_/_hc_eVcm_,3); // photon number density today
      nb = *Param["omega_b"]*3*100*1e3*100*1e3/_Mpc_SI_/_Mpc_SI_/(8*pi*_GN_cgs_* m_proton*1e9*eV2g); // baryon number density today

      result =  nb/ngamma;
      logger() << "Baryon to photon ratio (eta) today computed to be " << result << EOM;
    }

    /// Baryon-to-photon ratio
    void etaBBN_SM(double& result)
    {
      using namespace Pipes::etaBBN_SM;

      result = *Dep::eta0;
    }

    /// The total matter content today.
    void compute_Omega0_m(double &result)
    {
      using namespace Pipes::compute_Omega0_m;

      result =(*Dep::Omega0_b) + (*Dep::Omega0_cdm) + (*Dep::Omega0_ncdm);
    }

    /// The total baryon content today.
    void compute_Omega0_b(double &result)
    {
      using namespace Pipes::compute_Omega0_b;

      double h = *Dep::H0/100.;
      result =*Param["omega_b"]/h/h;
    }

    /// The total cold dark matter content today.
    void compute_Omega0_cdm(double &result)
    {
      using namespace Pipes::compute_Omega0_cdm;

      double h = *Dep::H0/100.;
      result =*Param["omega_cdm"]/h/h;
    }

    /// The total photon content today.
    void compute_Omega0_g(double &result)
    {
      using namespace Pipes::compute_Omega0_g;

      double h = *Dep::H0/100.;
      result = (4.*_sigmaB_SI_/_c_SI_*pow(*Param["T_cmb"],4.)) / (3.*_c_SI_*_c_SI_*1.e10*h*h/_Mpc_SI_/_Mpc_SI_/8./pi/_GN_SI_);
    }

    /// Number density of photons today
    void compute_n0_g(double &result)
    {
      using namespace Pipes::compute_n0_g;

      result = 2./pi/pi*zeta3 *pow(*Param["T_cmb"]*_kB_eV_over_K_,3.)/pow(_hP_eVs_*_c_SI_/2./pi,3)/100/100/100; // result per cm^3
    }

    /// The total ultrarelativistic content today.
    void compute_Omega0_ur(double &result)
    {
      using namespace Pipes::compute_Omega0_ur;

      double N_ur = *Dep::N_ur;
      double Omega0_g = *Dep::Omega0_g;
      result = (N_ur)*7./8.*pow(4./11.,4./3.)* Omega0_g;
    }

    /// (JR) delete when CLASS c interface is removed.
    /// @TODO: can this be deleted now? 
    void compute_Omega0_ncdm(double &result)
    {
      using namespace Pipes::compute_Omega0_ncdm;

      double mNu_tot_eV = *Dep::mNu_tot;
      double h = *Dep::H0/100.;

      result = mNu_tot_eV/(93.14*h*h);  // TODO: heads up: explicit assumption of T_ncdm = 0.71611 and T_cmb goes in here. Has to be generalised
    }

    /// The total relativistic content today.
    void compute_Omega0_r(double &result)
    {
      using namespace Pipes::compute_Omega0_r;

      result = *Dep::Omega0_g + (*Dep::Omega0_ur);
    }

    /*****************/
    /* Classy inputs */
    /*****************/

    /// Create a Python dictionary with the inputs that have to be passed to class.
    /// Setting parameters related to (massive) neutrinos & non-CDM components.
    void set_classy_NuMasses_Nur_input(pybind11::dict &result)
    {
      using namespace Pipes::set_classy_NuMasses_Nur_input;

      // Make sure dict is empty
      result.clear();

      // Set number of ultra relativistic species
      result["N_ur"] = *Dep::N_ur;

      // Get the neutrino masses
      // The units of StandardModel_SLHA2 are GeV; here we are using eV.
      std::vector<double> nuMasses{
        1e9*(*Param["mNu1"]), 1e9*(*Param["mNu2"]), 1e9*(*Param["mNu3"])
      };

      // Count the nonzero entries
      auto isNonZero = [](double i) {return i > 0.;};
      int N_ncdm = std::count_if(nuMasses.begin(), nuMasses.end(), isNonZero);

      // if nonzero, a mass & temperature for each species has to be
      // passed to class
      if (N_ncdm > 0.)
      {
        result["N_ncdm"] = N_ncdm;

        std::vector<double> m_ncdm(N_ncdm);
        std::copy_if(nuMasses.begin(), nuMasses.end(), m_ncdm.begin(), isNonZero);

        // NOTE: this explicitly assumed that all non-CDM components have the same temperature!!
        std::vector<double> T_ncdm(N_ncdm,*Dep::T_ncdm);

        // Create one string with m_ncdm masses and
        // T_ncdm temperatures, separated by commas.
        // This matches the input format of CLASS.
        std::ostringstream ss1, ss2;
        std::string separator;
        for (auto x : m_ncdm)
        {
          ss1 << separator << x;
          separator = ",";
        }
        separator = "";
        for (auto x : T_ncdm)
        {
          ss2 << separator << x;
          separator = ",";
        }

        result["m_ncdm"] = ss1.str();
        result["T_ncdm"] = ss2.str();
      }
    }

    /// Set the classy parameters for an LCDM run with a parametrised primordial power spectrum.
    void set_classy_parameters_parametrised_ps(pybind11::dict& result)
    {
      using namespace Pipes::set_classy_parameters_parametrised_ps;

      // Clean the input container
      result.clear();

      // Now need to pass the primordial power spectrum
      result["n_s"] = *Param["n_s"];
      result["ln10^{10}A_s"] = *Param["ln10A_s"];

      // Pass pivot scale of external spectrum to CLASS
      result["P_k_ini type"] = "analytic_Pk";
      result["k_pivot"] = *Dep::k_pivot;

      // If r = 0 only compute scalar modes, else tensor modes as well
      //
      // => Don't explicitly set "modes" to 's' since it defaults to it. If you set it here anyways
      // you won't be able to run CLASS when only requesting background quantities (e.g. for BAO & SNe likelihoods)
      // as the perturbations module won't run and therefore the entry "modes" won't be read.
      if(*Param["r"] == 0){}
      else
      {
        // Don't set 'r' to zero in CLASS dictionary, as it won't be read if no tensor modes are requested
        result["r"] = *Param["r"];
        result["modes"] = "t,s";
      }

      // Set helium abundance
      result["YHe"] = *Dep::helium_abundance;
    }

    /// Set the classy parameters for an LCDM run with an explicit non-parametric primordial power spectrum.
    void set_classy_parameters_primordial_ps(pybind11::dict& result)
    {
      using namespace Pipes::set_classy_parameters_primordial_ps;

      // Clean the input container
      result.clear();

      // Now need to pass the primordial power spectrum
      static Primordial_ps pps{};
      pps = *Dep::primordial_power_spectrum;
      result["modes"] = "t,s";
      result["P_k_ini type"] = "pointer_to_Pk";
      result["k_array"] = memaddress_to_uint(pps.get_k().data());
      result["pks_array"] = memaddress_to_uint(pps.get_P_s().data());
      result["pkt_array"] = memaddress_to_uint(pps.get_P_t().data());
      result["lnk_size" ] = pps.get_vec_size(); // don't hard code but somehow make consistent with multimode
      
      // Pass pivot scale of external spectrum to CLASS
      result["k_pivot"] = *Dep::k_pivot;

      // Set helium abundance
      result["YHe"] = *Dep::helium_abundance;
    }

    /// Create a Python dictionary with the standard inputs that have to be passed
    /// to CLASS: cosmological parameters ([H0/100*theta_s],omega_b,tau_reio,omega_cdm) & add
    /// model-dependent results for N_ur, neutrino masses & helium abundance.
    /// Also read in any extra input options from the YAML file to pass to CLASS.
    void set_classy_input_params(Classy_input &result)
    {
      using namespace Pipes::set_classy_input_params;

      if (ModelInUse("LCDM") and ModelInUse("LCDM_theta"))
      {
        CosmoBit_error().raise(LOCAL_INFO, "You have requested to scan both LCDM and LCDM_theta.\n"
                                           "This is not allowed. Please select one in your YAML file.");
      }

      // Make sure dict is empty
      result.clear();

      // Keep track if it is the first run: if so, perform some
      // extra consistency checks to make sure no contradicting
      // values are in the classy python input dictionary.
      static bool first_run = true;

      // Get the dictionary with inputs for the neutrino masses and merge it
      // into Classy_Input dictionary
      std::string common_keyszero = result.add_dict(*Dep::classy_MPLike_input);
      result.merge_input_dicts(*Dep::classy_NuMasses_Nur_input);
      result.merge_input_dicts(*Dep::classy_primordial_input);

      // @TODO: can this be deleted?
      // No check if something went wrong and some parameters were defined twice
      //if(common_keys != "")
      //{
      //  CosmoBit_error().raise(LOCAL_INFO, "The key(s) '" + common_keys + "' already "
      //          "exists in the CLASSY dictionary. You are probably trying to override a CLASS setting. Check that none "
      //          "of the parameters you pass through your YAML file through runOptions for the capability 'set_classy_input_params' "
      //          "is in contradiction with any settings made via the dependency resolution by CosmoBit in the function '"+__func__+"'.");
      //}

      // Standard cosmological parameters (common to all CDM-like models)
      result.add_entry("T_cmb"      , *Param["T_cmb"]);
      result.add_entry("omega_b"    , *Param["omega_b"]);
      result.add_entry("tau_reio"   , *Param["tau_reio"]);
      result.add_entry("omega_cdm"  , *Param["omega_cdm"]);

      // Depending on parametrisation, pass either Hubble or the acoustic scale
      if (ModelInUse("LCDM")) result.add_entry("H0", *Param["H0"]);
      else result.add_entry("100*theta_s", *Param["100theta_s"]);

      // TODO: need to test if class or exo_class in use! does not work -> (JR) should be fixed with classy implementation
      // -> (JR again) not sure if that is actually true.. need to test.
      if (ModelInUse("DecayingDM_general") || ModelInUse("AnnihilatingDM_general"))
      {
        // Add decaying/annihilating DM-specific options to Python dictionary passed to CLASS (consistency checks only executed in first run).
        std::string common_keystwo = result.add_dict(*Dep::classy_parameters_EnergyInjection);
      }

      // Other CLASS input direct from the YAML file.
      // Check if these are already contained in the input dictionary -- if so throw an error.
      // Only do this for the first run...
      static pybind11::dict yaml_input;
      if(first_run)
      {
        YAML::Node classy_dict;
        if (runOptions->hasKey("classy_dict"))
        {
          classy_dict = runOptions->getValue<YAML::Node>("classy_dict");
          for (auto it=classy_dict.begin(); it != classy_dict.end(); it++)
          {
            std::string name = it->first.as<std::string>();
            std::string value = it->second.as<std::string>();

            // Check if the key exists in the dictionary
            if (not result.has_key(name.c_str()))
            {
              yaml_input[name.c_str()] = value;
            }
            // If it does, throw an error, there's some YAML conflict going on.
            else
            {
              CosmoBit_error().raise(LOCAL_INFO, "The key '" + name + "' already "
                "exists in the CLASSY dictionary. You are probably trying to override a model parameter. If you really"
                "want to do this you should define an extra function to set the class parameters for the model you "
                "are considering.");
            }
          }
          // Make sure that user did not try to pass k_pivot, N_star or P_k_ini type through class dictionary.
          // These are fixed by capabilities to ensure consistent use throughout the code
          if (yaml_input.contains("k_pivot") || yaml_input.contains("N_star"))
          {
            CosmoBit_error().raise(LOCAL_INFO,
                  "You tried to pass 'k_pivot' and/or 'N_star' to CLASS. These values must \n"
                  "be set consistently throughout the code. N_pivot is set automatically\n"
                  "by the assumption of instant reheating, or as an explicit model parameter.\n"
                  "k_pivot can be set by adding\n "
                  "  - capability: k_pivot\n    function: set_k_pivot\n    options:\n"
                  "      k_pivot: 0.02\n"
                  "to the Rules section of your YAML file.");
          }
          if (yaml_input.contains("P_k_ini type"))
          {
            CosmoBit_error().raise(LOCAL_INFO,
              "GAMBIT will take care of setting all CLASS inputs regarding the primordial power spectrum consistently.\n"
              "Please remove the option 'P_k_ini type' for the capability 'classy_baseline_params'.");
          }
        }
      }

      // Add YAML options to python dictionary passed to CLASS; consistency checks only executed on first run
      result.merge_input_dicts(yaml_input);

      // If the Planck likelihood is used, add the following relevant input parameters to the CLASS dictionary:
      // - output: 'lCl, pCl, tCl'
      // - lensing: yes 
      // - non linear: halofit
      // - l_max_scalars: 2508

      // NOTE: this if performed *after* the YAML options are read in, since it would trigger an error of duplicated
      // keys in the event the user specifies one of the options themselves. The 'merge_pybind_dicts' routine will properly
      // deal with concatenating the output values and choosing the maximum passed value for l_max_scalars. Contradictions
      // in the lensing or non-linear choice will rightfully trigger an error.
      if (ModelInUse("cosmo_nuisance_Planck_lite") || ModelInUse("cosmo_nuisance_Planck_TTTEEE") || ModelInUse("cosmo_nuisance_Planck_TT"))
      {
        // add Planck-likelihood-specific options to python dictionary passed to CLASS; consistency checks only executed on first run
        result.merge_input_dicts(*Dep::classy_PlanckLike_input);
      }

      // Only want to do the gnarly stuff once!
      first_run = false;
    }

    /// Initialises the container within CosmoBit from classy. This holds
    /// an instance of the classy class Class() (Yep, I know...)
    /// which can be handed over to MontePython, or just used to compute
    /// some observables.
    void set_classy_input_no_MPLike(pybind11::dict& result)
    {
      using namespace Pipes::set_classy_input_no_MPLike;

      // Make sure nobody is trying to use MP downstream (prevents segfaults).
      if (Downstream::neededFor("MP_LogLikes"))
      {
        std::ostringstream ss;
        ss << "Sorry, you cannot use the function CosmoBit::set_classy_input with MontePython." << endl
           << "Please modify the Rules section of your YAML file to instead point to the function" << endl
           << "CosmoBit::set_classy_input_with_MPLike. An appropriate rule would look like this:" << endl << endl
           << "  - capability: classy_final_input" << endl
           << "    function: set_classy_input_with_MPLike" << endl;
        CosmoBit_error().raise(LOCAL_INFO, ss.str());
      }

      // If the above check passes, nothing to do here.
      pybind11::dict r;
      result = r;
    }

    /// Initialises the container within CosmoBit from classy, but designed specifically
    /// to be used when MontePython is in use. This will ensure additional outputs are
    /// computed by classy CLASS to be passed to MontePython:
    /// When initialising the MontePython Likelihood objects they add the output that needs to be computed by class
    /// to the input dictionary. We need to get these before starting the class run
    /// e.g. for Planck_SZ likelihood the entries {'output': ' mPk ', 'P_k_max_h/Mpc': '1.'} need to be added
    /// to compute all needed observables, these entries are collected here.
    void set_classy_input_with_MPLike(pybind11::dict& result)
    {
      using namespace Pipes::set_classy_input_with_MPLike;
      static bool first = true;

      // Only do this this first time through, and if something actually needs info from MP downstream
      if (first and Downstream::neededFor("MP_LogLikes"))
      {
        // Get extra cosmo_arguments from MP (gives a dictionary with output values that need
        // to be set for the class run)
        pybind11::dict tmp_dict = std::get<0>(*Dep::MP_objects).attr("cosmo_arguments");
        // Stringify all values in the dictionary and strip off leading and trailing whitespaces
        for (auto it: tmp_dict)
        {
          std::string key = (pybind11::str(it.first)).cast<std::string>();
          std::string val = (pybind11::str(it.second)).cast<std::string>();
          boost::algorithm::trim(val);
          result[key.c_str()] = val.c_str();
        }
        logger() << LogTags::debug << "Extra cosmo_arguments needed from MP Likelihoods: ";
        logger() << pybind11::repr(result) << EOM;
        first = false;
      }
    }

    /// Set the parameters for exoCLASS for a scenario with annihilating dark matter.
    void set_classy_parameters_EnergyInjection_AnnihilatingDM(pybind11::dict &result)
    {
      using namespace Pipes::set_classy_parameters_EnergyInjection_AnnihilatingDM;

      // Make sure nothing from previous run is contained
      result.clear();

      // Set relevant inputs for the scenario of s-wave annihilating DM
      const ModelParameters& NP_params = *Dep::AnnihilatingDM_general_parameters;
      result["DM_annihilation_cross_section"] = NP_params.at("sigmav");
      result["DM_annihilation_mass"] = NP_params.at("mass");

      // Get the results from the DarkAges tables that hold extra information to be passed to the CLASS thermodynamics structure
      static DarkAges::Energy_injection_efficiency_table fz;
      static DarkAges::Energy_injection_efficiency_table cached_fz;

      fz = *Dep::energy_injection_efficiency;
      bool f_eff_mode = fz.f_eff_mode;

      // Flag passed to CLASS to signal that the energy_deposition_function is coming from GAMBIT
      // we patched exoCLASS to accept this. An alternative way without patching would be to write the tables to disk &
      // just have CLASS read in the file. To avoid the repeated file writing & deleting we pass pointers to the vector/arrays
      // to CLASS instead
      if (f_eff_mode)
      {
        result["f_eff_type"] = "pointer_to_fz_eff";
      }
      else
      {
        result["f_eff_type"] = "pointer_to_fz_channel";
      }

      // Set the lengths of the input tables (since we are passing pointers to arrays CLASS has to know how long they are)
      result["energyinj_coef_num_lines"] = fz.redshift.size();

      // Add the pointers to arrays class needs to know about to input dictionary
      // NOTE: memory addresses are passed as strings (the Python wrapper for CLASS 
      // converts every entry to a string internally so we need to do that for the 
      // memory addresses, before Python casts them to something else)
      result["energyinj_coef_z"] = memaddress_to_uint(fz.redshift.data());
      if (f_eff_mode)
      {
        result["energyinj_coef_tot"] = memaddress_to_uint(fz.f_eff.data());
      }
      else
      {
        result["energyinj_coef_heat"] = memaddress_to_uint(fz.f_heat.data());
        result["energyinj_coef_lya"] = memaddress_to_uint(fz.f_lya.data());
        result["energyinj_coef_ionH"] = memaddress_to_uint(fz.f_hion.data());
        result["energyinj_coef_ionHe"] = memaddress_to_uint(fz.f_heion.data());
        result["energyinj_coef_lowE"] = memaddress_to_uint(fz.f_lowe.data());
      }

      // Check if the table has changed compared to the previous iteration.
      // If so, notify class by adding {"EnergyInjection_changed":"yes"}
      // to the dictionary.
      // The classy frontend will just look for the key - the value is not important here.
      if (fz != cached_fz)
        result["EnergyInjection_changed"] = "yes";

      // Copy fz to cache
      cached_fz = fz;
    }

    /// Set the parameters for exoCLASS for a scenario with decaying dark matter.
    void set_classy_parameters_EnergyInjection_DecayingDM(pybind11::dict &result)
    {
      using namespace Pipes::set_classy_parameters_EnergyInjection_DecayingDM;

      // Make sure nothing from previous run is contained
      result.clear();

      // Set relevant inputs for the scenario of decaying DM
      const ModelParameters& NP_params = *Dep::DecayingDM_general_parameters;
      result["DM_decay_tau"] = NP_params.at("lifetime");
      result["DM_decay_fraction"] = NP_params.at("fraction");

      // Get the results from the DarkAges tables that hold extra information to be passed to the CLASS thermodynamics structure
      static DarkAges::Energy_injection_efficiency_table fz;
      static DarkAges::Energy_injection_efficiency_table cached_fz;

      fz = *Dep::energy_injection_efficiency;
      bool f_eff_mode = fz.f_eff_mode;

      // Flag passed to CLASS to signal that the energy_deposition_function is coming from GAMBIT
      // (exoCLASS has been patched to accept this). An alternative way without patching would be 
      // to write the tables to disk & just have CLASS read in the file. To avoid the repeated 
      // file writing & deleting we pass pointers to the vector/arrays to CLASS instead.
      if (f_eff_mode)
      {
        result["f_eff_type"] = "pointer_to_fz_eff";
      }
      else
      {
        result["f_eff_type"] = "pointer_to_fz_channel";
      }

      // Set the lengths of the input tables (since we are passing pointers to arrays CLASS has to know how long they are)
      result["energyinj_coef_num_lines"] = fz.redshift.size();

      // Add the pointers to arrays class needs to know about to input dictionary
      // NOTE: memory addresses are passed as strings (the Python wrapper for CLASS 
      // converts every entry to a string internally so we need to do that for the 
      // memory addresses, before Python casts them to something else)
      result["energyinj_coef_z"] = memaddress_to_uint(fz.redshift.data());
      if (f_eff_mode)
      {
        result["energyinj_coef_tot"] = memaddress_to_uint(fz.f_eff.data());
      }
      else
      {
        result["energyinj_coef_heat"] = memaddress_to_uint(fz.f_heat.data());
        result["energyinj_coef_lya"] = memaddress_to_uint(fz.f_lya.data());
        result["energyinj_coef_ionH"] = memaddress_to_uint(fz.f_hion.data());
        result["energyinj_coef_ionHe"] = memaddress_to_uint(fz.f_heion.data());
        result["energyinj_coef_lowE"] = memaddress_to_uint(fz.f_lowe.data());
      }

      // Check if the table has changed compared to the previous iteration.
      // If not, notify class by adding {"EnergyInjection_changed":"yes"}
      // to the dictionary.
      // The classy frontend will just look for the key - the value is not important here.
      if (fz != cached_fz)
        result["EnergyInjection_changed"] = "yes";

      // Copy fz to cache
      cached_fz = fz;
    }

    /// Add all inputs for CLASS needed to produce the correct output to be
    /// able to compute the Planck CMB likelihoods
    void set_classy_PlanckLike_input(pybind11::dict &result)
    {
      using namespace Pipes::set_classy_PlanckLike_input;

      // Make sure nothing from previous run is contained
      result.clear();

      static std::ostringstream output;
      static std::ostringstream l_max_scalars;

      static bool first = true;
      if(first)
      {
        int lmax = -1;
        bool needs_tCl = false;
        bool needs_pCl = false;

        // Get requirements of the loaded likelihoods in the plc backend
        BEreq::plc_required_Cl(lmax,needs_tCl,needs_pCl);

        // Prepare the classy input for "output"
        // -- The likelihoods need the lensed Cl such that lensing is required everytime
        output << "lCl";
        // -- Are additional Cl, other to Cl_phiphi, required?
        if (needs_tCl)
          output << ", tCl";
        if (needs_pCl)
          output << ", pCl";

        // Prepare the classy input for "l_max_scalars"
        l_max_scalars << lmax;

        first = false;
      }

      result["lensing"] = "yes";
      result["non linear"] = "halofit";
      result["output"] = output.str();
      result["l_max_scalars"] = l_max_scalars.str();
    }

    /* Classy getter functions */

    /// Hubble
    void get_H0_classy(double &result)
    {
      using namespace Pipes::get_H0_classy;

      // Rescale by c [km/s]
      result = _c_SI_*BEreq::class_get_H0()/1000;
    }

    /// Energy densities *today* (Omega0)

    /// Matter
    void get_Omega0_m_classy(double& result)
    {
      using namespace Pipes::get_Omega0_m_classy;

      result = BEreq::class_get_Omega0_m();
    }

    /// Radiation
    void get_Omega0_r_classy(double& result)
    {
      using namespace Pipes::get_Omega0_r_classy;

      result = BEreq::class_get_Omega0_r();
    }

    /// Ultra-relativistic
    void get_Omega0_ur_classy(double& result)
    {
      using namespace Pipes::get_Omega0_ur_classy;

      result = BEreq::class_get_Omega0_ur();
    }

    /// Non-cold dark matter
    void get_Omega0_ncdm_tot_classy(double& result)
    {
      using namespace Pipes::get_Omega0_ncdm_tot_classy;

      result = BEreq::class_get_Omega0_ncdm_tot();
    }

    /// Sigma8
    void get_Sigma8_classy(double& result)
    {
      using namespace Pipes::get_Sigma8_classy;

      double sigma8 = BEreq::class_get_sigma8();
      double Omega0_m = *Dep::Omega0_m;

      result = sigma8*pow(Omega0_m/0.3, 0.5);
    }

    /// Effective number of neutrino species
    /// (mostly for cross-checking!)
    void get_Neff_classy(double& result)
    {
      using namespace Pipes::get_Neff_classy;

      result = BEreq::class_get_Neff();
    }

    /// Comoving sound horizon at baryon drag epoch
    void get_rs_drag_classy(double& result)
    {
      using namespace Pipes::get_rs_drag_classy;

      result = BEreq::class_get_rs();
    }

    /// -----------
    /// Background Likelihoods (H0 + Sigma8)
    /// -----------

    /// Compute a Gaussian lnL for the Hubble constant from a given DataFile.
    /// The file must contain an ASCII table with a single row entry, 
    /// containing the mean value and (1 sigma) error.
    void compute_H0_LogLike(double &result)
    {
      using namespace Pipes::compute_H0_LogLike;

      const str filename = runOptions->getValue<std::string>("DataFile");
      const str path_to_file = GAMBIT_DIR "/CosmoBit/data/H0/" + runOptions->getValue<std::string>("DataFile");
      static ASCIItableReader data(path_to_file);
      static bool read_data = false;

      if(read_data == false)
      {
        logger() << "H0 data read from file '"<<filename<<"'." << EOM;
        std::vector<std::string> colnames = initVector<std::string>("mean", "sigma");
        data.setcolnames(colnames);

        if(data.getnrow() != 1)
        {
          std::ostringstream err;
          err << data.getnrow() << " data points for H0 measurement data read from '"<< filename << "'.\n Only one expected.";
          CosmoBit_error().raise(LOCAL_INFO, err.str());
        }
        read_data = true;
      }
      result = -0.5 * pow(*Dep::H0 - data["mean"][0],2)/ pow(data["sigma"][0],2);
    }

    /// Compute a Gaussian lnL for sigma8 from a given DataFile.
    /// The quantity constrained is: sigma8*(Omega_m/Omega_m_ref)^Omega_m_index (= bestfit +/- sigma)
    /// The file must contain an ASCII table with entries:s 
    /// 1. Omega_m_ref 2. Omega_m_index 3. bestfit 4. sigma 
    void compute_sigma8_LogLike(double &result)
    {
      using namespace Pipes::compute_sigma8_LogLike;

      double sigma8,Omega0_m; // cosmo distances and sound horizon at drag,Hubble, theoretical prediction
      double chi2 =0;

      const str filename = runOptions->getValue<std::string>("DataFile");
      const str path_to_file = GAMBIT_DIR "/CosmoBit/data/sigma8/" + runOptions->getValue<std::string>("DataFile");
      static ASCIItableReader data(path_to_file);
      static bool read_data = false;
      static int nrow;

      sigma8 = BEreq::class_get_sigma8();
      Omega0_m = *Dep::Omega0_m;

      if(read_data == false)
      {
        logger() << "Sigma8 data read from file '"<<filename<<"'." << EOM;
        std::vector<std::string> colnames = initVector<std::string>("Omega_m_ref","Omega_m_index", "bestfit","sigma");
        data.setcolnames(colnames);
        nrow=data.getnrow();

        read_data = true;
      }

      for(int ie=0;ie<nrow;ie++)
      {
        chi2 += pow(( sigma8*pow(Omega0_m/data["Omega_m_ref"][ie], data["Omega_m_index"][ie]) - data["bestfit"][ie])/ data["sigma"][ie],2);
      }
      result = -0.5*chi2;
    }

    /***************/
    /* MontePython */
    /***************/

    /// Function to fill the mcmc_parameters dictionary of MontePython's Data object
    /// with current values of nuisance parameters.
    void set_parameter_dict_for_MPLike(pybind11::dict & result)
    {
      using namespace Pipes::set_parameter_dict_for_MPLike;
      using namespace pybind11::literals;

      // The loop has to be executed for every parameter point. It takes about 0.00023s -> ~4 minutes for 1e6 points
      for (auto it=Param.begin(); it != Param.end(); it++)
      {
        std::string name = it->first;
        double value = *Param[name];

        // Check if any models are scanned for which we had to rename the nuisance parameters due to:
        //    a) parameters having the same name  -> e.g. 'epsilon' and 'sigma_NL'
        //    b) parameter names containing symbols that can't be used in macros -> e.g. "^" in 'beta_0^Euclid'

        // a) have to rename parameters epsilon_ska, epsilon_euclid,.. to "epsilon" as they are implemented in MontePython
        if (name.find("epsilon") != std::string::npos){name="epsilon";}

        // a) have to rename parameters sigma_NL_ska, sigma_NL_euclid,.. to "sigma_NL" as they are implemented in MontePython
        else if (name.find("sigma_NL") != std::string::npos){name="sigma_NL";}

        // b) get the "^" characters back into the parameter names
        //   -> beta_x<experiment> has to be beta_x^<experiment> where x = 0 or 1 and <experiment> = Euclid, SKA1 or SKA2
        //if(ModelInUse("cosmo_nuisance_euclid_pk") or ModelInUse("cosmo_nuisance_ska"))
        else if (name.find("beta_")!= std::string::npos){name = name.insert(6,"^");}

        result[name.c_str()] = pybind11::dict("current"_a=value,"scale"_a=1.); // scale always 1 in GAMBIT
      }
    }

    /// Function to fill the mcmc_parameters dictionary of MontePython's Data object with an empty dictionary.
    /// This version of the capability 'parameter_dict_for_MPLike' is used when no Likelihood with nuisance
    /// parameters are in use, and just passes an empty Python dictionary
    void pass_empty_parameter_dict_for_MPLike(pybind11::dict & result)
    {
      using namespace Pipes::pass_empty_parameter_dict_for_MPLike;
      pybind11::dict r;
      result = r;
      // Nothing to do here.
    }

    /// Create the MontePython data and likelihood objects, determining which experiments are in use in the process
    void create_MP_objects(MPLike_objects_container &result)
    {
      using namespace Pipes::create_MP_objects;
      static map_str_pyobj likelihoods;
      static bool first = true;

      // Determine which likelihoods to compute and initialise the relevant MontePython objects.
      if (first)
      {
        // Set up some references for easier reading
        auto& data = std::get<0>(result);
        auto& experiments = std::get<1>(result);
        auto& likelihoods = std::get<2>(result);

        // Get list of likelihoods implemented in MP
        std::vector<str> avail_likes = BEreq::get_MP_available_likelihoods();

        // Get the list of the experiments and datafiles given in the YAML file as sub-capabilities.
        // Using the default datafile can be achieved by leaving the datafile out, or setting it to "default".
        YAML::Node subcaps = Downstream::subcaps->getNode();
        std::vector<YAML::Node> empties;
        for (const auto& x : subcaps) if (x.second.IsNull()) empties.push_back(x.first);
        for (const auto& x : empties) subcaps[x] = "default";
        if (subcaps.IsNull())
        {
          if (Downstream::neededFor("MP_LogLikes"))
          {
            std::ostringstream ss;
            ss << "No sub-capabilities found when attempting to create MontePython objects." << endl
               << "This can happen because you either forgot to choose any experiments," << endl
               << "or because you used incorrect syntax to choose them as sub-capabilities." << endl
               << "You can do this in the relevant entry of the ObsLikes section of your YAML file," << endl
               << "by setting sub_capabilities as a scalar (if you only want one experiment), e.g." << endl
               << "    sub_capabilities: bao_smallz_2014" << endl
               << "or as a sequence (if you don't need to specify data files), e.g." << endl
               << "    sub_capabilities:" << endl
               << "      - bao_smallz_2014" << endl
               << "      - Pantheon" << endl
               << "or even as a map (if you want to specify data files), e.g." << endl
               << "    sub_capabilities:" << endl
               << "      bao_smallz_2014: default" << endl
               << "      Pantheon: default" << endl;
            CosmoBit_error().raise(LOCAL_INFO, ss.str());
          }
        }
        else experiments = subcaps.as<map_str_str>();

        // Check that all the requested likelihoods can actually be provided by MP
        for (const auto& x : experiments)
        {
          if (std::find(avail_likes.begin(), avail_likes.end(), x.first) == avail_likes.end())
          {
            str errmsg = "Likelihood '" + x.first + "' is not implemented in MontePython. Check for typos or implement it.\nLikelihoods currently available are:\n";
            for (const auto& value : avail_likes) errmsg += ("\t"+value+"\n");
            CosmoBit_error().raise(LOCAL_INFO, errmsg);
          }
          logger() << LogTags::debug << "Read MontePythonLike option "<< x.first << ", using data file " << x.second << EOM;
        }

        // MPLike_data_container should only be created and set once, when calculating the first point.
        // After that it has to be kept alive since it contains a vector with the initialised MPLike Likelihood objects.
        data = BEreq::create_MP_data_object(experiments);

        // Add current parameters to data object to enable check if all nuisance parameters are
        // scanned upon initialisation of likelihood objects
        data.attr("mcmc_parameters") = *Dep::parameter_dict_for_MPLike;
        likelihoods = BEreq::create_MP_likelihood_objects(data, experiments);

        // It's been nice, but let's not do this again.
        first = false;
      }
    }


    /// Computes lnL for each experiment initialised in MontePython
    void compute_MP_LogLikes(map_str_dbl & result)
    {
      using namespace Pipes::compute_MP_LogLikes;
      static pybind11::object data = std::get<0>(*Dep::MP_objects);
      static const map_str_str& experiments = std::get<1>(*Dep::MP_objects);
      static const map_str_pyobj& likelihoods = std::get<2>(*Dep::MP_objects);
      static const MPLike_data_container mplike_cont(data, likelihoods);

      // Pass current values of nuisance parameters to data.mcmc_parameters dictionary for likelihood computation in MP
      mplike_cont.data.attr("mcmc_parameters") = *Dep::parameter_dict_for_MPLike;

      // Create instance of classy class Class
      pybind11::object cosmo = BEreq::get_classy_cosmo_object();

      // Loop through the list of experiments, and query the lnL from the MontePython backend.
      for (sspair it : experiments)
      {
        // Likelihood names are keys of experiment map (str, str map mapping likelihood name to .data file)
        double logLike = BEreq::get_MP_loglike(mplike_cont, cosmo, it.first);
        result[it.first] = logLike;
        logger() << "(compute_MP_LogLikes):  name: " << it.first << "\tvalue: " << logLike << EOM;
      }
    }

    /// Computes the combined lnL from the set of experiments
    /// given to MontePython.
    void compute_MP_combined_LogLike(double& result)
    {
      using namespace Pipes::compute_MP_combined_LogLike;

      // Get likelihoods computed by MontePython
      map_str_dbl MP_lnLs = *Dep::MP_LogLikes;

      // Retrieve the sub-capabilities requested in the YAML file
      std::vector<str> subcaps = Downstream::subcaps->getNames();

      // Iterate through map of doubles and return one big fat double,
      // selecting only those entries specified as sub-capabilities.
      double lnL = 0.;
      logger() << LogTags::debug << "(compute_MP_combined_LogLike):";
      for (const auto &p : MP_lnLs)
      {
        if (std::find(subcaps.begin(), subcaps.end(), p.first) != subcaps.end())
        {
          logger() << endl << "  name: "  << p.first << "\tvalue: " << p.second;
          lnL += p.second;
        }
      }
      logger() << EOM;
      result = lnL;
    }

    /// Get correlation coefficients and uncorrelated likelihood
    /// of MP likelihood "bao_correlations".
    /// Warning: this routine is specific to this likelihood, don't use for anything else!
    void get_bao_like_correlation(map_str_dbl& result)
    {
      using namespace Pipes::get_bao_like_correlation;

      // This function has a dependency on MP_LogLikes even though it is not directly
      // needed in the calculation. However, through this dependency we make sure that
      // MP was called before this function is executed -> don't remove it!

      // Get map containing python likelihood objects
      static const map_str_pyobj& likelihoods = std::get<2>(*Dep::MP_objects);

      // Check if "bao_correlations" likelihood was computed, if so
      // retrieve correlation coefficients and uncorrelated likelihood value
      if(likelihoods.find("bao_correlations") != likelihoods.end())
      {
          result["uncorrelated_loglike"] = likelihoods.at("bao_correlations").attr("uncorrelated_loglike").cast<double>();
          pybind11::list corr_coeffs =  likelihoods.at("bao_correlations").attr("correlation_coeffs");
          result["correlation_coeffs_0"] = corr_coeffs[0].cast<double>();
          result["correlation_coeffs_1"] = corr_coeffs[1].cast<double>();
          result["correlation_coeffs_2"] = corr_coeffs[2].cast<double>();
      }
      else
      {
          str errmsg = "Likelihood 'bao_correlations' was not requested in the YAML file, but you are asking for\n";
          errmsg += "the correlation coefficients from this likelihood. Either remove 'bao_like_correlation' from the ObsLikes section\n";
          errmsg += "in your YAML file or include the computation of the 'bao_correlations' likelihood by adding:\n\n";
          errmsg += "  - purpose:      LogLike\n";
          errmsg += "    capability:   MP_Combined_LogLike\n";
          errmsg += "    sub_capabilities:\n";
          errmsg += "      - bao_correlations\n\n";
          errmsg += "to the YAML file.";
          CosmoBit_error().raise(LOCAL_INFO, errmsg);
      }
    }

  } // namespace CosmoBit

} // namespace Gambit