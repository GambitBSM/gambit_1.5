//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Implementation of VectorDM routines 
///  (file format is based on SingletDM.cpp and 
///   MSSM.cpp).
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ankit Beniwal
///          (ankit.beniwal@adelaide.edu.au)
///  \date Aug 2016
///
///  *********************************************

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/Elements/virtual_higgs.hpp"
#include "gambit/DarkBit/DarkBit_rollcall.hpp"
#include "gambit/Utils/ascii_table_reader.hpp"
#include "boost/make_shared.hpp"
#include "gambit/DarkBit/DarkBit_utils.hpp"

namespace Gambit
{

  namespace DarkBit
  {

    void DarkMatter_ID_VectorDM(std::string & result) { result = "V"; }

    /// Direct detection couplings for the VectorDM model.
    void DD_couplings_VectorDM(DM_nucleon_couplings &result)
    {
      using namespace Pipes::DD_couplings_VectorDM;
      const Spectrum& spec = *Dep::VectorDM_spectrum;
      const SubSpectrum& he = spec.get_HE();
      double mass = spec.get(Par::Pole_Mass,"V");
      double lambda = he.get(Par::dimensionless,"lambda_hV");    
      double mh = spec.get(Par::Pole_Mass,"h0_1");

      // Expressions taken from Cline et al. (2013, PRD 88:055025, arXiv:1306.4710)
      double fp = 2./9. + 7./9.*(*Param["fpu"] + *Param["fpd"] + *Param["fps"]);
      double fn = 2./9. + 7./9.*(*Param["fnu"] + *Param["fnd"] + *Param["fns"]);

      result.gps = lambda*fp*m_proton/pow(mh,2)/mass/2;
      result.gns = lambda*fn*m_neutron/pow(mh,2)/mass/2;
      result.gpa = 0;  // Only SI cross-section
      result.gna = 0;

      logger() << LogTags::debug << "Vector DM DD couplings:" << std::endl;
      logger() << " gps = " << result.gps << std::endl;
      logger() << " gns = " << result.gns << std::endl;
      logger() << " gpa = " << result.gpa << std::endl;
      logger() << " gna = " << result.gna << EOM;

    } // function DD_couplings_VectorDM

    std::map<std::string, daFunk::Funk> get_f_vs_mass(std::string filename)
    {
      // Higgs branching ratios and total width Gamma [GeV], as function of
      // mass [GeV] (90 - 150 GeV)
      ASCIItableReader table(filename);
      std::vector<std::string> colnames =
        initVector<std::string>("mass", "bb", "tautau", "mumu",
            "ss", "cc", "tt", "gg", "gammagamma", "Zgamma",
            "WW", "ZZ", "Gamma");
      table.setcolnames(colnames);

      std::map<std::string, daFunk::Funk> f_vs_mass;
      for (auto it = colnames.begin(); it != colnames.end(); it++)
      {
        f_vs_mass[*it] = daFunk::interp("mass", table["mass"], table[*it]);
      }
      return f_vs_mass;
    }

    /// Set up process catalog for the VectorDM model.
    /// Uses micrOmega's functions to compute sigma-v and branching ratios.
    void TH_ProcessCatalog_VectorDM(DarkBit::TH_ProcessCatalog &result)
    {
      using namespace Pipes::TH_ProcessCatalog_VectorDM;
      using std::vector;
      using std::string;

      std::string DMid = *Dep::DarkMatter_ID;
      if ( DMid != "V" )
      {
        invalid_point().raise("TH_ProcessCatalog_VectorDM requires DMid to be V.");
      }

      // Initialize Higgs decay tables (static, hence only once)
      static std::map<string, daFunk::Funk> f_vs_mass =
        get_f_vs_mass("Elements/data/Higgs_decay_1101.0593.dat");

      // Initialize empty catalog 
      TH_ProcessCatalog catalog;

      ///////////////////////////////////////
      // Import particle masses and couplings
      ///////////////////////////////////////

      // Convenience macros
      #define getSMmass(Name, spinX2)                                           \
       catalog.particleProperties.insert(std::pair<string, TH_ParticleProperty> \
       (Name , TH_ParticleProperty(SM.get(Par::Pole_Mass,Name), spinX2)));    
      #define addParticle(Name, Mass, spinX2)                                   \
       catalog.particleProperties.insert(std::pair<string, TH_ParticleProperty> \
       (Name , TH_ParticleProperty(Mass, spinX2)));

      // Import Spectrum objects
      const Spectrum& spec = *Dep::VectorDM_spectrum;
      const SubSpectrum& he = spec.get_HE();
      const SubSpectrum& SM = spec.get_LE();
      const SMInputs& SMI   = spec.get_SMInputs();
  
      // Get SM pole masses
      getSMmass("e-_1",     1)
      getSMmass("e+_1",     1)
      getSMmass("e-_2",     1)
      getSMmass("e+_2",     1)
      getSMmass("e-_3",     1)
      getSMmass("e+_3",     1)
      getSMmass("Z0",     2)
      getSMmass("W+",     2)
      getSMmass("W-",     2)
      getSMmass("g",      2)
      getSMmass("gamma",  2)
      getSMmass("u_3",      1)
      getSMmass("ubar_3",   1)
      getSMmass("d_3",      1)
      getSMmass("dbar_3",   1)

      // Pole masses not available for the light quarks.
      addParticle("u_1"   , SMI.mU,  1) // mu(2 GeV)^MS-bar, not pole mass
      addParticle("ubar_1", SMI.mU,  1) // mu(2 GeV)^MS-bar, not pole mass
      addParticle("d_1"   , SMI.mD,  1) // md(2 GeV)^MS-bar, not pole mass
      addParticle("dbar_1", SMI.mD,  1) // md(2 GeV)^MS-bar, not pole mass
      addParticle("u_2"   , SMI.mCmC,1) // mc(mc)^MS-bar, not pole mass
      addParticle("ubar_2", SMI.mCmC,1) // mc(mc)^MS-bar, not pole mass
      addParticle("d_2"   , SMI.mS,  1) // ms(2 GeV)^MS-bar, not pole mass
      addParticle("dbar_2", SMI.mS,  1) // ms(2 GeV)^MS-bar, not pole mass

      // Masses for neutrino flavour eigenstates. Set to zero.
      // (presently not required)
      addParticle("nu_e",     0.0, 1)
      addParticle("nubar_e",  0.0, 1)
      addParticle("nu_mu",    0.0, 1)
      addParticle("nubar_mu", 0.0, 1)
      addParticle("nu_tau",   0.0, 1)
      addParticle("nubar_tau",0.0, 1)

      // Higgs-sector masses
      double mV = spec.get(Par::Pole_Mass,"V");
      double mH = spec.get(Par::Pole_Mass,"h0_1");
      addParticle("V",        mV, 0)  // Vector DM
      addParticle("h0_1",     mH, 0)  // SM-like Higgs
      addParticle("pi0",   meson_masses.pi0,       0)
      addParticle("pi+",   meson_masses.pi_plus,   0)
      addParticle("pi-",   meson_masses.pi_minus,  0)
      addParticle("eta",   meson_masses.eta,       0)
      addParticle("rho0",  meson_masses.rho0,      1)
      addParticle("rho+",  meson_masses.rho_plus,  1)
      addParticle("rho-",  meson_masses.rho_minus, 1)
      addParticle("omega", meson_masses.omega,     1)

      // Get rid of convenience macros
      #undef getSMmass
      #undef addParticle

      ////////////////////////////////////////////////////////////////////
      // Import two-body annihilation processes from MicrOmegas v3.6.9.2
      ////////////////////////////////////////////////////////////////////

      // Set of possible final state particles. Used to determine which decays to import.
      std::set<string> annFinalStates;

      // Initialize main annihilation process
      TH_Process process_ann(DMid, DMid);
      
      // Helper variables
      int err, key = 4, NZ = 250;
      double SpA[NZ], SpE[NZ], SpP[NZ];
      double *SpNe = NULL, *SpNm = NULL, *SpNl = NULL;
      double m_1, m_2, sigmav_total, min_prop = 1e-6;
         
      // Calculate the total sigmav using the calcSpectrum function from MicrOmegas
      sigmav_total = BEreq::calcSpectrum(byVal(key), byVal(SpA), byVal(SpE),
                    byVal(SpP), byVal(SpNe), byVal(SpNm), byVal(SpNl), &err);
      logger() << "Total sigmav = " << sigmav_total << " cm^3/s" << std::endl;

      // A convenience macro for setting up 2-body annihilations using MicrOmega's functions
      #define SETUP_MO_PROCESS(NAME, MO_PRTCL_NAME, P1, P2)                                     \
      /* Check if process is kinematically allowed */                                           \
      m_1 = catalog.getParticleProperty(STRINGIFY(P1)).mass;                                    \
      m_2 = catalog.getParticleProperty(STRINGIFY(P2)).mass;                                    \
      if(m_1 + m_2 < 2*mV)                                                                      \
      {                                                                                         \
      for (int i = 0; ((*BEreq::vSigmaCh)+i)->weight > min_prop; i++)                         \
        {                                                                                       \
          /* Calculate sigmav for the input channel */                                          \
          if (strcmp(((*BEreq::vSigmaCh)+i)->prtcl[2],STRINGIFY(MO_PRTCL_NAME)) == 0)         \
          {                                                                                     \
            double CAT(sigma_,NAME) = (((*BEreq::vSigmaCh)+i)->weight)*sigmav_total;          \
            logger() << DMid << "," << DMid << " -> " << ((*BEreq::vSigmaCh)+i)->prtcl[2]     \
               << " " << ((*BEreq::vSigmaCh)+i)->prtcl[3] << "\t"                             \
               << ((*BEreq::vSigmaCh)+i)->weight << std::endl;                                \
            /* Create associated kinematical functions */                                       \
            daFunk::Funk CAT(kinematicFunction_,NAME) = daFunk::cnst(CAT(sigma_,NAME), "v");    \
            /* Create channel identifier string */                                              \
            std::vector<std::string> CAT(finalStates_,NAME);                                    \
            CAT(finalStates_,NAME).push_back(STRINGIFY(P1));                                    \
            CAT(finalStates_,NAME).push_back(STRINGIFY(P2));                                    \
            /* Create channel and push it into channel list of process */                       \
            TH_Channel CAT(channel_,NAME)(CAT(finalStates_,NAME),CAT(kinematicFunction_,NAME)); \
            process_ann.channelList.push_back(CAT(channel_,NAME));                              \
            annFinalStates.insert(STRINGIFY(P1));                                               \
            annFinalStates.insert(STRINGIFY(P2));                                               \
          }                                                                                     \
        }                                                                                       \
      }

      // Include SM final states from DM annihilation
      SETUP_MO_PROCESS(WW, W+, W+, W-);
      SETUP_MO_PROCESS(ZZ, Z, Z0, Z0);
      SETUP_MO_PROCESS(HH, H, h0_1, h0_1);
      SETUP_MO_PROCESS(AA, A, gamma, gamma);

      // down-type quarks
      SETUP_MO_PROCESS(dD, d, d_1, dbar_1);
      SETUP_MO_PROCESS(sS, s, d_2, dbar_2);
      SETUP_MO_PROCESS(bB, b, d_3, dbar_3);

      // up-type quarks
      SETUP_MO_PROCESS(uU, u, u_1, ubar_1);
      SETUP_MO_PROCESS(cC, c, u_2, ubar_2);
      SETUP_MO_PROCESS(tT, t, u_3, ubar_3);

      // leptons
      SETUP_MO_PROCESS(eE, e, e+_1, e-_1);
      SETUP_MO_PROCESS(mM, m, e+_2, e-_2);
      SETUP_MO_PROCESS(tauTau, l, e+_3, e-_3);

      // Get rid of convenience macro
      #undef SETUP_MO_PROCESS

      /////////////////////////////
      // Import Decay information
      /////////////////////////////

      // Import decay table from DecayBit
      const DecayTable* tbl = &(*Dep::decay_rates);

      // Save Higgs width for later
      double gammaH = tbl->at("h0_1").width_in_GeV;

      // Set of imported decays
      std::set<string> importedDecays;

      // Minimum branching ratio to include
      double minBranching = 0;

      // Import relevant decays (only Higgs and subsequent decays)
      using DarkBit_utils::ImportDecays;
      ImportDecays("h0_1", catalog, importedDecays, tbl, minBranching,
          daFunk::vec<std::string>("Z0", "W+", "W-", "e+_2", "e-_2", "e+_3", "e-_3"));

      // Add process to previous list
      catalog.processList.push_back(process_ann);

      // Validate
      catalog.validate();

      // Return the finished process catalog
      result = catalog;

    } // function TH_ProcessCatalog_VectorDM
  }
}
