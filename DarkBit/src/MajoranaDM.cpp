//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Implementation of MajoranaDM routines.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ankit Beniwal
///          (ankit.beniwal@adelaide.edu.au)
///  \date Oct 2016
///  \date Jun, Sep 2017
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

    void DarkMatter_ID_MajoranaDM(std::string & result) { result = "X"; }

    /// Direct detection couplings for the MajoranaDM model.
    void DD_couplings_MajoranaDM(DM_nucleon_couplings &result)
    {
      using namespace Pipes::DD_couplings_MajoranaDM;
      const Spectrum& spec = *Dep::MajoranaDM_spectrum;
      const SubSpectrum& he = spec.get_HE();
      //double mass = spec.get(Par::Pole_Mass,"X");
      double lambda = he.get(Par::dimensionless,"lX");
      double cosXI = he.get(Par::dimensionless,"cosXI");
      double sinXI = sqrt(1-pow(cosXI,2));
      double mh = spec.get(Par::Pole_Mass,"h0_1");

      // Expressions taken from Cline et al. (2013, PRD 88:055025, arXiv:1306.4710)
      double fp = 2./9. + 7./9.*(*Param["fpu"] + *Param["fpd"] + *Param["fps"]);
      double fn = 2./9. + 7./9.*(*Param["fnu"] + *Param["fnd"] + *Param["fns"]);

      // SI scalar and pseudoscalar couplings
      result.gps = lambda*fp*m_proton*cosXI/pow(mh,2);
      result.gns = lambda*fn*m_neutron*cosXI/pow(mh,2);
      result.gpa = lambda*fp*m_proton*sinXI/pow(mh,2);
      result.gna = lambda*fn*m_neutron*sinXI/pow(mh,2);

      logger() << LogTags::debug << "Majorana DM DD couplings:" << std::endl;
      logger() << " gps = " << result.gps << std::endl;
      logger() << " gns = " << result.gns << std::endl;
      logger() << " gpa = " << result.gpa << std::endl;
      logger() << " gna = " << result.gna << EOM;

    } // function DD_couplings_MajoranaDM

    /// Set up process catalog for the MajoranaDM model.
    void TH_ProcessCatalog_MajoranaDM(DarkBit::TH_ProcessCatalog &result)
    {
      using namespace Pipes::TH_ProcessCatalog_MajoranaDM;
      using std::vector;
      using std::string;

      std::string DMid = *Dep::DarkMatter_ID;
      if ( DMid != "X" )
      {
        invalid_point().raise("TH_ProcessCatalog_MajoranaDM requires DMid to be X.");
      }

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
      const Spectrum& spec = *Dep::MajoranaDM_spectrum;
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
      double mX = spec.get(Par::Pole_Mass,"X");
      double mH = spec.get(Par::Pole_Mass,"h0_1");
      addParticle("X",        mX, 0)  // Majorana DM
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
      // Import two-body annihilation processes from micrOmegas_3.6.9.2
      ////////////////////////////////////////////////////////////////////

      // Set of possible final state particles
      std::set<string> annFinalStates;

      // Initialize main annihilation process
      TH_Process process_ann(DMid, DMid);

      // Helper variables
      int err, key = 1;
      double *SpA = NULL, *SpE = NULL, *SpP = NULL;
      double *SpNe = NULL, *SpNm = NULL, *SpNl = NULL;
      double m_1, m_2, sigmav_total, min_prop = 1e-6;

      // Calculate the total sigmav using the calcSpectrum function in micrOmegas_3.6.9.2
      sigmav_total = BEreq::calcSpectrum(byVal(key), byVal(SpA), byVal(SpE),
                    byVal(SpP), byVal(SpNe), byVal(SpNm), byVal(SpNl), &err);
      logger() << "Total zero-velocity annihilation cross section = " << sigmav_total << " cm^3/s" << std::endl;

      // Convenience macros for setting up 2-body annihilations using micrOmega's functions
      #define SETUP_KINEMATIC_PROCESS_MO(NAME, MO_PRTCL_NAME, P1, P2)                                   \
      m_1 = catalog.getParticleProperty(STRINGIFY(P1)).mass;                                            \
      m_2 = catalog.getParticleProperty(STRINGIFY(P2)).mass;                                            \
      if(2*mX > m_1 + m_2)                                                                              \
      {                                                                                                 \
      for (int i = 0; ((*BEreq::vSigmaCh)+i)->weight > min_prop; i++)                                   \
        {                                                                                               \
          /* Calculate sigmav for the input channel */                                                  \
          if (strcmp(((*BEreq::vSigmaCh)+i)->prtcl[2],STRINGIFY(MO_PRTCL_NAME)) == 0)                   \
          {                                                                                             \
            double CAT(sigma_,NAME) = (((*BEreq::vSigmaCh)+i)->weight)*sigmav_total;                    \
            logger() << "  BR(" << DMid << " + " << DMid << " -> "                                      \
                     << ((*BEreq::vSigmaCh)+i)->prtcl[2]                                                \
                     << " + " << ((*BEreq::vSigmaCh)+i)->prtcl[3] << ") = "                             \
                     << ((*BEreq::vSigmaCh)+i)->weight << std::endl;                                    \
            /* Create associated kinematical functions */                                               \
            daFunk::Funk CAT(kinematicFunction_,NAME) = daFunk::cnst(CAT(sigma_,NAME), "v");            \
            /* Create channel identifier string */                                                      \
            std::vector<std::string> CAT(finalStates_,NAME);                                            \
            CAT(finalStates_,NAME).push_back(STRINGIFY(P1));                                            \
            CAT(finalStates_,NAME).push_back(STRINGIFY(P2));                                            \
            /* Create channel and push it into channel list of process */                               \
            TH_Channel CAT(channel_,NAME)(CAT(finalStates_,NAME),CAT(kinematicFunction_,NAME));         \
            process_ann.channelList.push_back(CAT(channel_,NAME));                                      \
            annFinalStates.insert(STRINGIFY(P1));                                                       \
            annFinalStates.insert(STRINGIFY(P2));                                                       \
          }                                                                                             \
        }                                                                                               \
      }

      // Include SM final states from DM annihilation
      SETUP_KINEMATIC_PROCESS_MO(WW, W+, W+, W-);
      SETUP_KINEMATIC_PROCESS_MO(ZZ, Z, Z0, Z0);
      SETUP_KINEMATIC_PROCESS_MO(HH, H, h0_1, h0_1);
      SETUP_KINEMATIC_PROCESS_MO(AA, A, gamma, gamma);

      // down-type quarks
      SETUP_KINEMATIC_PROCESS_MO(dD, d, d_1, dbar_1);
      SETUP_KINEMATIC_PROCESS_MO(sS, s, d_2, dbar_2);
      SETUP_KINEMATIC_PROCESS_MO(bB, b, d_3, dbar_3);

      // up-type quarks
      SETUP_KINEMATIC_PROCESS_MO(uU, u, u_1, ubar_1);
      SETUP_KINEMATIC_PROCESS_MO(cC, c, u_2, ubar_2);
      SETUP_KINEMATIC_PROCESS_MO(tT, t, u_3, ubar_3);

      // leptons
      SETUP_KINEMATIC_PROCESS_MO(eE, e, e+_1, e-_1);
      SETUP_KINEMATIC_PROCESS_MO(mM, m, e+_2, e-_2);
      SETUP_KINEMATIC_PROCESS_MO(tauTau, l, e+_3, e-_3);

      // Get rid of convenience macro
      #undef SETUP_KINEMATIC_PROCESS_MO

      // Add process to previous list
      catalog.processList.push_back(process_ann);

      // Validate
      catalog.validate();

      // Return the finished process catalog
      result = catalog;

    } // function TH_ProcessCatalog_MajoranaDM
  }
}
