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
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date Nov 2017
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
  
    class MajoranaDM
    {
      public:
        /// Initialize SingletDM object (branching ratios etc)
        MajoranaDM(
            TH_ProcessCatalog* const catalog,
            double gammaH,
            double vev,
            double alpha_strong)
        : Gamma_mh(gammaH), v0 (vev),
          alpha_s (alpha_strong)
        {
          mh   = catalog->getParticleProperty("h0_1").mass;
          mb   = catalog->getParticleProperty("d_3").mass;
          mc   = catalog->getParticleProperty("u_2").mass;
          mtau = catalog->getParticleProperty("e-_3").mass;
          mt   = catalog->getParticleProperty("u_3").mass;
          mZ0  = catalog->getParticleProperty("Z0").mass;
          mW   = catalog->getParticleProperty("W+").mass;
        };
        ~MajoranaDM() {}

        /// Helper function (Breit-Wigner)
        double Dh2 (double s)
        {
          return 1/((s-mh*mh)*(s-mh*mh)+mh*mh*Gamma_mh*Gamma_mh);
        }

        /*! \brief Returns <sigma v> in cm3/s for given channel, velocity and
         *         model parameters.
         *
         * channel: bb, tautau, mumu, ss, cc, tt, gg, gammagamma, Zgamma, WW,
         * ZZ, hh
         */
        double sv(std::string channel, double lambda, double mass, double v)
        {
          // Note: Valid for mass > 45 GeV
          double s = 4*mass*mass/(1-v*v/4);
          double sqrt_s = sqrt(s);
          if ( sqrt_s < 90 )
          {
            piped_invalid_point.request(
                "MajoranaDM sigmav called with sqrt_s < 90 GeV.");
            return 0;
          }

          if ( channel == "hh" )
          {
            if ( sqrt_s > mh*2 )
            {
              double GeV2tocm3s1 = gev2cm2*s2cm;
              return sv_hh(lambda, mass, v)*GeV2tocm3s1;
            }
            else return 0;
          }

          if ( channel == "bb" and sqrt_s < mb*2 ) return 0;
          if ( channel == "cc" and sqrt_s < mc*2  ) return 0;
          if ( channel == "tautau" and sqrt_s < mtau*2 ) return 0;
          if ( channel == "tt" and sqrt_s < mt*2 ) return 0;
          if ( channel == "ZZ" and sqrt_s < mZ0*2) return 0;
          if ( channel == "WW" and sqrt_s < mW*2) return 0;

          if ( sqrt_s < 300 )
          {
            double br = virtual_SMHiggs_widths(channel,sqrt_s);
            double Gamma_s = virtual_SMHiggs_widths("Gamma",sqrt_s);
            double GeV2tocm3s1 = gev2cm2*s2cm;

            // Explicitly close channel for off-shell top quarks
            if ( channel == "tt" and sqrt_s < mt*2) return 0;

            double res = 2*lambda*lambda*v0*v0/
              sqrt_s*Dh2(s)*Gamma_s*GeV2tocm3s1*br;
            return res;
          }
          else
          {
            if ( channel == "bb" ) return sv_ff(lambda, mass, v, mb, true);
            if ( channel == "cc" ) return sv_ff(lambda, mass, v, mc, false);
            if ( channel == "tautau" ) return sv_ff(lambda, mass, v, mtau, false);
            if ( channel == "tt" ) return sv_ff(lambda, mass, v, mt, false);
            if ( channel == "ZZ" ) return sv_ZZ(lambda, mass, v);
            if ( channel == "WW" ) return sv_WW(lambda, mass, v);
          }
          return 0;
        }

        // Annihilation into W bosons.
        double sv_WW(double lambda, double mass, double v)
        {
          double s = 4*mass*mass/(1-v*v/4);
          double x = pow(mW,2)/s;
          double GeV2tocm3s1 = gev2cm2*s2cm;
          return pow(lambda,2)*s/8/M_PI*sqrt(1-4*x)*Dh2(s)*(1-4*x+12*pow(x,2))
            *GeV2tocm3s1;
        }

        // Annihilation into Z bosons.
        double sv_ZZ(double lambda, double mass, double v)
        {
          double s = 4*mass*mass/(1-v*v/4);
          double x = pow(mZ0,2)/s;
          double GeV2tocm3s1 = gev2cm2*s2cm;
          return pow(lambda,2)*s/16/M_PI*sqrt(1-4*x)*Dh2(s)*(1-4*x+12*pow(x,2))
            * GeV2tocm3s1;
        }

        // Annihilation into fermions
        double sv_ff(
            double lambda, double mass, double v, double mf, bool is_quark)
        {
          double s = 4*mass*mass/(1-v*v/4);
          double vf = sqrt(1-4*pow(mf,2)/s);
          double Xf = 1;
          if ( is_quark ) Xf = 3 *
            (1+(3/2*log(pow(mf,2)/s)+9/4)*4*alpha_s/3/M_PI);
          double GeV2tocm3s1 = gev2cm2*s2cm;
          return pow(lambda,2)*
            pow(mf,2)/4/M_PI*Xf*pow(vf,3) * Dh2(s) * GeV2tocm3s1;
        }

        /// Annihilation into hh
        double sv_hh(double lambda, double mass, double v)
        {
          double s = 4*mass*mass/(1-v*v/4);  // v is relative velocity
          double vh = sqrt(1-4*mh*mh/s);  // vh and vs are lab velocities
          // Hardcoded lower velocity avoids nan results
          double vs = std::max(v/2, 1e-6);
          double tp = pow(mass,2)+pow(mh,2)-0.5*s*(1-vs*vh);
          double tm = pow(mass,2)+pow(mh,2)-0.5*s*(1+vs*vh);

          double aR = 1+3*mh*mh*(s-mh*mh)*Dh2(s);
          double aI = 3*mh*mh*sqrt(s)*Gamma_mh*Dh2(s);

          return pow(lambda,2)/16/M_PI/pow(s,2)/vs *
            (
             (pow(aR,2)+pow(aI,2))*s*vh*vs
             +4*lambda*pow(v0,2)*(aR-lambda*pow(v0,2)/(s-2*pow(mh,2)))
             *log(std::abs(pow(mass,2)-tp)/std::abs(pow(mass,2)-tm))
             +(2*pow(lambda,2)*pow(v0,4)*s*vh*vs)
             /(pow(mass,2)-tm)/(pow(mass,2)-tp));
        }

      private:
        double Gamma_mh, mh, v0, alpha_s, mb, mc, mtau, mt, mZ0, mW;
    };

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
