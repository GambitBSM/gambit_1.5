// -*- C++ -*-
#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include "gambit/ColliderBit/analyses/Cutflow.hpp"
#include "gambit/ColliderBit/CMSEfficiencies.hpp"
// #include "Eigen/Eigen"

namespace Gambit {
  namespace ColliderBit {

    using namespace std;
    using namespace HEPUtils;


    /// @brief CMS Run 2 monojet analysis (no W/Z region) with 36/fb of data
    ///
    /// @todo Add W/Z region with AKT8 jets and 2/1 n-subjettiness ratio cut
    ///
    class Analysis_CMS_13TeV_MONOJET_36invfb : public HEPUtilsAnalysis {
    public:

      static const size_t NUMSR = 22;
      double _srnums[NUMSR];
      Cutflow _cutflow;

      Analysis_CMS_13TeV_MONOJET_36invfb()
      // : _cutflow("CMS monojet 13 TeV", {"Njet >= 3", "HT > 300", "HTmiss > 300", "Nmuon = 0", "Nelectron = 0", "Nhadron = 0 (no-op)", "Dphi_htmiss_j1", "Dphi_htmiss_j2", "Dphi_htmiss_j3", "Dphi_htmiss_j4"})
      {
        for (double& n : _srnums) n = 0;
        set_luminosity(35.9);
      }


      void analyze(const Event* event) {

        HEPUtilsAnalysis::analyze(event);
        // _cutflow.fillinit();

        // Require large MET
        const P4 pmiss = event->missingmom();
        const double met = pmiss.pT();
        if (met < 250) return; //< VETO

        // Record a trigger weight; we can aggregate this rather than wastefully random-vetoing
        const double trigweight = (met < 350) ? 0.97 : 1.0;

        // Veto on isolated leptons and photons
        for (const Particle* e : event->electrons()) if (e->pT() > 10) return; //< VETO
        for (const Particle* m : event->muons()) if (m->pT() > 10) return; //< VETO
        for (const Particle* t : event->taus()) if (t->pT() > 18) return; //< VETO
        for (const Particle* y : event->photons()) if (y->pT() > 15 && y->abseta() < 2.5) return; //< VETO

        // Get jets
        vector<const Jet*> jets4;
        for (const Jet* jet : event->jets())
          if (jet->pT() > 20) jets4.push_back(jet);

        // Veto if there are any b-tagged jets (reduce top background)
        for (const Jet* jet : jets4) {
          if (jet->abseta() > 2.4) continue;
          const double btag_rate = jet->btag() ? 0.8 : jet->ctag() ? 0.4 : 0.1;
          if (rand01() < btag_rate) return; //< VETO
        }

        // Get the 4 leading jets > 3 GeV, and veto if pTmiss is too close to them
        for (size_t i = 0; i < 4; ++i) {
          if (i >= jets4.size()) break;
          if (jets4[i]->pT() < 30) break;
          if (deltaPhi(jets4[i]->mom(), pmiss)) return; //< VETO
        }

        // Now the signal regions, but we'll just look at the monojet one
        if (jets4.empty()) return;
        if (jets4[0]->pT() < 100*GeV || jets4[0]->abseta() > 2.4) return;

        // Identify the ptmiss bin and fill the counter
        const static vector<double> metedges = {250, 280, 310, 340, 370, 400, 430, 470, 510, 550, 590,
                                                640, 690, 740, 790, 840, 900, 960, 1020, 1090, 1160, 1250};
        const int i_sr = binIndex(met, metedges, true);
        if (i_sr >= 0) _srnums[i_sr] += trigweight;

      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        HEPUtilsAnalysis::add(other);
        Analysis_CMS_13TeV_MONOJET_36invfb* specificOther = dynamic_cast<Analysis_CMS_13TeV_MONOJET_36invfb*>(other);
        for (size_t i = 0; i < NUMSR; ++i)
          _srnums[i] += specificOther->_srnums[i];
      }


      /// Register results objects with the results for each SR; obs & bkg numbers from the CONF note
      void collect_results() {
        //cout << _cutflow << endl;
        static const string ANAME = "CMS_13TeV_MONOJET_36invfb";
        static const double OBSNUM[NUMSR] = {
          // 162, 130, 97.8, 84.8, 65.2, 53.5, 53.9, 41.4, 34.3, 28.1, 27.5,
          // 20.4, 16.6, 12.5, 8.94, 10.1, 6.62, 5.19, 4.35, 2.84, 3.44, 6.39
          136865, 74340, 42540, 25316, 15653, 10092, 8298, 4906, 2987, 2032, 1514,
          926, 557, 316, 233, 172, 101, 65, 46, 26, 31, 29
        };
        static const double BKGNUM[NUMSR] = {
          134500, 73400, 42320, 25490, 15430, 10160, 8480, 4865, 2970, 1915, 1506,
          844, 526, 325, 223, 169, 107, 88.1, 52.8, 25.0, 25.5, 26.9
        };
        static const double BKGERR[NUMSR] = {
          3.0, 3.0, 2.4, 2.1, 1.9, 1.8, 1.8, 1.5, 1.4, 1.2, 1.2,
          1.1, 0.9, 0.8, 0.72, 0.7, 0.61, 0.54, 0.52, 0.43, 0.38, 0.58
        };
        for (size_t ibin = 0; ibin < NUMSR; ++ibin) {
          stringstream ss; ss << "sr-" << ibin;
          add_result(SignalRegionData(ANAME, ss.str(), OBSNUM[ibin], {_srnums[ibin],  0.}, {BKGNUM[ibin], BKGERR[ibin]}));
        }
      }

    };


    // Factory fn
    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_MONOJET_36invfb)


  }
}
