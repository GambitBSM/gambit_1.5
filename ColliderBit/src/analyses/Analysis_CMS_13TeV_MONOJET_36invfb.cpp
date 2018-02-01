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

        // Register signal region data
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

        // Covariance
        static const double BKGCOV[NUMSR][NUMSR] = {
          {  9.00e+00,  8.73e+00,  6.19e+00,  5.35e+00,  4.62e+00,  3.67e+00,  3.40e+00,  2.61e+00,  1.93e+00,  1.58e+00,  1.40e+00,  1.16e+00,  7.56e-01,  6.48e-01,  4.32e-01,  5.67e-01,  3.29e-01,  1.30e-01,  0.00e+00,  1.81e-01,  4.56e-02, -6.96e-02 },
          {  8.73e+00,  9.00e+00,  6.12e+00,  5.42e+00,  4.62e+00,  3.62e+00,  3.35e+00,  2.48e+00,  1.93e+00,  1.55e+00,  1.30e+00,  1.12e+00,  7.29e-01,  6.48e-01,  3.89e-01,  5.25e-01,  3.48e-01,  1.62e-01,  3.12e-02,  2.06e-01,  7.98e-02, -5.22e-02 },
          {  6.19e+00,  6.12e+00,  5.76e+00,  4.54e+00,  3.97e+00,  3.37e+00,  3.33e+00,  2.48e+00,  1.98e+00,  1.61e+00,  1.50e+00,  1.27e+00,  7.99e-01,  6.53e-01,  5.36e-01,  4.70e-01,  3.37e-01,  1.94e-01,  8.74e-02,  1.55e-01,  6.38e-02,  1.39e-02 },
          {  5.35e+00,  5.42e+00,  4.54e+00,  4.41e+00,  3.47e+00,  2.99e+00,  2.84e+00,  2.14e+00,  1.76e+00,  1.44e+00,  1.31e+00,  1.06e+00,  7.94e-01,  5.71e-01,  4.54e-01,  4.70e-01,  3.07e-01,  1.59e-01,  3.28e-02,  1.17e-01,  4.79e-02, -2.44e-02 },
          {  4.62e+00,  4.62e+00,  3.97e+00,  3.47e+00,  3.61e+00,  2.67e+00,  2.53e+00,  1.82e+00,  1.54e+00,  1.28e+00,  1.25e+00,  9.61e-01,  7.01e-01,  5.47e-01,  3.42e-01,  4.92e-01,  2.67e-01,  1.03e-01,  2.96e-02,  1.39e-01,  9.39e-02, -1.10e-02 },
          {  3.67e+00,  3.62e+00,  3.37e+00,  2.99e+00,  2.67e+00,  3.24e+00,  2.40e+00,  1.78e+00,  1.56e+00,  1.12e+00,  1.12e+00,  1.05e+00,  6.64e-01,  6.19e-01,  3.24e-01,  4.16e-01,  2.20e-01,  1.17e-01,  1.40e-01,  1.86e-01,  6.16e-02,  1.04e-02 },
          {  3.40e+00,  3.35e+00,  3.33e+00,  2.84e+00,  2.53e+00,  2.40e+00,  3.24e+00,  1.86e+00,  1.61e+00,  1.32e+00,  1.19e+00,  1.09e+00,  7.61e-01,  6.34e-01,  5.31e-01,  4.16e-01,  3.07e-01,  1.26e-01,  1.40e-01,  9.29e-02,  3.42e-02,  8.35e-02 },
          {  2.61e+00,  2.48e+00,  2.48e+00,  2.14e+00,  1.82e+00,  1.78e+00,  1.86e+00,  2.25e+00,  1.18e+00,  9.90e-01,  1.01e+00,  8.25e-01,  4.59e-01,  4.92e-01,  3.13e-01,  3.67e-01,  2.47e-01,  8.10e-02,  7.02e-02,  3.23e-02,  7.41e-02,  2.61e-02 },
          {  1.93e+00,  1.93e+00,  1.98e+00,  1.76e+00,  1.54e+00,  1.56e+00,  1.61e+00,  1.18e+00,  1.96e+00,  9.58e-01,  9.58e-01,  8.16e-01,  3.91e-01,  4.59e-01,  3.23e-01,  3.23e-01,  1.37e-01,  1.29e-01,  2.91e-02,  1.02e-01,  1.06e-01,  3.25e-02 },
          {  1.58e+00,  1.55e+00,  1.61e+00,  1.44e+00,  1.28e+00,  1.12e+00,  1.32e+00,  9.90e-01,  9.58e-01,  1.44e+00,  7.06e-01,  6.73e-01,  3.78e-01,  3.55e-01,  2.59e-01,  2.60e-01,  1.17e-01,  7.78e-02,  6.24e-03,  4.64e-02,  4.10e-02,  3.48e-02 },
          {  1.40e+00,  1.30e+00,  1.50e+00,  1.31e+00,  1.25e+00,  1.12e+00,  1.19e+00,  1.01e+00,  9.58e-01,  7.06e-01,  1.44e+00,  6.07e-01,  3.78e-01,  4.03e-01,  2.76e-01,  2.60e-01,  1.46e-01,  1.30e-02, -6.24e-03,  9.29e-02,  1.82e-02, -6.96e-03 },
          {  1.15e+00,  1.12e+00,  1.27e+00,  1.06e+00,  9.61e-01,  1.05e+00,  1.09e+00,  8.25e-01,  8.16e-01,  6.73e-01,  6.07e-01,  1.21e+00,  3.37e-01,  3.70e-01,  2.85e-01,  1.62e-01,  1.68e-01,  1.25e-01,  5.72e-02,  2.37e-02,  4.60e-02,  8.93e-02 },
          {  7.56e-01,  7.29e-01,  7.99e-01,  7.94e-01,  7.01e-01,  6.64e-01,  7.61e-01,  4.59e-01,  3.91e-01,  3.78e-01,  3.78e-01,  3.37e-01,  8.10e-01,  2.23e-01,  2.01e-01,  2.20e-01,  1.76e-01,  5.83e-02,  5.62e-02,  0.00e+00,  1.37e-02,  2.09e-02 },
          {  6.48e-01,  6.48e-01,  6.53e-01,  5.71e-01,  5.47e-01,  6.19e-01,  6.34e-01,  4.92e-01,  4.59e-01,  3.55e-01,  4.03e-01,  3.70e-01,  2.23e-01,  6.40e-01,  1.61e-01,  1.18e-01,  9.76e-02,  2.16e-02,  4.16e-03,  5.16e-02,  2.13e-02,  4.18e-02 },
          {  4.32e-01,  3.89e-01,  5.36e-01,  4.54e-01,  3.42e-01,  3.24e-01,  5.31e-01,  3.13e-01,  3.23e-01,  2.59e-01,  2.76e-01,  2.85e-01,  2.01e-01,  1.61e-01,  5.18e-01,  8.06e-02,  9.66e-02,  6.22e-02,  3.37e-02, -3.10e-03,  1.64e-02,  4.18e-02 },
          {  5.67e-01,  5.25e-01,  4.70e-01,  4.70e-01,  4.92e-01,  4.16e-01,  4.16e-01,  3.67e-01,  3.23e-01,  2.60e-01,  2.60e-01,  1.62e-01,  2.20e-01,  1.18e-01,  8.06e-02,  4.90e-01,  1.71e-02, -1.13e-02, -3.64e-03,  2.11e-02,  3.46e-02, -1.22e-02 },
          {  3.29e-01,  3.48e-01,  3.37e-01,  3.07e-01,  2.67e-01,  2.20e-01,  3.07e-01,  2.47e-01,  1.37e-01,  1.17e-01,  1.46e-01,  1.68e-01,  1.76e-01,  9.76e-02,  9.66e-02,  1.71e-02,  3.72e-01,  5.27e-02,  5.08e-02, -2.89e-02,  2.09e-02,  3.18e-02 },
          {  1.30e-01,  1.62e-01,  1.94e-01,  1.59e-01,  1.03e-01,  1.17e-01,  1.26e-01,  8.10e-02,  1.29e-01,  7.78e-02,  1.30e-02,  1.25e-01,  5.83e-02,  2.16e-02,  6.22e-02, -1.13e-02,  5.27e-02,  2.92e-01,  8.42e-03, -2.09e-02,  8.21e-03, -9.40e-03 },
          {  0.00e+00,  3.12e-02,  8.74e-02,  3.28e-02,  2.96e-02,  1.40e-01,  1.40e-01,  7.02e-02,  2.91e-02,  6.24e-03, -6.24e-03,  5.72e-02,  5.62e-02,  4.16e-03,  3.37e-02, -3.64e-03,  5.08e-02,  8.42e-03,  2.70e-01,  1.79e-02,  5.93e-03,  4.52e-02 },
          {  1.81e-01,  2.06e-01,  1.55e-01,  1.17e-01,  1.39e-01,  1.86e-01,  9.29e-02,  3.23e-02,  1.02e-01,  4.64e-02,  9.29e-02,  2.37e-02,  0.00e+00,  5.16e-02, -3.10e-03,  2.11e-02, -2.89e-02, -2.09e-02,  1.79e-02,  1.85e-01,  3.27e-03,  2.24e-02 },
          {  4.56e-02,  7.98e-02,  6.38e-02,  4.79e-02,  9.39e-02,  6.16e-02,  3.42e-02,  7.41e-02,  1.06e-01,  4.10e-02,  1.82e-02,  4.60e-02,  1.37e-02,  2.13e-02,  1.64e-02,  3.46e-02,  2.09e-02,  8.21e-03,  5.93e-03,  3.27e-03,  1.44e-01,  1.76e-02 },
          { -6.96e-02, -5.22e-02,  1.39e-02, -2.44e-02, -1.10e-02,  1.04e-02,  8.35e-02,  2.61e-02,  3.25e-02,  3.48e-02, -6.96e-03,  8.93e-02,  2.09e-02,  4.18e-02,  4.18e-02, -1.22e-02,  3.18e-02, -9.40e-03,  4.52e-02,  2.24e-02,  1.76e-02,  3.36e-01 }
        };
        set_covariance({{71.6875, 32.1512},{32.1512, 34.5625}});

      }

    };


    // Factory fn
    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_MONOJET_36invfb)


  }
}
