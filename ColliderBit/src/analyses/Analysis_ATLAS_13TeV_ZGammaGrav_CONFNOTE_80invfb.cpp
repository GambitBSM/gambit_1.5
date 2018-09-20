#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "gambit/ColliderBit/mt2_bisect.h"
using namespace std;

namespace Gambit {
  namespace ColliderBit {


    /// @brief ATLAS ZH(->photon+gravitino) (79.8 fb^-1)
    ///
    /// Based on:
    ///  - https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2018-019/
    ///  - https://cds.cern.ch/record/2621481/files/ATLAS-CONF-2018-019.pdf
    ///
    /// @author Andy Buckley
    ///
    /// @warning NOT YET VALIDATED!!!
    ///
    class Analysis_ATLAS_13TeV_ZGammaGrav_CONFNOTE_80invfb : public HEPUtilsAnalysis {
    public:

      Analysis_ATLAS_13TeV_ZGammaGrav_CONFNOTE_80invfb() {
        set_analysis_name(ATLAS_13TeV_ZGammaGrav_CONFNOTE_80invfb);
        set_luminosity(79.8);
        nsig = 0;

        // NCUTS= 66;
        // for (int i=0;i<NCUTS;i++){
        //   cutFlowVector.push_back(0);
        //   cutFlowVector_str.push_back("");
        // }
      }


      void analyze(const Event* event) {
        HEPUtilsAnalysis::analyze(event);

        // Baseline electrons
        ParticlePtrs baselineElectrons;
        for (Particle* e : event->electrons()) {
          const bool crack = e->abseta() > 1.37 && e->abseta() < 1.52;
          if (e->pT() > 10. && e->abseta() < 2.47 && !crack)
            baselineElectrons.push_back(e);
        }
        ATLAS::applyMediumIDElectronSelection(baselineElectrons);

        // Baseline muons
        // NB. medium muon ID for pT > 10 ~ 99%: https://cds.cern.ch/record/2047831/files/ATL-PHYS-PUB-2015-037.pdf
        ParticlePtrs baselineMuons;
        for (Particle* m : event->muons())
          if (m->pT() > 10. && m->abseta() < 2.7 && random_bool(0.99))
            baselineMuons.push_back(m);

        // Photons
        ParticlePtrs baselinePhotons;
        for (Particle* y : event->photons())
          if (y->pT() > 20.)
            baselinePhotons.push_back(y);
        ATLAS::applyPhotonEfficiencyR2(baselinePhotons);
        /// @todo Need to do some explicit isolation? Unlike leptons, there will be significant hadronic photons

        // Jets
        JetPtrs jets;
        for (Jet* j : event->jets())
          if (j->pT() > 20. && j->absrap() < 4.4)
            jets.push_back(j);

        // Overlap removal
        removeOverlap(jets, baselineElectrons, 0.2);
        removeOverlap(baselineElectrons, jets, 0.4);
        removeOverlap(baselineElectrons, jets, 0.4);
        removeOverlap(baselinePhotons, baselineElectrons, 0.4);
        removeOverlap(baselinePhotons, baselineMuons, 0.4);
        removeOverlap(baselineJets, baselinePhotons, 0.4);

        // Put objects in pT order
        sortByPt(jets);
        sortByPt(baselineElectrons);
        sortByPt(baselineMuons);
        sortByPt(baselinePhotons);

        // Missing energy
        const double met = event->met();
        const P4 pmiss = event->missingmom();


        /////////////////


        // There must be a prompt photon
        if (baselinePhotons.empty()) return;

        // Find the Z system
        if (baselineElectrons.size() + baselineMuons.size() != 2) return; //< must be exactly two leptons
        if (!baselineElectrons.empty() && !baselineMuons.empty()) return; //< the two leptons must be same-flavour
        const ParticlePtrs& leps = baselineElectrons.empty() ? baselineMuons : baselineElectrons;

        // Check lepton pTs and require small delta_phi
        if (leps[0]->pT() < 25 || leps[1]->pT() < 20) return;
        if (deltaPhi(leps[0]->mom(), leps[1]->mom()) > 1.4) return;
        // The dilepton mass must be within 10 GeV of the Z mass
        const P4 pZ = leps[0]->mom() + leps[1]->mom();
        if (fabs(pZ.m() - 91.2) > 10.) return;

        // MET and jet requirements
        if (met < 95) return;
        if (!jets.empty() && jets[0]->pT() > 30) return;

        // Require separation of the Z and the MET+photon(s) systems
        const P4 pYMET = pmiss + baselinePhotons.[0]->mom() +
          (baselinePhotons.size() > 1 ? baselinePhotons[1]->mom() : P4());
        if (deltaPhi(pZ, pYMET) < 2.8) return;
        if (fabs(pZ.pT()-pYMET.pT())/pYMET.pT() > 0.2) return;

        // Signal count
        nsig += 1;

      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        HEPUtilsAnalysis::add(other);

        Analysis_ATLAS_13TeV_PhotonGGM_36invfb* specificOther
          = dynamic_cast<Analysis_ATLAS_13TeV_PhotonGGM_36invfb*>(other);

        // // Here we will add the subclass member variables:
        // if (NCUTS != specificOther->NCUTS) NCUTS = specificOther->NCUTS;
        // for (int j=0; j<NCUTS; j++) {
        //   cutFlowVector[j] += specificOther->cutFlowVector[j];
        //   cutFlowVector_str[j] = specificOther->cutFlowVector_str[j];
        // }

        nsig += specificOther->nsig;
      }


      void collect_results() {

        SignalRegionData results;
        results_SRaa_SL.sr_label = "SR";
        results_SRaa_SL.n_observed = 3.;
        results_SRaa_SL.n_background = 2.1;
        results_SRaa_SL.background_sys = 0.5;
        results_SRaa_SL.signal_sys = 0.;
        results_SRaa_SL.n_signal = nsig;
        add_result(results);

      }


      void clear() {
        nsig = 0;
        // std::fill(cutFlowVector.begin(), cutFlowVector.end(), 0);
      }


    private:

      // Numbers passing cuts
      int nsig = 0;

      // // Cut flow
      // vector<int> cutFlowVector;
      // vector<string> cutFlowVector_str;
      // int NCUTS;


      /// Overlap removal -- discard from first list if within deltaRMax of any from the second list
      template<typename MOMS1, typename MOMS2>
      void removeOverlap(MOMS1& momstofilter, const MOMS2& momstocmp, double deltaRMax) {
        ifilter_reject(momstofilter, [&](const MOMS1::value_type& mom1) {
            for (const MOMS2::value_type& mom2 : momstocmp) {
              const double dR = mom1.deltaR_eta(mom2);
              if (dR < deltaRMax) return true;
            }
            return false;
          }, false);
      }


    };



    DEFINE_ANALYSIS_FACTORY(ATLAS_13TeV_ZGammaGrav_CONFNOTE_80invfb);


  }
}
