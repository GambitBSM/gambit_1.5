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
    class Analysis_ATLAS_13TeV_ZGammaGrav_CONFNOTE_80invfb : public HEPUtilsAnalysis {
    public:

      Analysis_ATLAS_13TeV_ZGammaGrav_CONFNOTE_80invfb() {
        set_analysis_name(ATLAS_13TeV_ZGammaGrav_CONFNOTE_80invfb);
        set_luminosity(79.8);

        // num_SRaa_SL=0;
        // num_SRaa_SH=0;
        // num_SRaa_WL=0;
        // num_SRaa_WH=0;
        // num_SRaj_L=0;
        // num_SRaj_L200=0;
        // num_SRaj_H=0;

        // NCUTS= 66;
        // for (int i=0;i<NCUTS;i++){
        //   cutFlowVector.push_back(0);
        //   cutFlowVector_str.push_back("");
        // }
      }


      void analyze(const Event* event) {
        HEPUtilsAnalysis::analyze(event);

        // Missing energy
        double met = event->met();
        P4 ptot = event->missingmom();

        // Baseline lepton objects
        ParticlePtrs blElectrons, blMuons;              // Used for SR-2body and SR-3body
        ParticlePtrs baselineElectrons, baselineMuons;  // Used for SR-4body
        for (Particle* electron : event->electrons()) {
          bool crack = (electron->abseta() > 1.37) && (electron->abseta() < 1.52);
          if (electron->pT() > 25. && electron->abseta() < 2.47 && !crack) baselineElectrons.push_back(electron);
        }
        ATLAS::applyTightIDElectronSelection(baselineElectrons);

        const std::vector<double>  a = {0,10.};
        const std::vector<double>  b = {0,10000.};
        const vector<double> cMu={0.89};
        HEPUtils::BinnedFn2D<double> _eff2dMu(a,b,cMu);
        for (Particle* muon : event->muons()) {
          bool hasTrig=has_tag(_eff2dMu, muon->eta(), muon->pT());
          if (muon->pT() > 25. && muon->abseta() < 2.7 && hasTrig) baselineMuons.push_back(muon);
        }

        // Photons
        ParticlePtrs baselinePhotons;
        for (Particle* photon : event->photons()) {
          bool crack = (photon->abseta() > 1.37) && (photon->abseta() < 1.52);
          if (photon->pT() > 25. && photon->abseta() < 2.37 && !crack) baselinePhotons.push_back(photon);
        }

        // Jets
        JetPtrs jets28, jets28_nophooverlap, jets25;
        for (Jet* jet : event->jets()) {
          if (jet->pT() > 30. && fabs(jet->eta()) < 2.8) {
            jets28.push_back(jet);
            jets28_nophooverlap.push_back(jet);
          }
        }

        // Overlap removal
        jetLeptonOverlapRemoval(jets28,baselineElectrons,0.2);
        jetLeptonOverlapRemoval(jets28_nophooverlap,baselineElectrons,0.2);
        jetLeptonOverlapRemoval(jets28,baselinePhotons,0.2); // JetLepton still works for photons
        leptonJetOverlapRemoval(baselineElectrons,jets28_nophooverlap,0.4);
        leptonJetOverlapRemoval(baselineElectrons,jets28,0.4);
        leptonJetOverlapRemoval(baselinePhotons,jets28,0.4); // LeptonJet still works for photons
        leptonJetOverlapRemoval(baselineMuons,jets28,0.4);
        leptonJetOverlapRemoval(baselineMuons,jets28_nophooverlap,0.4);

        // Make |eta| < 2.5 jets
        for (Jet* jet : jets28){
          if (jet->abseta() < 2.5) jets25.push_back(jet);
        }

        //Put signal jets／leptons in pT order
        std::sort(jets28.begin(), jets28.end(), cmpJetsByPt);
        std::sort(jets25.begin(), jets25.end(), cmpJetsByPt);
        std::sort(baselineElectrons.begin(), baselineElectrons.end(), cmpParticlesByPt);

        // Function used to get b jets
        JetPtrs bJets25, bJets28;
        const std::vector<double> c = {0.77};
        HEPUtils::BinnedFn2D<double> _eff2d(a,b,c);
        for (Jet* jet : jets25) {
          bool hasTag=has_tag(_eff2d, jet->eta(), jet->pT());
          if (jet->btag() && hasTag) bJets25.push_back(jet);
        }
        for (Jet* jet : jets28) {
          bool hasTag=has_tag(_eff2d, jet->eta(), jet->pT());
          if (jet->btag() && hasTag) bJets28.push_back(jet);
        }

        // Pre-selection
        bool preSelection2a=false;
        if (baselinePhotons.size()==2 && baselinePhotons[0]->pT() > 75. && baselinePhotons[1]->pT() > 75.)
          preSelection2a=true;
        bool preSelectionSRLaj = false;
        if (baselinePhotons.size()==1 && baselinePhotons[0]->pT() > 145.)
          preSelectionSRLaj=true;
        bool preSelectionSRHaj = false;
        if (baselinePhotons.size()==1 && baselinePhotons[0]->pT() > 400.)
          preSelectionSRHaj=true;

        // Useful variables
        // "Photon-enhanced" HT, with no overlap removal of photons-jets
        double HT = 0.;
        for (Particle* photon : baselinePhotons) HT += photon->pT();
        for (Particle* electron : baselineElectrons) HT += electron->pT();
        for (Particle* muon : baselineMuons) HT += muon->pT();
        for (Jet* jet : jets28_nophooverlap) HT += jet->pT();

        double meff = met;
        for (Particle* photon : baselinePhotons) meff += photon->pT();
        for (Particle* electron : baselineElectrons) meff += electron->pT();
        for (Particle* muon : baselineMuons) meff += muon->pT();
        // Note that meff is only used for aj signal regions -> |jet eta| < 2.5
        for (Jet* jet : jets25) meff += jet->pT();

        // dphimin(a,met)
        double dphimin_amet = 999.;
        for (Particle* photon : baselinePhotons) {
          double dphi_tmp = ptot.deltaPhi(photon->mom());
          if (dphi_tmp < dphimin_amet)dphimin_amet = dphi_tmp;
        }
        double dphimin_j25met = 999.;
        for (Jet* jet : jets25) {
          double dphi_tmp = ptot.deltaPhi(jet->mom());
          if (dphi_tmp < dphimin_j25met)dphimin_j25met = dphi_tmp;
        }
        double dphimin_j28met = 999.;
        for (Jet* jet : jets28) {
          double dphi_tmp = ptot.deltaPhi(jet->mom());
          if (dphi_tmp < dphimin_j28met)dphimin_j28met = dphi_tmp;
        }

        // RT4
        // Only used in aj regions -> use |jet eta| < 2.5
        double RT4 = 0.;
        if (jets25.size() > 3) {
          RT4 = jets25[0]->pT() + jets25[1]->pT() + jets25[2]->pT() + jets25[3]->pT();
        }
        double denom=0.;
        for (Jet* jet : jets25) denom += jet->pT();
        RT4=RT4/denom;

        // Multiplicities
        int nLep = baselineElectrons.size() + baselineMuons.size();
        int nJets25 = jets25.size();

        // All variables are now done
        // Increment signal region variables
        // 2a regions
        if (preSelection2a && met > 150. && HT > 2750 && dphimin_j28met > 0.5)num_SRaa_SL++;
        if (preSelection2a && met > 250. && HT > 2000 && dphimin_j28met > 0.5 && dphimin_amet > 0.5)num_SRaa_SH++;
        if (preSelection2a && met > 150. && HT > 1500 && dphimin_j28met > 0.5)num_SRaa_WL++;
        if (preSelection2a && met > 250. && HT > 1000 && dphimin_j28met > 0.5 && dphimin_amet > 0.5)num_SRaa_WH++;

        // aj regions
        if (preSelectionSRLaj && nJets25 >=5 && nLep == 0 && met > 300. && meff > 2000. && RT4 < 0.90 && dphimin_j25met > 0.5 && dphimin_amet > 0.5)num_SRaj_L++;
        if (preSelectionSRLaj && nJets25 >=5 && nLep == 0 && met > 200. && meff > 2000. && RT4 < 0.90 && dphimin_j25met > 0.5 && dphimin_amet > 0.5)num_SRaj_L200++;
        if (preSelectionSRHaj && nJets25 >=3 && nLep == 0 && met > 400. && meff > 2400. && dphimin_j25met > 0.5 && dphimin_amet > 0.5)num_SRaj_H++;


        /////////////////////////////////////////////////


        /// @todo Init these ONCE, in the constructor

        cutFlowVector_str[0] = "Total ";
        // --------
        cutFlowVector_str[1] = "SR2A--trigger && 2 OS lepton";
        cutFlowVector_str[2] = "SR2ASF--Same flavour";
        cutFlowVector_str[3] = "SR2ASF--mll>111GeV";
        cutFlowVector_str[4] = "SR2ASF--n_{b-jets}=0";
        cutFlowVector_str[5] = "SR2ASF--R_{2l2j}>0.3";
        cutFlowVector_str[6] = "SR2ASF--Delta x<0.07";
        cutFlowVector_str[7] = "SR2ASF--120<MT2<140";
        cutFlowVector_str[8] = "SR2ASF--140<MT2<160";
        cutFlowVector_str[9] = "SR2ASF--160<MT2<180";
        cutFlowVector_str[10] = "SR2ASF--180<MT2";
        // --------
        cutFlowVector_str[11] = "SR2ADF--Different falvour";
        cutFlowVector_str[12] = "SR2ADF--mll>111GeV(only SF)";
        cutFlowVector_str[13] = "SR2ADF--n_{b-jets}=0";
        cutFlowVector_str[14] = "SR2ADF--R_{2l2j}>0.3(only SF)";
        cutFlowVector_str[15] = "SR2ADF--Delta x<0.07";
        cutFlowVector_str[16] = "SR2ADF--120<MT2<140";
        cutFlowVector_str[17] = "SR2ADF--140<MT2<160";
        cutFlowVector_str[18] = "SR2ADF--160<MT2<180";
        cutFlowVector_str[19] = "SR2ADF--180<MT2";
        // --------
        cutFlowVector_str[20] = "SR2BC--trigger && 2 OS lepton";
        cutFlowVector_str[21] = "SR2BSF--Same flavour";
        cutFlowVector_str[22] = "SR2BSF--mll>111GeV or mll<71GeV";
        cutFlowVector_str[23] = "SR2BSF--n_{b-jets}>0 && n_{jets}>1";
        cutFlowVector_str[24] = "SR2BSF--Delta phi_{boost}<1.5";
        cutFlowVector_str[25] = "SR2BSF--120<MT2<140";
        cutFlowVector_str[26] = "SR2BSF--140<MT2";
        // --------
        cutFlowVector_str[27] = "SR2BDF--Different flavour";
        cutFlowVector_str[28] = "SR2BDF--mll>111GeV or mll<71GeV(only SF)";
        cutFlowVector_str[29] = "SR2BDF--n_{b-jets}>0 && n_{jets}>1";
        cutFlowVector_str[30] = "SR2BDF--Delta phi_{boost}<1.5";
        cutFlowVector_str[31] = "SR2BDF--120<MT2<140";
        cutFlowVector_str[32] = "SR2BDF--140<MT2";
        // --------
        cutFlowVector_str[33] = "SR2CSF--n_{b-jets}>0 && n_{jets}>1";
        cutFlowVector_str[34] = "SR2CSF--n_{jets}>2";
        cutFlowVector_str[35] = "SR2CSF--R_{2l}>1.2";
        cutFlowVector_str[36] = "SR2CSF--E_T^{miss}>200GeV";
        cutFlowVector_str[37] = "SR2CSF--110<MT2";
        // --------
        cutFlowVector_str[38] = "SR2CDF--n_{b-jets}>0 && n_{jets}>1";
        cutFlowVector_str[39] = "SR2CDF--n_{jets}>2";
        cutFlowVector_str[40] = "SR2CDF--R_{2l}>1.2";
        cutFlowVector_str[41] = "SR2CDF--E_T^{miss}>200GeV";
        cutFlowVector_str[42] = "SR2CDF--110<MT2";
        // --------
        cutFlowVector_str[57] = "SR4b--MET trigger && 2 OS Leptons ";
        cutFlowVector_str[58] = "SR4b--m_{ll}>10GeV ";
        cutFlowVector_str[59] = "SR4b--PT(l1)<80GeV && PT(l1)<35GeV ";
        cutFlowVector_str[60] = "SR4b--n_{jets}>2 ";
        cutFlowVector_str[61] = "SR4b--PT(j1)>150GeV ";
        cutFlowVector_str[62] = "SR4b--PT(j3)/MET<0.14 ";
        cutFlowVector_str[63] = "SR4b--R_{2l4j}>0.35 ";
        cutFlowVector_str[64] = "SR4b--R_{2l}>12 ";
        cutFlowVector_str[65] = "SR4b--veto on j1 and j2 ";
      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        HEPUtilsAnalysis::add(other);

        Analysis_ATLAS_13TeV_PhotonGGM_36invfb* specificOther
          = dynamic_cast<Analysis_ATLAS_13TeV_PhotonGGM_36invfb*>(other);

        // Here we will add the subclass member variables:
        if (NCUTS != specificOther->NCUTS) NCUTS = specificOther->NCUTS;
        for (int j=0; j<NCUTS; j++) {
          cutFlowVector[j] += specificOther->cutFlowVector[j];
          cutFlowVector_str[j] = specificOther->cutFlowVector_str[j];
        }

        //num_SRaj_L, num_SRaj_L200, num_SRaj_H;
        num_SRaa_SL += specificOther->num_SRaa_SL;
        num_SRaa_SH += specificOther->num_SRaa_SH;
        num_SRaa_WL += specificOther->num_SRaa_WL;
        num_SRaa_WH += specificOther->num_SRaa_WH;
        num_SRaj_L += specificOther->num_SRaj_L;
        num_SRaj_L200 += specificOther->num_SRaj_L200;
        num_SRaj_H += specificOther->num_SRaj_H;
      }


      void collect_results() {

        SignalRegionData results_SRaa_SL;
        results_SRaa_SL.sr_label = "SRaa_SL";
        results_SRaa_SL.n_observed = 0.;
        results_SRaa_SL.n_background = 0.50;
        results_SRaa_SL.background_sys = 0.30;
        results_SRaa_SL.signal_sys = 0.;
        results_SRaa_SL.n_signal = num_SRaa_SL;
        add_result(results_SRaa_SL);

        SignalRegionData results_SRaa_SH;
        results_SRaa_SH.sr_label = "SRaa_SH";
        results_SRaa_SH.n_observed = 0.;
        results_SRaa_SH.n_background = 0.48;
        results_SRaa_SH.background_sys = 0.30;
        results_SRaa_SH.signal_sys = 0.;
        results_SRaa_SH.n_signal = num_SRaa_SH;
        add_result(results_SRaa_SH);

        SignalRegionData results_SRaa_WL;
        results_SRaa_WL.sr_label = "SRaa_WL";
        results_SRaa_WL.n_observed = 6.;
        results_SRaa_WL.n_background = 3.7;
        results_SRaa_WL.background_sys = 1.1;
        results_SRaa_WL.signal_sys = 0.;
        results_SRaa_WL.n_signal = num_SRaa_WL;
        add_result(results_SRaa_WL);

        SignalRegionData results_SRaa_WH;
        results_SRaa_WH.sr_label = "SRaa_WH";
        results_SRaa_WH.n_observed = 1.;
        results_SRaa_WH.n_background = 2.05;
        results_SRaa_WH.background_sys = 0.65;
        results_SRaa_WH.signal_sys = 0.;
        results_SRaa_WH.n_signal = num_SRaa_WH;
        add_result(results_SRaa_WH);

        SignalRegionData results_SRaj_L;
        results_SRaj_L.sr_label = "SRaj_L";
        results_SRaj_L.n_observed = 4.;
        results_SRaj_L.n_background = 1.33;
        results_SRaj_L.background_sys = 0.54;
        results_SRaj_L.signal_sys = 0.;
        results_SRaj_L.n_signal = num_SRaj_L;
        add_result(results_SRaj_L);

        SignalRegionData results_SRaj_L200;
        results_SRaj_L200.sr_label = "SRaj_L200";
        results_SRaj_L200.n_observed = 8.;
        results_SRaj_L200.n_background = 2.68;
        results_SRaj_L200.background_sys = 0.64;
        results_SRaj_L200.signal_sys = 0.;
        results_SRaj_L200.n_signal = num_SRaj_L200;
        add_result(results_SRaj_L200);

        SignalRegionData results_SRaj_H;
        results_SRaj_H.sr_label = "SRaj_H";
        results_SRaj_H.n_observed = 3.;
        results_SRaj_H.n_background = 1.14;
        results_SRaj_H.background_sys = 0.61;
        results_SRaj_H.signal_sys = 0.;
        results_SRaj_H.n_signal = num_SRaj_H;
        add_result(results_SRaj_H);
      }


      void clear() {
        num_SRaa_SL=0;
        num_SRaa_SH=0;
        num_SRaa_WL=0;
        num_SRaa_WH=0;
        num_SRaj_L=0;
        num_SRaj_L200=0;
        num_SRaj_H=0;
        std::fill(cutFlowVector.begin(), cutFlowVector.end(), 0);
      }


    private:

      // Numbers passing cuts
      int num_SRaa_SL, num_SRaa_SH, num_SRaa_WL, num_SRaa_WH, num_SRaj_L, num_SRaj_L200, num_SRaj_H;

      // // Cut flow
      // vector<int> cutFlowVector;
      // vector<string> cutFlowVector_str;
      // int NCUTS;

      /// Jet overlap removal -- discards jets if they are within deltaRMax of a lepton
      void jetLeptonOverlapRemoval(JetPtrs& jetvec, ParticlePtrs& lepvec, double deltaRMax) {
        JetPtrs survivors;
        for (unsigned int itjet = 0; itjet < jetvec.size(); itjet++) {
          bool overlap = false;
          P4& jetmom = jetvec.at(itjet)->mom();
          for (unsigned int itlep = 0; itlep < lepvec.size(); itlep++) {
            P4& lepmom = lepvec.at(itlep)->mom();
            double dR = jetmom.deltaR_eta(lepmom);
            if (fabs(dR) < deltaRMax) {
              overlap = true;
              break;
            }
          }
          if (overlap) continue;
          survivors.push_back(jetvec.at(itjet));
        }
        jetvec = survivors;
      }

      /// Lepton overlap removal -- discards leptons if they are within deltaRMax of a jet
      void leptonJetOverlapRemoval(ParticlePtrs& lepvec, JetPtrs& jetvec, double deltaRMax) {
        ParticlePtrs survivors;
        for (unsigned int itlep = 0; itlep < lepvec.size(); itlep++) {
          bool overlap = false;
          P4& lepmom = lepvec.at(itlep)->mom();
          for (unsigned int itjet= 0; itjet < jetvec.size(); itjet++) {
            P4& jetmom = jetvec.at(itjet)->mom();
            double dR = jetmom.deltaR_eta(lepmom);
            if (fabs(dR) < deltaRMax) {
              overlap = true;
              break;
            }
          }
          if (overlap) continue;
          survivors.push_back(lepvec.at(itlep));
        }
        lepvec = survivors;
      }


    };



    DEFINE_ANALYSIS_FACTORY(ATLAS_13TeV_ZGammaGrav_CONFNOTE_80invfb);


  }
}
