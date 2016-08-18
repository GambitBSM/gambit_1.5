// -*- C++ -*-
#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "Eigen/Eigen"

namespace Gambit {
  namespace ColliderBit {

    using namespace std;
    using namespace HEPUtils;


    /// @brief ATLAS Run 2 0-lepton jet+MET SUSY analysis, with 13/fb of data
    ///
    /// Based on:
    ///   https://cds.cern.ch/record/2206252
    ///   https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2016-078/
    ///
    /// Recursive jigsaw reconstruction signal regions are currently not included
    ///
    class Analysis_ATLAS_13TeV_0LEP_13invfb : public HEPUtilsAnalysis {
    public:

      // Numbers passing cuts
      static const size_t NUMSR = 13;
      double _srnums[13] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};


      // vector<int> cutFlowVector;
      // vector<string> cutFlowVector_str;
      // size_t NCUTS; //=16;


      Analysis_ATLAS_13TeV_0LEP_13invfb() {
        set_luminosity(13.3);
        // NCUTS=60;
        // for (size_t i=0;i<NCUTS;i++){
        //   cutFlowVector.push_back(0);
        //   cutFlowVector_str.push_back("");
        // }
      }


      void analyze(const Event* event) {

        HEPUtilsAnalysis::analyze(event);


        // Missing energy
        const P4 pmiss = event->missingmom();
        const double met = event->met();


        // Get baseline jets
        vector<Jet*> baselineJets;
        for (Jet* jet : event->jets())
          if (jet->pT() > 20. && jet->abseta() < 2.8)
            baselineJets.push_back(jet);
        /// @todo Drop b-tag if pT < 50 GeV or |eta| > 2.5

        // Get baseline electrons
        vector<Particle*> baselineElectrons;
        for (Particle* electron : event->electrons())
          if (electron->pT() > 10. && electron->abseta() < 2.47)
            baselineElectrons.push_back(electron);

        // Get baseline muons
        vector<Particle*> baselineMuons;
        for (Particle* muon : event->muons())
          if (muon->pT() > 10. && muon->abseta() < 2.7)
            baselineMuons.push_back(muon);


        // Remove any |eta| < 0.2 jet within dR = 0.2 of an electron
        /// @todo Unless b-tagged
        vector<const Jet*> signalJets;
        for (const Jet* j : baselineJets)
          if (j->abseta() > 2.8 ||
              all_of(baselineElectrons.begin(), baselineElectrons.end(),
                     [&](const Particle* e){ return deltaR_rap(*e, *j) > 0.2; }))
            signalJets.push_back(j);

        // Remove electrons with dR = 0.4 of surviving |eta| < 2.8 jets
        /// @todo Actually only within 0.2--0.4
        vector<const Particle*> signalElectrons;
        for (const Particle* e : baselineElectrons)
          if (all_of(signalJets.begin(), signalJets.end(),
                     [&](const Jet* j){ return j->abseta() > 2.8 || deltaR_rap(*e, *j) > 0.4; }))
            signalElectrons.push_back(e);
        // Apply electron ID selection
        /// @todo Use *loose* electron selection
        ATLAS::applyMediumIDElectronSelection(signalElectrons);

        // Remove muons with dR = 0.4 of surviving |eta| < 2.8 jets
        /// @todo Note says that dR is in rap rather than eta
        /// @todo Actually only within 0.2--0.4
        /// @todo Within 0.2, discard the *jet* based on jet track vs. muon criteria
        vector<const Particle*> signalMuons;
        for (const Particle* m : baselineMuons)
          if (all_of(signalJets.begin(), signalJets.end(),
                     [&](const Jet* j){ return j->abseta() > 2.8 || deltaR_rap(*m, *j) > 0.4; }))
            signalMuons.push_back(m);

        // The subset of jets with pT > 50 GeV is used for several calculations
        vector<const Jet*> signalJets50;
        for (const Jet* j : signalJets)
          if (j->pT() > 50) signalJets50.push_back(j);


        ////////////////////////////////
        // Calculate common variables and cuts

        // Multiplicities
        const size_t nElectrons = signalElectrons.size();
        const size_t nMuons = signalMuons.size();
        const size_t nJets = signalJets.size();
        const size_t nJets50 = signalJets50.size();

        // HT-related quantities (calculated over all >20 GeV jets)
        double sumptj = 0;
        for (const Jet* j : signalJets) sumptj += j->pT();
        const double HT = sumptj;
        const double sqrtHT = sqrt(HT);
        const double met_sqrtHT = met/sqrtHT;

        // Meff-related quantities (calculated over >50 GeV jets only)
        double sumptj50_4 = 0, sumptj50_5 = 0, sumptj50_6 = 0, sumptj50_incl = 0;
        for (size_t i = 0; i < signalJets50.size(); ++i) {
          const Jet* j = signalJets50[i];
          if (i < 4) sumptj50_4 += j->pT();
          if (i < 5) sumptj50_5 += j->pT();
          if (i < 6) sumptj50_6 += j->pT();
          sumptj50_incl += j->pT();
        }
        const double meff_4 = met + sumptj50_4;
        const double meff_5 = met + sumptj50_5;
        const double meff_6 = met + sumptj50_6;
        const double meff_incl = met + sumptj50_incl;
        const double met_meff_4 = met / meff_4;
        const double met_meff_5 = met / meff_5;
        const double met_meff_6 = met / meff_6;

        // Jet |eta|s
        double etamax_2 = 0;
        for (size_t i = 0; i < min(2lu,signalJets50.size()); ++i)
          etamax_2 = max(etamax_2, signalJets[i]->abseta());
        double etamax_4 = etamax_2;
        for (size_t i = 2; i < min(4lu,signalJets50.size()); ++i)
          etamax_4 = max(etamax_4, signalJets[i]->abseta());
        double etamax_6 = etamax_4;
        for (size_t i = 4; i < min(6lu,signalJets50.size()); ++i)
          etamax_6 = max(etamax_6, signalJets[i]->abseta());

        // Jet--MET dphis
        double dphimin_123 = DBL_MAX, dphimin_more = DBL_MAX;
        for (size_t i = 0; i < min(3lu,signalJets50.size()); ++i)
          dphimin_123 = min(dphimin_123, acos(cos(signalJets50[i]->phi() - pmiss.phi())));
        for (size_t i = 3; i < signalJets50.size(); ++i)
          dphimin_more = min(dphimin_more, acos(cos(signalJets50[i]->phi() - pmiss.phi())));

        // Jet aplanarity
        /// @todo Computed over all jets, all >50 jets, or 4,5,6 jets? Currently using all (> 20) jets
        Eigen::Matrix3d momtensor = Eigen::Matrix3d::Zero();
        double norm = 0;
        for (const Jet* jet : signalJets) {
          const P4& p4 = jet->mom();
          norm += p4.p2();
          for (size_t i = 0; i < 3; ++i) {
            const double pi = (i == 0) ? p4.px() : (i == 1) ? p4.py() : p4.pz();
            for (size_t j = 0; j < 3; ++j) {
              const double pj = (j == 0) ? p4.px() : (j == 1) ? p4.py() : p4.pz();
              momtensor(i,j) += pi*pj;
            }
          }
        }
        momtensor *= norm;
        const double mineigenvalue = momtensor.eigenvalues().real().minCoeff();
        const double aplanarity = 1.5 * mineigenvalue;


        ////////////////////////////////
        // Fill signal regions


        const bool leptonCut = (nElectrons == 0 && nMuons == 0);
        const bool metCut = (met > 250.);
        if (nJets50 >= 2 && leptonCut && metCut) {

          // 2 jet regions
          if (dphimin_123 > 0.8 && dphimin_more > 0.4) {
            if (signalJets[1]->pT() > 200 && etamax_2 < 0.8) { //< implicit pT[0] cut
              if (met_sqrtHT > 14 && meff_incl >  800) _srnums[0] += 1;
            }
            if (signalJets[1]->pT() > 250 && etamax_2 < 1.2) { //< implicit pT[0] cut
              if (met_sqrtHT > 16 && meff_incl > 1200) _srnums[1] += 1;
              if (met_sqrtHT > 18 && meff_incl > 1600) _srnums[2] += 1;
              if (met_sqrtHT > 20 && meff_incl > 2000) _srnums[3] += 1;
            }
          }

          // 3 jet region
          if (nJets50 >= 3 && dphimin_123 > 0.4 && dphimin_more > 0.2) {
            if (signalJets[0]->pT() > 600 && signalJets[2]->pT() > 50) { //< implicit pT[1] cut
              if (met_sqrtHT > 16 && meff_incl > 1200) _srnums[4] += 1;
            }
          }

          // 4 jet regions (note implicit pT[1,2] cuts)
          if (nJets >= 4 && dphimin_123 > 0.4 && dphimin_more > 0.4 && signalJets[0]->pT() > 200 && aplanarity > 0.04) {
            if (signalJets[3]->pT() > 100 && etamax_4 < 1.2 && met_meff_4 > 0.25 && meff_incl > 1000) _srnums[5] += 1;
            if (signalJets[3]->pT() > 100 && etamax_4 < 2.0 && met_meff_4 > 0.25 && meff_incl > 1400) _srnums[6] += 1;
            if (signalJets[3]->pT() > 100 && etamax_4 < 2.0 && met_meff_4 > 0.20 && meff_incl > 1800) _srnums[7] += 1;
            if (signalJets[3]->pT() > 150 && etamax_4 < 2.0 && met_meff_4 > 0.20 && meff_incl > 2200) _srnums[8] += 1;
            if (signalJets[3]->pT() > 150 &&                   met_meff_4 > 0.20 && meff_incl > 2600) _srnums[9] += 1;
          }

          // 5 jet region (note implicit pT[1,2,3] cuts)
          if (nJets >= 5 && dphimin_123 > 0.4 && dphimin_more > 0.2 && signalJets[0]->pT() > 500) {
            if (signalJets[4]->pT() > 50 && met_meff_5 > 0.3 && meff_incl > 1400) _srnums[10] += 1;
          }

          // 6 jet regions (note implicit pT[1,2,3,4] cuts)
          if (nJets >= 6 && dphimin_123 > 0.4 && dphimin_more > 0.2 && signalJets[0]->pT() > 200 && aplanarity > 0.08) {
            if (signalJets[5]->pT() >  50 && etamax_6 < 2.0 && met_meff_6 > 0.20 && meff_incl > 1800) _srnums[11] += 1;
            if (signalJets[5]->pT() > 100 &&                   met_meff_6 > 0.15 && meff_incl > 2200) _srnums[12] += 1;
          }

        }

        // cutFlowVector_str[0] = "No cuts ";
        // cutFlowVector_str[1] = "2j: MET > 160 GeV and jet pT ";
        // cutFlowVector_str[2] = "2j: dPhiMin > 0.4 ";
        // cutFlowVector_str[3] = "2j: met/sqrt(HT) > 15 ";
        // cutFlowVector_str[4] = "2j: meff_incl > 1200 ";
        // cutFlowVector_str[5] = "2j: meff_incl > 1600 ";
        // cutFlowVector_str[6] = "3j: MET > 160 and jet pT ";
        // cutFlowVector_str[7] = "3j: dPhiMin > 0.4 ";
        // cutFlowVector_str[8] = "3j: met/meff3j > 0.3 ";
        // cutFlowVector_str[9] = "3j: met/meff_incl > 2200. ";
        // cutFlowVector_str[10] = "4jlm: MET > 160 and jet pT ";
        // cutFlowVector_str[11] = "4jlm: dPhiMin > 0.4 ";
        // cutFlowVector_str[12] = "4jlm: dPhiMin2 > 0.2 ";
        // cutFlowVector_str[13] = "4jlm: met/sqrt(HT) > 10 ";
        // cutFlowVector_str[14] = "4jlm: meff incl > 700 ";
        // cutFlowVector_str[15] = "4jl: meff incl > 1000 ";
        // cutFlowVector_str[16] = "4jt: met/meff4j > 0.25 ";
        // cutFlowVector_str[17] = "4jt: meff incl > 2200 ";
        // cutFlowVector_str[18] = "5j: MET > 160 and jet pT ";
        // cutFlowVector_str[19] = "5j: dPhiMin > 0.4 ";
        // cutFlowVector_str[20] = "5j: dPhiMin2 > 0.2 ";
        // cutFlowVector_str[21] = "5j: met/meff5j > 0.2 ";
        // cutFlowVector_str[22] = "5j: meff incl > 1200. ";
        // cutFlowVector_str[23] = "6jl: MET >  160 and jet pT  ";
        // cutFlowVector_str[24] = "6jl: dPhiMin > 0.4 ";
        // cutFlowVector_str[25] = "6jl: dPhiMin2 > 0.2 ";
        // cutFlowVector_str[26] = "6jl: met/meff6j > 0.2 ";
        // cutFlowVector_str[27] = "6jl: meff incl > 900. ";
        // cutFlowVector_str[28] = "6jt: met/meff6j > 0.25 ";
        // cutFlowVector_str[29] = "6jt: meff incl > 1500. ";

        // for (size_t j=0;j<NCUTS;j++){
        //   if(
        //      (j==0) ||

        //      (j==1 && signalJets.size()>1 && signalJets[0]->pT()>130. && signalJets[1]->pT()>60. && metCut && leptonCut) ||

        //      (j==2 && signalJets.size()>1 && signalJets[0]->pT()>130. && signalJets[1]->pT()>60. && metCut && dPhiMin2j>0.4 && leptonCut) ||

        //      (j==3 && signalJets.size()>1 && signalJets[0]->pT()>130. && signalJets[1]->pT()>60. && metCut && dPhiMin2j>0.4 && met/sqrt(HT)>15. && leptonCut) ||

        //      (j==4 && signalJets.size()>1 && signalJets[0]->pT()>130. && signalJets[1]->pT()>60. && metCut && dPhiMin2j>0.4 && met/sqrt(HT)>15. && leptonCut && meff_incl>1200.) ||

        //      (j==5 && signalJets.size()>1 && signalJets[0]->pT()>130. && signalJets[1]->pT()>60. && metCut && dPhiMin2j>0.4 && met/sqrt(HT)>15. && leptonCut && meff_incl>1600.) ||

        //      (j==6 && signalJets.size()>2 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && metCut && leptonCut) ||

        //      (j==7 && signalJets.size()>2 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && metCut && leptonCut && dPhiMin3j > 0.4) ||

        //      (j==8 && signalJets.size()>2 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && metCut && leptonCut && dPhiMin3j > 0.4 && met/meff3j>0.3) ||

        //      (j==9 && signalJets.size()>2 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && metCut && leptonCut && dPhiMin3j > 0.4 && met/meff3j>0.3 && meff_incl>2200.) ||

        //      (j==10 && signalJets.size() > 3 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && metCut) ||

        //      (j==11 && signalJets.size() > 3 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && metCut && dPhiMin4 > 0.4) ||

        //      (j==12 && signalJets.size() > 3 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && metCut && dPhiMin4 > 0.4 && dPhiMin2 > 0.2) ||

        //      (j==13 && signalJets.size() > 3 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && metCut && dPhiMin4 > 0.4 && dPhiMin2 > 0.2 && met/sqrt(HT) > 10.) ||

        //      (j==14 && signalJets.size() > 3 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && metCut && dPhiMin4 > 0.4 && dPhiMin2 > 0.2 && met/sqrt(HT) > 10. && meff_incl > 700.) ||

        //      (j==15 && signalJets.size() > 3 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && metCut && dPhiMin4 > 0.4 && dPhiMin2 > 0.2 && met/sqrt(HT) > 10. && meff_incl > 1000.) ||

        //      (j==16 && signalJets.size() > 3 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && metCut && dPhiMin4 > 0.4 && dPhiMin2 > 0.2 && met/meff4j>0.25) ||

        //      (j==17 && signalJets.size() > 3 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && metCut && dPhiMin4 > 0.4 && dPhiMin2 > 0.2 && met/meff4j>0.25 && meff_incl > 2200.) ||

        //      //Start 5j signal regions

        //      (j==18 && signalJets.size() > 4 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && signalJets.at(4)->pT()>60. && metCut && leptonCut) ||

        //      (j==19 && signalJets.size() > 4 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && signalJets.at(4)->pT()>60. && metCut && leptonCut && dPhiMin4 > 0.4) ||

        //      (j==20 && signalJets.size() > 4 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && signalJets.at(4)->pT()>60. && metCut && leptonCut && dPhiMin4 > 0.4 && dPhiMin2 > 0.2) ||

        //      (j==21 && signalJets.size() > 4 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && signalJets.at(4)->pT()>60. && metCut && leptonCut && dPhiMin4 > 0.4 && dPhiMin2 > 0.2 && met/meff4j > 0.25) ||

        //      (j==22 && signalJets.size() > 4 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && signalJets.at(4)->pT()>60. && metCut && leptonCut && dPhiMin4 > 0.4 && dPhiMin2 > 0.2 && met/meff4j > 0.25 && meff_incl > 1200.) ||

        //      //Start 6jl region

        //      (j==23 && signalJets.size() > 5 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && signalJets.at(4)->pT()>60. && signalJets.at(5)->pT()>60. && leptonCut && metCut) ||

        //      (j==24 && signalJets.size() > 5 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && signalJets.at(4)->pT()>60. && signalJets.at(5)->pT()>60. && leptonCut && metCut && dPhiMin4>0.4) ||

        //      (j==25 && signalJets.size() > 5 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && signalJets.at(4)->pT()>60. && signalJets.at(5)->pT()>60. && leptonCut && metCut && dPhiMin4>0.4 && dPhiMin2 > 0.2) ||

        //      (j==26 && signalJets.size() > 5 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && signalJets.at(4)->pT()>60. && signalJets.at(5)->pT()>60. && leptonCut && metCut && dPhiMin4>0.4 && dPhiMin2 > 0.2 && met/meff6j>0.2) ||

        //      (j==27 && signalJets.size() > 5 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && signalJets.at(4)->pT()>60. && signalJets.at(5)->pT()>60. && leptonCut && metCut && dPhiMin4>0.4 && dPhiMin2 > 0.2 && met/meff6j>0.2 && meff_incl > 900.) ||

        //      (j==28 && signalJets.size() > 5 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && signalJets.at(4)->pT()>60. && signalJets.at(5)->pT()>60. && leptonCut && metCut && dPhiMin4>0.4 && dPhiMin2 > 0.2 && met/meff6j>0.2 && meff_incl > 900. && met/meff6j>0.25) ||

        //      (j==29 && signalJets.size() > 5 && signalJets.at(0)->pT()>130. && signalJets.at(1)->pT()>60. && signalJets.at(2)->pT()>60. && signalJets.at(3)->pT()>60. && signalJets.at(4)->pT()>60. && signalJets.at(5)->pT()>60. && leptonCut && metCut && dPhiMin4>0.4 && dPhiMin2 > 0.2 && met/meff6j>0.2 && meff_incl > 900. && met/meff6j>0.25 && meff_incl>1500.)

        //      ){

        //     cutFlowVector[j]++;

        //   }
        // }

      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        HEPUtilsAnalysis::add(other);

        Analysis_ATLAS_13TeV_0LEP_13invfb* specificOther = dynamic_cast<Analysis_ATLAS_13TeV_0LEP_13invfb*>(other);

        // // Here we will add the subclass member variables:
        // if (NCUTS != specificOther->NCUTS) NCUTS = specificOther->NCUTS;
        // for (size_t j = 0; j < NCUTS; j++) {
        //   cutFlowVector[j] += specificOther->cutFlowVector[j];
        //   cutFlowVector_str[j] = specificOther->cutFlowVector_str[j];
        // }

        for (size_t i = 0; i < NUMSR; ++i)
          _srnums[i] += specificOther->_srnums[i];
      }


      void collect_results() {

        // cout << "------------------------------------------------------------------------------------------------------------------------------ " << endl;
        // cout << "CUT FLOW: ATLAS R2 0-lepton paper "<<endl;
        // cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        // cout << right << setw(40) << "CUT" << setw(20) << "RAW" << setw(20) << "SCALED" << setw(20) << "%" << setw(20) << "clean adj RAW"<< setw(20) << "clean adj %" << endl;
        // const double scale_by = 1;
        // for (size_t j=0; j<NCUTS; j++) {
        //   cout << right << setw(40) << cutFlowVector_str[j].c_str() << setw(20) << cutFlowVector[j] << setw(20) <<
        //     cutFlowVector[j]*scale_by << setw(20) << 100.*cutFlowVector[j]/cutFlowVector[0] << "%" << setw(20) <<
        //     cutFlowVector[j]*scale_by << setw(20) << 100.*cutFlowVector[j]/cutFlowVector[0]<< "%" << endl;
        // }
        // cout << "------------------------------------------------------------------------------------------------------------------------------ " << endl;


        // Now fill a results object with the results for each SR
        // Numbers are taken from CONF note
        static const string ANAME = "Analysis_ATLAS_13TeV_0LEP_13invfb";
        add_result(SignalRegionData(ANAME, "meff-2j-800", 650, {_srnums[0], 0.}, {610., 10.}));
        /// @todo More SRs
      }

    };


    // Factory fn
    DEFINE_ANALYSIS_FACTORY(ATLAS_13TeV_0LEP_13invfb)


  }
}
