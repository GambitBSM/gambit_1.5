///
///  \author Yang Zhang
///  \date 2018 Nov
///
///  *********************************************

// Based on http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/index.html
// Searches for pair production of charginos and top squarks in final states with two oppositely charged leptons in proton-proton collisions at √s= 13 TeV

#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>
#include <fstream>

#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include "gambit/ColliderBit/CMSEfficiencies.hpp"
#include "gambit/ColliderBit/mt2_bisect.h"
#include "gambit/ColliderBit/analyses/Cutflow.hpp"

using namespace std;

namespace Gambit {
  namespace ColliderBit {


    class Analysis_CMS_13TeV_2OSLEP_chargino_36invfb : public HEPUtilsAnalysis {
    public:
    
      // Counters for the number of accepted events for each signal region
      static const size_t NUMSR = 70; 
      double _srnums[NUMSR];
      Cutflow _cutflow;

      Analysis_CMS_13TeV_2OSLEP_chargino_36invfb():
      _cutflow("CMS 2-lep chargino->slepton 13 TeV", {"Two_OC_leptons", "Third_lepton_veto", "mll_20", "mll_mZ_15", "PTmiss_140"})
      {
        set_analysis_name("CMS_13TeV_2OSLEP_chargino_36invfb");
        set_luminosity(35.9);
        for (size_t i = 0; i < NUMSR; ++i) _srnums[i] = 0;
      }

      struct ptComparison {
        bool operator() (HEPUtils::Particle* i,HEPUtils::Particle* j) {return (i->pT()>j->pT());}
      } comparePt;
      
      void analyze(const HEPUtils::Event* event) 
      {
        // Baseline objects
        HEPUtilsAnalysis::analyze(event);
        double met = event->met();
        _cutflow.fillinit();
        
        // Apply electron efficiency and collect baseline electrons
        //@note Numbers digitized from https://twiki.cern.ch/twiki/pub/CMSPublic/SUSMoriond2017ObjectsEfficiency/2d_full_pteta_el_17010_ttbar.pdf
        //@note The efficiency map has been extended to cover the low-pT region, using the efficiencies from BuckFast (CMSEfficiencies.hpp)
        const vector<double> aEl={0., 0.8, 1.442, 1.556, 2., 2.4, DBL_MAX};   // Bin edges in eta
        const vector<double> bEl={0., 10., 20., 25., 30., 40., 50., DBL_MAX}; // Bin edges in pT. Assume flat efficiency above 200, where the CMS map stops.
        const vector<double> cEl={                 
                          // pT: (0,10),  (10,20),  (20,25),  (25,30),  (30,40),  (40,50),  (50,inf)     
                                    0.0,    0.95,    0.463,     0.537,    0.597,   0.648,    0.701,  // eta: (0, 0.8)
                                    0.0,    0.95,    0.424,     0.520,    0.603,   0.635,    0.700, // eta: (0.8, 1.4429
                                    0.0,    0.95,    0.041,     0.044,    0.053,   0.049,    0.053, // eta: (1.442, 1.556)
                                    0.0,    0.85,    0.271,     0.343,    0.439,   0.508,    0.557, // eta: (1.556, 2)
                                    0.0,    0.85,    0.262,     0.336,    0.411,   0.463,    0.513, // eta: (2, 2.5) 
                                    0.0,    0.0,     0.0,       0.0,      0.0,     0.0,      0.0,   // eta > 2.5
                                  };
        HEPUtils::BinnedFn2D<double> _eff2dEl(aEl,bEl,cEl);
        vector<HEPUtils::Particle*> baselineElectrons;
        for (HEPUtils::Particle* electron : event->electrons()) 
        {
          bool isEl=has_tag(_eff2dEl, electron->eta(), electron->pT());
          if (isEl && electron->pT()>15. && fabs(electron->eta())<2.4) baselineElectrons.push_back(electron);
        }


        // Apply muon efficiency and collect baseline muons
        //@note Numbers digitized from https://twiki.cern.ch/twiki/pub/CMSPublic/SUSMoriond2017ObjectsEfficiency/2d_full_pteta_mu_17010_ttbar.pdf
        //@note The efficiency map has been extended to cover the low-pT region, using the efficiencies from BuckFast (CMSEfficiencies.hpp)
        const vector<double> aMu={0., 0.9, 1.2, 2.1, 2.4, DBL_MAX};   // Bin edges in eta
        const vector<double> bMu={0., 10., 20., 25., 30., 40., 50., DBL_MAX};  // Bin edges in pT. Assume flat efficiency above 200, where the CMS map stops.
        const vector<double> cMu={
                           // pT:   (0,10),  (10,20),  (20,25),  (25,30),  (30,40),  (40,50),  (50,inf)     
                                     0.0,    0.950,    0.692,    0.759,    0.813,    0.852,    0.896,  // eta: (0, 0.9)
                                     0.0,    0.950,    0.679,    0.752,    0.803,    0.856,    0.897,  // eta: (0.9, 1.2)
                                     0.0,    0.950,    0.706,    0.757,    0.801,    0.836,    0.886,  // eta: (1.2, 2.1)
                                     0.0,    0.950,    0.647,    0.711,    0.753,    0.800,    0.806,  // eta: (2.1, 2.4)
                                     0.0,    0.0,      0.0,      0.0,      0.0,      0.0,      0.0,    // eta > 2.4
                                 };
        HEPUtils::BinnedFn2D<double> _eff2dMu(aMu,bMu,cMu);
        vector<HEPUtils::Particle*> baselineMuons;
        for (HEPUtils::Particle* muon : event->muons()) 
        {
          bool isMu=has_tag(_eff2dMu, muon->eta(), muon->pT());
          if (isMu && muon->pT()>15. && fabs(muon->eta())<2.4) baselineMuons.push_back(muon);
        }

        // Baseline jets
        vector<HEPUtils::Jet*> baselineJets;
        for (HEPUtils::Jet* jet : event->jets()) 
        {
          if (jet->pT()>20. &&fabs(jet->eta())<2.4) baselineJets.push_back(jet);
        }

        
        // Signal leptons = electrons + muons
        vector<HEPUtils::Particle*> signalLeptons;
        signalLeptons=baselineElectrons;
        signalLeptons.insert(signalLeptons.end(),baselineMuons.begin(),baselineMuons.end());
        sort(signalLeptons.begin(),signalLeptons.end(),comparePt);
        
        // Signal jets
        const vector<double> aBJet={0,10.};
        const vector<double> bBJet={0,30., 40., 50., 70., 80., 90., 100.,150., 200., 10000.};
        const vector<double> cBJet={0.63, 0.705, 0.745, 0.76, 0.775, 0.79,0.795, 0.805, 0.795, 0.76};
        HEPUtils::BinnedFn2D<double> _eff2d(aBJet,bBJet,cBJet);
        vector<HEPUtils::Jet*> signalJets,signalBJets;
        for (size_t iJet=0;iJet<baselineJets.size();iJet++) {
          bool hasTag=has_tag(_eff2d, baselineJets.at(iJet)->eta(), baselineJets.at(iJet)->pT());
          bool ovelap=false;
          for (size_t iEl=0;iEl<signalLeptons.size();iEl++) {
            if (signalLeptons.at(iEl)->mom().deltaR_eta(baselineJets.at(iJet)->mom())<0.4) ovelap=true;
          }
          if ( not ovelap ) signalJets.push_back(baselineJets.at(iJet));
          if ( (not ovelap) and baselineJets.at(iJet)->btag() and  hasTag) signalBJets.push_back(baselineJets.at(iJet));
        }
        int nj  = signalJets.size();
        int nbj = signalBJets.size();
        
        
        // Tow OC lepton
        if (signalLeptons.size()<2) return;
        if (signalLeptons[0]->pid()*signalLeptons[1]->pid()>0) return;
        if (signalLeptons[0]->pT()<25. or signalLeptons[1]->pT()<20.) return;
        _cutflow.fill(1);
        
        // Third lepton veto
        if (signalLeptons.size()>2) return;
        _cutflow.fill(2);
        
        // m_{ll} > 20 GeV
        double mll=(signalLeptons[0]->mom()+signalLeptons[1]->mom()).m();
        if (mll<20) return;
        _cutflow.fill(3);
        
        // |m_{ll}-m_Z| >15 GeV for ee and mumu events
        bool same_flavor = signalLeptons[0]->abspid() == signalLeptons[1]->abspid();        
        if (same_flavor and fabs(mll-91.2)<15 ) return;
        _cutflow.fill(4);

        // MET>140 GeV
        if (met<140) return;
        _cutflow.fill(5);
        
        // Mt2
        double pLep1[3] = {signalLeptons[0]->mass(), signalLeptons[0]->mom().px(), signalLeptons[0]->mom().py()};
        double pLep2[3] = {signalLeptons[1]->mass(), signalLeptons[1]->mom().px(), signalLeptons[1]->mom().py()};
        double pMiss[3] = {0., event->missingmom().px(), event->missingmom().py() };
        mt2_bisect::mt2 mt2_calc;
        mt2_calc.set_momenta(pLep1,pLep2,pMiss);
        mt2_calc.set_mn(0.0);
        double mT2 = mt2_calc.get_mt2();
        
        // In the order of x axis
        // http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/CMS-SUS-17-010_Figure-aux_001.png
        if ( not same_flavor ) {
            if (met<200) {
                if (nbj==0 and nj ==0) { //SR1_{0tag}^{0jet}
                    if      (mT2< 20) ++_srnums[ 0];
                    else if (mT2< 40) ++_srnums[ 1];
                    else if (mT2< 60) ++_srnums[ 2];
                    else if (mT2< 80) ++_srnums[ 3];
                    else if (mT2<100) ++_srnums[ 4];
                    else if (mT2<120) ++_srnums[ 5];
                    else              ++_srnums[ 6];
                }
                if (nbj==0 and nj >=1) { //SR1_{0tag}^{jets}
                    if      (mT2< 20) ++_srnums[14];
                    else if (mT2< 40) ++_srnums[15];
                    else if (mT2< 60) ++_srnums[16];
                    else if (mT2< 80) ++_srnums[17];
                    else if (mT2<100) ++_srnums[18];
                    else if (mT2<120) ++_srnums[19];
                    else              ++_srnums[20];      
                }
            } else if (met<300){
                if (nbj==0 and nj ==0) { //SR2_{0tag}^{0jet}
                    if      (mT2< 20) ++_srnums[28];
                    else if (mT2< 40) ++_srnums[29];
                    else if (mT2< 60) ++_srnums[30];
                    else if (mT2< 80) ++_srnums[31];
                    else if (mT2<100) ++_srnums[32];
                    else if (mT2<120) ++_srnums[33];
                    else              ++_srnums[34];
                }
                if (nbj==0 and nj >=1) { //SR2_{0tag}^{jets}
                    if      (mT2< 20) ++_srnums[42];
                    else if (mT2< 40) ++_srnums[43];
                    else if (mT2< 60) ++_srnums[44];
                    else if (mT2< 80) ++_srnums[45];
                    else if (mT2<100) ++_srnums[46];
                    else if (mT2<120) ++_srnums[47];
                    else              ++_srnums[48];      
                }
            } else {
                if (nbj==0 and nj >=0) { //SR3_{0tag}
                    if      (mT2< 20) ++_srnums[56];
                    else if (mT2< 40) ++_srnums[57];
                    else if (mT2< 60) ++_srnums[58];
                    else if (mT2< 80) ++_srnums[59];
                    else if (mT2<100) ++_srnums[60];
                    else if (mT2<120) ++_srnums[61];
                    else              ++_srnums[62];        
                }
            }
        } else {
            if (met<200) {
                if (nbj==0 and nj ==0) { //SR1_{0tag}^{0jet}
                    if      (mT2< 20) ++_srnums[ 7];
                    else if (mT2< 40) ++_srnums[ 8];
                    else if (mT2< 60) ++_srnums[ 9];
                    else if (mT2< 80) ++_srnums[10];
                    else if (mT2<100) ++_srnums[11];
                    else if (mT2<120) ++_srnums[12];
                    else              ++_srnums[13];
                }
                if (nbj==0 and nj >=1) { //SR1_{0tag}^{jets}
                    if      (mT2< 20) ++_srnums[21];
                    else if (mT2< 40) ++_srnums[22];
                    else if (mT2< 60) ++_srnums[23];
                    else if (mT2< 80) ++_srnums[24];
                    else if (mT2<100) ++_srnums[25];
                    else if (mT2<120) ++_srnums[26];
                    else              ++_srnums[27]; 
                         
                }
            } else if (met<300){
                if (nbj==0 and nj ==0) { //SR2_{0tag}^{0jet}
                    if      (mT2< 20) ++_srnums[35];
                    else if (mT2< 40) ++_srnums[36];
                    else if (mT2< 60) ++_srnums[37];
                    else if (mT2< 80) ++_srnums[38];
                    else if (mT2<100) ++_srnums[39];
                    else if (mT2<120) ++_srnums[40];
                    else              ++_srnums[41];
                }
                if (nbj==0 and nj >=1) { //SR2_{0tag}^{jets}
                    if      (mT2< 20) ++_srnums[49];
                    else if (mT2< 40) ++_srnums[50];
                    else if (mT2< 60) ++_srnums[51];
                    else if (mT2< 80) ++_srnums[52];
                    else if (mT2<100) ++_srnums[53];
                    else if (mT2<120) ++_srnums[54];
                    else              ++_srnums[55];
                }
            } else {
                if (nbj==0 and nj >=0) { //SR3_{0tag}
                    if      (mT2< 20) ++_srnums[63];
                    else if (mT2< 40) ++_srnums[64];
                    else if (mT2< 60) ++_srnums[65];
                    else if (mT2< 80) ++_srnums[66];
                    else if (mT2<100) ++_srnums[67];
                    else if (mT2<120) ++_srnums[68];
                    else              ++_srnums[69];
                }
            }
        }

      }    

      void add(BaseAnalysis* other) 
      {
        // The base class add function handles the signal region vector and total # events.
        
        HEPUtilsAnalysis::add(other);

        Analysis_CMS_13TeV_2OSLEP_chargino_36invfb* specificOther
                = dynamic_cast<Analysis_CMS_13TeV_2OSLEP_chargino_36invfb*>(other);

        // Here we will add the subclass member variables:
        for (size_t i = 0; i < NUMSR; ++i)
            _srnums[i] += specificOther->_srnums[i];
         
      }


      virtual void collect_results() 
      {
        cout << _cutflow << endl;
        for (size_t ibin = 0; ibin < NUMSR; ++ibin) {
            stringstream ss; ss << "sr-" << ibin;
            cout << "sr-" << ibin << "\t" << _srnums[ibin] << endl;
        }
        
//        // Observed event counts
//        static const double OBSNUM[SR_size_cov] = {
//          57., 29., 2., 0., 9., 5., 1.
//        };
//        // Background estimates
//        static const double BKGNUM[SR_size_cov] = {
//          54.9, 21.6, 6., 2.5, 7.6, 5.6, 1.3
//        };
//        // Background uncertainties, same-flavor signal regions
//        static const double BKGERR[SR_size_cov] = {
//          7., 5.6, 1.9, 0.9, 2.8, 1.6, 0.4,
//        };
//        
        
        
//        for (size_t ibin = 0; ibin < NUMSR; ++ibin) {
//          stringstream ss; ss << "sr-" << ibin;
//          add_result(SignalRegionData(ss.str(), OBSNUM[ibin], {_srnums[ibin],  0.}, {BKGNUM[ibin], BKGERR[ibin]}));
//        }
//        
        // Fake
        add_result(SignalRegionData("SR1", 57., {_srnums[0], 0.}, {54.9, 7.}));
        add_result(SignalRegionData("SR2", 29., {_srnums[1], 0.}, {21.6, 5.6}));
        add_result(SignalRegionData("SR3", 2.,  {_srnums[2], 0.}, {6., 1.9}));
        add_result(SignalRegionData("SR4", 0.,  {_srnums[3], 0.}, {2.5, 0.9}));
        add_result(SignalRegionData("SR5", 9.,  {_srnums[4], 0.}, {7.6, 2.8}));
        add_result(SignalRegionData("SR6", 5.,  {_srnums[5], 0.}, {5.6, 1.6}));
        add_result(SignalRegionData("SR7", 1.,  {_srnums[6], 0.}, {1.3, 0.4}));

//        // Covariance matrix
//        static const vector< vector<double> > BKGCOV = {
//          { 52.8, 12.7,  3.0,  1.2,  4.5,  5.1,  1.2},
//          { 12.7, 41.4,  3.6,  2.0,  2.5,  2.0,  0.7},
//          {  3.0,  3.6,  1.6,  0.6,  0.4,  0.3,  0.1},
//          {  1.2,  2.0,  0.6,  1.1,  0.3,  0.1,  0.1},
//          {  4.5,  2.5,  0.4,  0.3,  6.5,  1.8,  0.4},
//          {  5.1,  2.0,  0.3,  0.1,  1.8,  2.4,  0.4},
//          {  1.2,  0.7,  0.1,  0.1,  0.4,  0.4,  0.2},
//        };        

//        set_covariance(BKGCOV);
      }
      


    protected:
      void clear() {

        for(size_t i=0;i<NUMSR;i++) { _srnums[i]=0; }

      }

    };

    // Factory fn
    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_2OSLEP_chargino_36invfb)


  }
}
