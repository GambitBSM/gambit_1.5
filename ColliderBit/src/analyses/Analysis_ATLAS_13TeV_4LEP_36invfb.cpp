///
///  \author Anders Kvellestad
///  \date 2018 June
///  *********************************************

// Based on https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-21/

// - So far only the signal regions with zero taus have been implemented
// - Some cuts have not yet been implemented



#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>
#include <algorithm>
#include <fstream>

#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "gambit/ColliderBit/mt2_bisect.h"

using namespace std;

namespace Gambit
{
  namespace ColliderBit 
  {

    class Analysis_ATLAS_13TeV_4LEP_36invfb : public HEPUtilsAnalysis 
    {

    protected:

      // Counters for the number of accepted events for each signal region
      std::map<string,double> _numSR = {
        {"SR0A", 0},
        {"SR0B", 0},
        {"SR0C", 0},
        {"SR0D", 0},
      }; 

    private:

      vector<int> cutFlowVector1;
      vector<string> cutFlowVector1_str;
      size_t NCUTS1;
      // vector<double> cutFlowVector1ATLAS_200_100;
      // double xsec1ATLAS_200_100; 

      struct ptComparison 
      {
        bool operator() (HEPUtils::Particle* i,HEPUtils::Particle* j) {return (i->pT()>j->pT());}
      } comparePt;
      
      struct ptJetComparison 
      {
        bool operator() (HEPUtils::Jet* i,HEPUtils::Jet* j) {return (i->pT()>j->pT());}
      } compareJetPt;

      // Jet lepton overlap removal
      // Discards jets if they are within DeltaRMax of a lepton
      void JetLeptonOverlapRemoval(vector<HEPUtils::Jet*>& jets, vector<HEPUtils::Particle*>& leptons, double DeltaRMax) 
      {
        vector<HEPUtils::Jet*> survivors;
        for(HEPUtils::Jet* jet : jets)
        {
          bool overlap = false;
          for(HEPUtils::Particle* lepton : leptons) 
          {
            double dR = jet->mom().deltaR_eta(lepton->mom());
            if(fabs(dR) <= DeltaRMax) overlap = true;
          }
          if(!overlap) survivors.push_back(jet);
        }
        jets = survivors;
        return;
      }

      // Lepton jet overlap removal
      // Discards leptons if they are within DeltaRMax of a jet
      void LeptonJetOverlapRemoval(vector<HEPUtils::Particle*>& leptons, vector<HEPUtils::Jet*>& jets, double DeltaRMax) 
      {
        vector<HEPUtils::Particle*> survivors;
        for(HEPUtils::Particle* lepton : leptons)
        {
          bool overlap = false;
          for(HEPUtils::Jet* jet : jets) 
          {
            double dR = jet->mom().deltaR_eta(lepton->mom());
            if(fabs(dR) <= DeltaRMax) overlap = true;
          }
          if(!overlap) survivors.push_back(lepton);
        }
        leptons = survivors;
        return;
      }


      // Removes a lepton from the leptons1 vector if it forms an OS pair with a
      // lepton in leptons2 and the pair has a mass in the range (m_low, m_high).
      void removeOSPairsInMassRange(vector<HEPUtils::Particle*>& leptons1, vector<HEPUtils::Particle*>& leptons2, double m_low, double m_high)
      {
        vector<HEPUtils::Particle*> l1_survivors;
        for(HEPUtils::Particle* l1 : leptons1)
        {
          bool survived = true;
          for(HEPUtils::Particle* l2 : leptons2)
          {
            if(l2 == l1) continue;
            if (l1->pid()*l2->pid() < 0.) 
            {
              double m = (l1->mom() + l2->mom()).m();
              if ((m >= m_low) && (m <= m_high)) 
              {
                survived = false;
                break;
              }
            }
          }
          if(survived) l1_survivors.push_back(l1);
        }
        leptons1 = l1_survivors;
        return;
      }


    public:

      Analysis_ATLAS_13TeV_4LEP_36invfb() 
      {

        set_analysis_name("ATLAS_13TeV_4LEP_36invfb");
        set_luminosity(36.1);

        NCUTS1=22;

        // xsec1ATLAS_200_100=1807.4;
        for (size_t i=0;i<NCUTS1;i++)
        {
          cutFlowVector1.push_back(0);
          // cutFlowVector1ATLAS_200_100.push_back(0);
          cutFlowVector1_str.push_back("");
        }

      }

      void analyze(const HEPUtils::Event* event) 
      {
        HEPUtilsAnalysis::analyze(event);

        // Baseline objects
        vector<HEPUtils::Particle*> baselineElectrons;
        vector<HEPUtils::Particle*> baselineMuons;
        vector<HEPUtils::Jet*> baselineJets;
        double met = event->met();


        for (HEPUtils::Particle* electron : event->electrons()) 
        {
          if (electron->pT()>7. && electron->abseta()<2.47) baselineElectrons.push_back(electron);
        }
        ATLAS::applyLooseIDElectronSelectionR2(baselineElectrons);

        for (HEPUtils::Particle* muon : event->muons()) 
        {
          if (muon->pT()>5. && muon->abseta()<2.7) baselineMuons.push_back(muon);
        }
        // Missing: Apply "medium" ID criteria 

        for (HEPUtils::Jet* jet : event->jets()) 
        {
          if (jet->pT()>20. && jet->abseta()<2.8) baselineJets.push_back(jet);
        }
        // Missing: Some additional requirements for jets with pT < 60 and abseta < 2.4 (see paper)


        // Overlap removal
        // 1) Remove taus within DeltaR = 0.2 of an electron or muon

        // 2) Remove electron sharing an ID track with a muon

        // 3) Remove jets within DeltaR = 0.2 of electron
        JetLeptonOverlapRemoval(baselineJets, baselineElectrons, 0.2);

        // 4) Remove electrons within DeltaR = 0.4 of a jet
        LeptonJetOverlapRemoval(baselineElectrons, baselineJets, 0.4);

        // 5) Remove jets with < 3 assocated tracks if a muon is within DeltaR = 0.2 
        //    *or* the muon is a track in the jet.

        // 6) Remove muons within DeltaR = 0.4 of jet
        LeptonJetOverlapRemoval(baselineMuons, baselineJets, 0.4);

        // 7) Remove jets within DetlaR = 0.4 of a "medium" tau


        // Suppress low-mass particle decays
        vector<HEPUtils::Particle*> baselineLeptons;
        baselineLeptons = baselineElectrons;
        baselineLeptons.insert(baselineLeptons.end(), baselineMuons.begin(), baselineMuons.end());
        // - Remove low-mass OS pairs
        removeOSPairsInMassRange(baselineElectrons, baselineLeptons, 0.0, 4.0);
        removeOSPairsInMassRange(baselineMuons, baselineLeptons, 0.0, 4.0);
        // - Remove SFOS pairs in the mass range (8.4, 10.4) GeV
        removeOSPairsInMassRange(baselineElectrons, baselineElectrons, 8.4, 10.4);
        removeOSPairsInMassRange(baselineMuons, baselineMuons, 8.4, 10.4);


        // Signal objects
        vector<HEPUtils::Jet*> signalJets = baselineJets;
        vector<HEPUtils::Particle*> signalElectrons = baselineElectrons;
        vector<HEPUtils::Particle*> signalMuons = baselineMuons;
        vector<HEPUtils::Particle*> signalLeptons;
        signalLeptons = signalElectrons;
        signalLeptons.insert(signalLeptons.end(), signalMuons.begin(), signalMuons.end());

        // Missing: pT-dependent isolation criteria for signal leptons (see paper)

        // Sort by pT
        sort(signalJets.begin(), signalJets.end(), compareJetPt); 
        sort(signalLeptons.begin(), signalLeptons.end(), comparePt);

        // Count signal leptons and jets
        size_t nSignalElectrons = signalElectrons.size();
        size_t nSignalMuons = signalMuons.size();
        size_t nSignalLeptons = signalLeptons.size();
        size_t nSignalJets = signalJets.size();

        // Get OS and SFOS pairs
        vector<vector<HEPUtils::Particle*>> SFOSpairs = getSFOSpair(signalLeptons);
        vector<vector<HEPUtils::Particle*>> OSpairs = getOSpair(signalLeptons);


        // Z requirements
        vector<double> SFOSpair_masses;
        for (vector<HEPUtils::Particle*> pair : SFOSpairs)
        {
          SFOSpair_masses.push_back( (pair.at(0)->mom() + pair.at(1)->mom()).m() );
        }
        std::sort(SFOSpair_masses.begin(), SFOSpair_masses.end(), std::greater<double>());

        bool Z1 = false;
        bool Z2 = false;
        bool Zlike = false;
        for(double m : SFOSpair_masses)
        {
          if (!Z1 && (m > 81.2) && (m < 101.2))
          {
            Z1 = true;
          }          
          else if (Z1 && (m > 61.2) && (m < 101.2)) 
          {
            Z2 = true;
          }
        }
        if (Z1) Zlike = true;
        // Missing: Also check Z-like combinations of SFOS+L and SFOS+SFOS (see paper)


        // Effective mass (met + pT of all signal leptons + pT of all jets with pT>40 GeV)
        double meff = met;
        for (HEPUtils::Particle* l : signalLeptons)
        {
          meff += l->pT();
        }
        for (HEPUtils::Jet* jet : signalJets) 
        {
          if(jet->pT()>40.) meff += jet->pT();
        }


        // Signal Regions

        // SR0A
        // Missing: check for zero taus
        if (nSignalLeptons >= 4 && !Zlike && meff > 600.) _numSR["SR0A"]++;
        // if (nSignalLeptons >= 4 && !Zlike && meff > 600.)
        // {
        //   cout << "DEBUG: " << "--- Got event for SR0A ---" << endl;
        //   cout << "DEBUG: " << "  leptons: " << nSignalLeptons << ", electrons: " << nSignalElectrons << ", muons: " << nSignalMuons << endl;
        //   cout << "DEBUG: " << "  jets: " << nSignalJets << endl;
        //   cout << "DEBUG: " << "  meff = " << meff << endl;
        //   cout << "DEBUG: " << "  nSFOSpairs = " << SFOSpairs.size() << endl;
        //   for (double mass : SFOSpair_masses)
        //   {
        //     cout << "DEBUG: " << "  pair mass: " << mass << endl;
        //   }

        //   _numSR["SR0A"]++;
        // }

        // SR0B
        // Missing: check for zero taus
        if (nSignalLeptons >= 4 && !Zlike && meff > 1100.) _numSR["SR0B"]++;
        // if (nSignalLeptons >= 4 && !Zlike && meff > 1100.)
        // {
        //   cout << "DEBUG: " << "--- Got event for SR0B ---" << endl;
        //   cout << "DEBUG: " << "  leptons: " << nSignalLeptons << ", electrons: " << nSignalElectrons << ", muons: " << nSignalMuons << endl;
        //   cout << "DEBUG: " << "  jets: " << nSignalJets << endl;
        //   cout << "DEBUG: " << "  meff = " << meff << endl;
        //   cout << "DEBUG: " << "  nSFOSpairs = " << SFOSpairs.size() << endl;
        //   for (double mass : SFOSpair_masses)
        //   {
        //     cout << "DEBUG: " << "  pair mass: " << mass << endl;
        //   }

        //   _numSR["SR0B"]++;
        // }


        // SR0C
        // Missing: check for zero taus
        if (nSignalLeptons >= 4 && Z1 && Z2 && met > 50.) _numSR["SR0C"]++;
        // if (nSignalLeptons >= 4 && Z1 && Z2 && met > 50.)
        // {
        //   cout << "DEBUG: " << "--- Got event for SR0C ---" << endl;
        //   cout << "DEBUG: " << "  leptons: " << nSignalLeptons << ", electrons: " << nSignalElectrons << ", muons: " << nSignalMuons << endl;
        //   cout << "DEBUG: " << "  jets: " << nSignalJets << endl;
        //   cout << "DEBUG: " << "  met = " << met << endl;
        //   cout << "DEBUG: " << "  nSFOSpairs = " << SFOSpairs.size() << endl;
        //   for (double mass : SFOSpair_masses)
        //   {
        //     cout << "DEBUG: " << "  pair mass: " << mass << endl;
        //   }

        //   _numSR["SR0C"]++;
        // }

        // SR0D
        // Missing: check for zero taus
        if (nSignalLeptons >= 4 && Z1 && Z2 && met > 100.) _numSR["SR0D"]++;
        // if (nSignalLeptons >= 4 && Z1 && Z2 && met > 100.)
        // {
        //   cout << "DEBUG: " << "--- Got event for SR0D ---" << endl;
        //   cout << "DEBUG: " << "  leptons: " << nSignalLeptons << ", electrons: " << nSignalElectrons << ", muons: " << nSignalMuons << endl;
        //   cout << "DEBUG: " << "  jets: " << nSignalJets << endl;
        //   cout << "DEBUG: " << "  met = " << met << endl;
        //   cout << "DEBUG: " << "  nSFOSpairs = " << SFOSpairs.size() << endl;
        //   for (double mass : SFOSpair_masses)
        //   {
        //     cout << "DEBUG: " << "  pair mass: " << mass << endl;
        //   }

        //   _numSR["SR0D"]++;
        // }

        // Missing: signal regions SR1 (3L1T) and SR2 (2L2T)
      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        
        HEPUtilsAnalysis::add(other);

        Analysis_ATLAS_13TeV_4LEP_36invfb* specificOther
                = dynamic_cast<Analysis_ATLAS_13TeV_4LEP_36invfb*>(other);

        // Here we will add the subclass member variables:
        if (NCUTS1 != specificOther->NCUTS1) NCUTS1 = specificOther->NCUTS1;
        for (size_t j = 0; j < NCUTS1; j++) {
          cutFlowVector1[j] += specificOther->cutFlowVector1[j];
          cutFlowVector1_str[j] = specificOther->cutFlowVector1_str[j];
        }

        for (auto& el : _numSR) { 
          el.second += specificOther->_numSR[el.first];
        }

      }

      // This function can be overridden by the derived SR-specific classes
      virtual void collect_results() {

        // add_result(SignalRegionData("SR label", n_obs, {s, s_sys}, {b, b_sys}));
        add_result(SignalRegionData("SR0A", 13., {_numSR["SR0A"], 0.}, {10.2, 2.1}));
        add_result(SignalRegionData("SR0B",  2., {_numSR["SR0B"], 0.}, {1.31, 0.24}));
        add_result(SignalRegionData("SR0C", 47., {_numSR["SR0C"], 0.}, {37., 9.}));
        add_result(SignalRegionData("SR0D", 10., {_numSR["SR0D"], 0.}, {4.1, 0.7}));
      }

      // Construct collection of SFOS pairs
      vector<vector<HEPUtils::Particle*>> getSFOSpair(vector<HEPUtils::Particle*> leptons) {
        vector<vector<HEPUtils::Particle*>> SFOSpair_container;
        for (size_t iLe1=0; iLe1<leptons.size(); iLe1++) {
          for (size_t iLe2=iLe1+1; iLe2<leptons.size(); iLe2++) {
            if (leptons.at(iLe1)->abspid()==leptons.at(iLe2)->abspid() && leptons.at(iLe1)->pid()!=leptons.at(iLe2)->pid()) {
              vector<HEPUtils::Particle*> SFOSpair;
              SFOSpair.push_back(leptons.at(iLe1));
              SFOSpair.push_back(leptons.at(iLe2));
              SFOSpair_container.push_back(SFOSpair);
            }
          }
        }
        return SFOSpair_container;
      }

      // Construct collection of OS pairs
      vector<vector<HEPUtils::Particle*>> getOSpair(vector<HEPUtils::Particle*> leptons) {
        vector<vector<HEPUtils::Particle*>> OSpair_container;
        for (size_t iLe1=0;iLe1<leptons.size();iLe1++) {
          for (size_t iLe2=iLe1+1; iLe2<leptons.size(); iLe2++) {
            if (leptons.at(iLe1)->pid()*leptons.at(iLe2)->pid()<0.) {
              vector<HEPUtils::Particle*> OSpair;
              OSpair.push_back(leptons.at(iLe1));
              OSpair.push_back(leptons.at(iLe2));
              OSpair_container.push_back(OSpair);
            }
          }
        }
        return OSpair_container;
      }


    protected:
      void clear() {
        for (auto& el : _numSR) { el.second = 0.;}
        std::fill(cutFlowVector1.begin(), cutFlowVector1.end(), 0);
      }

    };

    // Factory fn
    DEFINE_ANALYSIS_FACTORY(ATLAS_13TeV_4LEP_36invfb)


  }
}
