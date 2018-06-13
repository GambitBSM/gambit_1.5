///
///  \author Are Raklev
///  \date 2018 June
///
///  *********************************************

// Based on 1806.04030

#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>

#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
//#include "gambit/ColliderBit/mt2_bisect.h"
//#include "gambit/ColliderBit/lester_mt2_bisect.h"

using namespace std;

// Known issues:

namespace Gambit {
  namespace ColliderBit {

    
    class Analysis_ATLAS_13TeV_3b_24invfb : public HEPUtilsAnalysis {
    private:

      // Variables that hold the number of events passing signal region cuts

      double _meff160_ETmiss0, _meff160_ETmiss20;

      vector<int> cutFlowVector;
      vector<string> cutFlowVector_str;
      size_t NCUTS;
      vector<double> cutFlowVectorATLAS_130;

    public:

      static bool sortByPT(HEPUtils::Jet* jet1, HEPUtils::Jet* jet2) { return (jet1->pT() > jet2->pT()); }
      
      Analysis_ATLAS_13TeV_3b_24invfb() {

        set_analysis_name("ATLAS_13TeV_3b_24invfb");
        set_luminosity(24.3);

        // Set number of events passing cuts to zero upon initialisation
        _meff160_ETmiss0=0; _meff160_ETmiss20=0;
        
        NCUTS=9;

        for(size_t i=0;i<NCUTS;i++){
          cutFlowVector.push_back(0);
          cutFlowVectorATLAS_130.push_back(0);
          cutFlowVector_str.push_back("");
        }

      }

      // The following section copied from Analysis_ATLAS_1LEPStop_20invfb.cpp
      void JetLeptonOverlapRemoval(vector<HEPUtils::Jet*> &jetvec, vector<HEPUtils::Particle*> &lepvec, double DeltaRMax) {
        //Routine to do jet-lepton check
        //Discards jets if they are within DeltaRMax of a lepton

        vector<HEPUtils::Jet*> Survivors;

        for(unsigned int itjet = 0; itjet < jetvec.size(); itjet++) {
          bool overlap = false;
          HEPUtils::P4 jetmom=jetvec.at(itjet)->mom();
          for(unsigned int itlep = 0; itlep < lepvec.size(); itlep++) {
            HEPUtils::P4 lepmom=lepvec.at(itlep)->mom();
            double dR;

            dR=jetmom.deltaR_eta(lepmom);

            if(fabs(dR) <= DeltaRMax) overlap=true;
          }
          if(overlap) continue;
          Survivors.push_back(jetvec.at(itjet));
        }
        jetvec=Survivors;

        return;
      }

      void LeptonJetOverlapRemoval(vector<HEPUtils::Particle*> &lepvec, vector<HEPUtils::Jet*> &jetvec) {
        //Routine to do lepton-jet check
        //Discards leptons if they are within dR of a jet as defined in analysis paper

        vector<HEPUtils::Particle*> Survivors;

        for(unsigned int itlep = 0; itlep < lepvec.size(); itlep++) {
          bool overlap = false;
          HEPUtils::P4 lepmom=lepvec.at(itlep)->mom();
          for(unsigned int itjet= 0; itjet < jetvec.size(); itjet++) {
            HEPUtils::P4 jetmom=jetvec.at(itjet)->mom();
            double dR;
            double DeltaRMax = std::max(0.1,std::min(0.4, 0.04 + 10 / lepmom.pT()));
            dR=jetmom.deltaR_eta(lepmom);

            if(fabs(dR) <= DeltaRMax) overlap=true;
          }
          if(overlap) continue;
          Survivors.push_back(lepvec.at(itlep));
        }
        lepvec=Survivors;

        return;
      }

      void SpecialLeptonJetOverlapRemoval(vector<HEPUtils::Particle*> &lepvec, vector<HEPUtils::Jet*> &jetvec) {
        //Routine to do lepton-jet check
        //Discards leptons if they are within dR of a jet as defined in analysis paper

        vector<HEPUtils::Particle*> Survivors;

        for(unsigned int itlep = 0; itlep < lepvec.size(); itlep++) {
          bool overlap = false;
          HEPUtils::P4 lepmom=lepvec.at(itlep)->mom();
          for(unsigned int itjet= 0; itjet < jetvec.size(); itjet++) {
            HEPUtils::P4 jetmom=jetvec.at(itjet)->mom();
            double dR;
            double DeltaRMax = std::min(0.4, 0.04 + 10/lepvec[itlep]->pT());
            dR=jetmom.deltaR_eta(lepmom);

            if(fabs(dR) <= DeltaRMax) overlap=true;
          }
          if(overlap) continue;
          Survivors.push_back(lepvec.at(itlep));
        }
        lepvec=Survivors;

        return;
      }


      void analyze(const HEPUtils::Event* event) {
        HEPUtilsAnalysis::analyze(event);

        // Get the missing energy and momentum in the event
        HEPUtils::P4 metVec = event->missingmom();
        double met = event->met();

        // Now define vectors of baseline objects, including:
        // - retrieval of electron, muon and jets from the event
        // - application of basic pT and eta cuts
        vector<HEPUtils::Particle*> electrons;
        for (HEPUtils::Particle* electron : event->electrons()) {
          if (electron->pT() > 10.
              && fabs(electron->eta()) < 2.47)
            electrons.push_back(electron);
        }
        vector<HEPUtils::Particle*> muons;
        for (HEPUtils::Particle* muon : event->muons()) {
          if (muon->pT() > 10.
              && fabs(muon->eta()) < 2.5)
            muons.push_back(muon);
        }

        vector<HEPUtils::Jet*> candJets;
        for (HEPUtils::Jet* jet : event->jets()) {
          if (jet->pT() > 25. && fabs(jet->eta()) < 2.5)
            candJets.push_back(jet);
        }

   	    // Jets
        vector<HEPUtils::Jet*> bJets;
        vector<HEPUtils::Jet*> nonBJets;

        // Find b-jets
        // TODO: implement misstag
        const std::vector<double>  a = {0,10.};
        const std::vector<double>  b = {0,10000.};
        const std::vector<double> c = {0.70}; // set b-tag efficiency to 70%, light jet mistag 1/381
        HEPUtils::BinnedFn2D<double> _eff2d(a,b,c);
        for (HEPUtils::Jet* jet : candJets) {
          bool hasTag=has_tag(_eff2d, jet->eta(), jet->pT());
          if(jet->btag() && hasTag){
            bJets.push_back(jet);
          } else {
            nonBJets.push_back(jet);
          }
        }

        // TODO: Check overlap removal
//        JetLeptonOverlapRemoval(nonBJets,electrons,0.2);
//        LeptonJetOverlapRemoval(electrons,nonBJets);
//        LeptonJetOverlapRemoval(electrons,bJets);
//        JetLeptonOverlapRemoval(nonBJets,muons,0.2);
//        SpecialLeptonJetOverlapRemoval(muons,nonBJets);
//        SpecialLeptonJetOverlapRemoval(muons,bJets);

        size_t nbJets = bJets.size();
        size_t nMuons=muons.size();
        size_t nElectrons=electrons.size();
        size_t nLeptons = nElectrons+nMuons;

        // Effective mass
        double meff = met;
        for(HEPUtils::Jet* jet : candJets){
          meff += jet->pT();
        }
        
	// Now order the jet collections by pT
	//std::sort(signalJets35.begin(), signalJets35.end(), sortByPT);
	//std::sort(signalBJets35.begin(), signalBJets35.end(), sortByPT);
	//std::sort(signalJets20.begin(), signalJets20.end(), sortByPT);
	//std::sort(signalBJets20.begin(), signalBJets20.end(), sortByPT);
		    
        // Increment cutFlowVector elements
        cutFlowVector_str[0]  = "No cuts ";
        cutFlowVector_str[1]  = "Trigger, 4 jets ($p_T > 40$ GeV, 2 b-tags)";
	    cutFlowVector_str[2]  = "$\\ge 4$ b-tags";
        cutFlowVector_str[3]  = "$\\ge 2$ Higgses ";
        cutFlowVector_str[4]  = "Lepton veto";
        cutFlowVector_str[5]  = "$X_{Wt} > 1.8$";
        cutFlowVector_str[6]  = "$X_{hh}^{SR} < 1.6$";
        cutFlowVector_str[7]  = "low-SR-MET0meff440";
        cutFlowVector_str[8]  = "low-SR-MET150meff440";

        // Cut flow from paper
        cutFlowVectorATLAS_130[0] = 169015.8;
        cutFlowVectorATLAS_130[1] =  11206.7;
        cutFlowVectorATLAS_130[2] =   1250.8;
        cutFlowVectorATLAS_130[3] =   1015.9;
        cutFlowVectorATLAS_130[4] =   1015.9;
        cutFlowVectorATLAS_130[5] =    961.9;
        cutFlowVectorATLAS_130[6] =    559.8;
        cutFlowVectorATLAS_130[7] =    217.4;
        cutFlowVectorATLAS_130[8] =      0.0;
 
        // Apply cuts to each signal region
        for(size_t j=0;j<NCUTS;j++){
          if(
             (j==0) ||

             (j==1 && nbJets > 1) ||

             (j==2 && nbJets > 3) ||

             (j==3 && nbJets > 3) ||
	     
             (j==4 && nbJets > 3 && nLeptons == 0) ||

             (j==5 && nbJets > 3 && nLeptons == 0) ||

             (j==6 && nbJets > 3 && nLeptons == 0) ||

             (j==7 && nbJets > 3 && nLeptons == 0 && meff > 440.) ||

             (j==8 && nbJets > 3 && nLeptons == 0 && meff > 440. && met > 150.)
	     
             ) cutFlowVector[j]++;
        }

        // Now increment signal region variables
        if(meff > 160. ) _meff160_ETmiss0++;
        if(met > 20. && meff > 160.) _meff160_ETmiss20++;

        return;
        
      } // End of analyze


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region number and total # events combination across threads
        HEPUtilsAnalysis::add(other);

        Analysis_ATLAS_13TeV_3b_24invfb* specificOther
          = dynamic_cast<Analysis_ATLAS_13TeV_3b_24invfb*>(other);

        // Here we will add the subclass member variables:
        if (NCUTS != specificOther->NCUTS) NCUTS = specificOther->NCUTS;
        for (size_t j=0; j<NCUTS; j++) {
          cutFlowVector[j] += specificOther->cutFlowVector[j];
          cutFlowVector_str[j] = specificOther->cutFlowVector_str[j];
        }
        _meff160_ETmiss0   += specificOther->_meff160_ETmiss0;
        _meff160_ETmiss20  += specificOther->_meff160_ETmiss20;

      }


      void collect_results() {

        // DEBUG
        double L = 24.3;
        double xsec = 6955.;
        cout << "DEBUG:" << endl;
        for (size_t i=0; i<NCUTS; i++)
        {
          double ATLAS_abs = cutFlowVectorATLAS_130[i];
        
          double eff = (double)cutFlowVector[i] / (double)cutFlowVector[0];
        
          double GAMBIT_scaled = eff * xsec * L;
        
          double ratio = GAMBIT_scaled/ATLAS_abs;
          cout << "DEBUG 1: i: " << i << ":   " << setprecision(4) << ATLAS_abs << "\t" << GAMBIT_scaled << "\t" << "\t" << ratio << "\t\t" << cutFlowVector_str[i] << endl;
        }
        cout << "DEBUG:" << endl;
        
        // Now fill a results object with the results for each SR
        SignalRegionData results_meff160_ETmiss0;
        results_meff160_ETmiss0.sr_label = "meff160_ETmiss0"; // label must be unique for each signal region
        results_meff160_ETmiss0.n_observed = 20.;             // set number of observed events (in LHC paper)
        results_meff160_ETmiss0.n_background = 16.21;         // set number of predicted background events (in LHC paper)
        results_meff160_ETmiss0.background_sys = 0.11;        // set background uncertainty (in LHC paper)
        results_meff160_ETmiss0.signal_sys = 0.;              // set signal uncertainty
        results_meff160_ETmiss0.n_signal = _meff160_ETmiss0;  // set this to number of signal events incremented in the analysis above
        add_result(results_meff160_ETmiss0);

        SignalRegionData results_meff160_ETmiss20;
        results_meff160_ETmiss20.sr_label = "meff160_ETmiss20";
        results_meff160_ETmiss20.n_observed = 3.;
        results_meff160_ETmiss20.n_background = 0.6503;
        results_meff160_ETmiss20.background_sys = 0.0747;
        results_meff160_ETmiss20.signal_sys = 0.;
        results_meff160_ETmiss20.n_signal = _meff160_ETmiss20;
        add_result(results_meff160_ETmiss20);

        return;
      }

      void clear() {

        _meff160_ETmiss20 = 0; _meff160_ETmiss20 = 0;
	
        std::fill(cutFlowVector.begin(), cutFlowVector.end(), 0);
      }



    };

    DEFINE_ANALYSIS_FACTORY(ATLAS_13TeV_3b_24invfb)

  }
}
