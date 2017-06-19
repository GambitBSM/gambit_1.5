///
///  \author Rose Kudzman-Blais
///  \date 2017 May
///
///  *********************************************

#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>
#include <fstream>
#include <ctime>

#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "gambit/ColliderBit/analyses/Perf_Plot.hpp"

using namespace std;

namespace Gambit {
  namespace ColliderBit {

    class Analysis_ATLAS_1LEPbb_20invfb : public HEPUtilsAnalysis {
    private:

      // Numbers passing cuts
      double _numSRA, _numSRB; 
      vector<int> cutFlowVector;
      vector<string> cutFlowVector_str;
      size_t NCUTS;

      Perf_Plot* plots;	
      ofstream cutflowFile;
      string analysisRunName;

    public:

      struct particleComparison {
  	bool operator() (HEPUtils::Particle* i,HEPUtils::Particle* j) {return (i->pT()>j->pT());}
      }	compareParticlePt;

      struct jetComparison {
  	bool operator() (HEPUtils::Jet* i,HEPUtils::Jet* j) {return (i->pT()>j->pT());}
      }	compareJetPt;

      Analysis_ATLAS_1LEPbb_20invfb() {

        _numSRA=0;
        _numSRB=0;

        NCUTS=8;
        set_luminosity(20.3);

        for (size_t i=0;i<NCUTS;i++){
          cutFlowVector.push_back(0);
          cutFlowVector_str.push_back("");
        }

	//time_t now = time(0);
	//tm *ltm = localtime(&now);
	analysisRunName = "ATLAS_1LEPbb_20invfb_130_0";
	//analysisRunName.append(to_string(ltm->tm_sec));
	vector<const char*> variables = {"met","mct","mbb","mt","j0pt","lpt"};
	plots = new Perf_Plot(analysisRunName+"_mbb", &variables);

      }


      void analyze(const HEPUtils::Event* event) {
	  
	HEPUtilsAnalysis::analyze(event);

        // Missing energy
        double met = event->met();

        // Baseline objects
        vector<HEPUtils::Particle*> baselineElectrons;
        for (HEPUtils::Particle* electron : event->electrons()) {
          if (electron->pT() > 10. && fabs(electron->eta()) < 2.47) baselineElectrons.push_back(electron);
        }
        vector<HEPUtils::Particle*> baselineMuons;
        for (HEPUtils::Particle* muon : event->muons()) {
          if (muon->pT() > 10. && fabs(muon->eta()) < 2.4) baselineMuons.push_back(muon);
        }
        vector<HEPUtils::Jet*> baselineJets;
        for (HEPUtils::Jet* jet : event->jets()) {
          if (jet->pT() > 20. && fabs(jet->eta()) < 4.5) baselineJets.push_back(jet);
        }

        //Overlap procedure
        vector<HEPUtils::Particle*> overlapElectrons;
	vector<HEPUtils::Particle*> overlapMuons;
	vector<HEPUtils::Jet*> overlapJets;

        //Remove any jet within dR=0.2 of an electrons
        for (size_t iJet=0;iJet<baselineJets.size();iJet++) {
          bool overlap=false;
          HEPUtils::P4 jetVec=baselineJets.at(iJet)->mom();
          for (size_t iEl=0;iEl<baselineElectrons.size();iEl++) {
            HEPUtils::P4 elVec=baselineElectrons.at(iEl)->mom();
            if (fabs(elVec.deltaR_eta(jetVec))<0.2)overlap=true;
          }
          if (!overlap)overlapJets.push_back(baselineJets.at(iJet));
        }

        //Remove electrons with dR=0.4 of surviving jets
        for (size_t iEl=0;iEl<baselineElectrons.size();iEl++) {
          bool overlap=false;
          HEPUtils::P4 elVec=baselineElectrons.at(iEl)->mom();
          for (size_t iJet=0;iJet<overlapJets.size();iJet++) {
            HEPUtils::P4 jetVec=overlapJets.at(iJet)->mom();
            if (fabs(elVec.deltaR_eta(jetVec))<0.4)overlap=true;
          }
          if (!overlap)overlapElectrons.push_back(baselineElectrons.at(iEl));
        }

        //Remove muons with dR=0.4 of surviving jets
        for (size_t iMu=0;iMu<baselineMuons.size();iMu++) {
          bool overlap=false;
          HEPUtils::P4 muVec=baselineMuons.at(iMu)->mom();
          for (size_t iJet=0;iJet<overlapJets.size();iJet++) {
            HEPUtils::P4 jetVec=overlapJets.at(iJet)->mom();
            if (fabs(muVec.deltaR_eta(jetVec))<0.4)overlap=true;
          }
          if (!overlap)overlapMuons.push_back(baselineMuons.at(iMu));
        }

	//Reject events with muons and electrons within dR=0.1
	for (size_t iEl=0;iEl<overlapElectrons.size();iEl++) {
	  HEPUtils::P4 elVec=overlapElectrons.at(iEl)->mom();
	  for (size_t iMu=0;iMu<overlapMuons.size();iMu++) {
	    HEPUtils::P4 muVec=overlapMuons.at(iMu)->mom();
	    if(fabs(elVec.deltaR_eta(muVec))<0.1) {
		overlapElectrons.clear();
		overlapMuons.clear();
		overlapJets.clear();
	    }
	  }
	}

  	//Reject events with muons within dR=0.05 of each other
        for (size_t iMu1=0;iMu1<overlapMuons.size();iMu1++) {
          HEPUtils::P4 muVec1=overlapMuons.at(iMu1)->mom();
          for (size_t iMu2=0;iMu2<overlapMuons.size();iMu2++) {
            HEPUtils::P4 muVec2=overlapMuons.at(iMu2)->mom();
            if(fabs(muVec1.deltaR_eta(muVec2))<0.05) {
                overlapElectrons.clear();
                overlapMuons.clear();
                overlapJets.clear();
            }
          }
        }


        // Signal requirements
        vector<HEPUtils::Particle*> signalLeptons;
        vector<HEPUtils::Jet*> signalJets;   
	vector<HEPUtils::Jet*> signalBJets;

	// Electrons
	for (size_t iEl=0;iEl<overlapElectrons.size();iEl++) {
	  if (overlapElectrons.at(iEl)->pT() > 25. && fabs(overlapElectrons.at(iEl)->eta()) < 2.47)signalLeptons.push_back(overlapElectrons.at(iEl));
        }
        
        //Muons
        for (size_t iMu=0;iMu<overlapMuons.size();iMu++) {
          if (overlapMuons.at(iMu)->pT() > 25. && fabs(overlapMuons.at(iMu)->eta()) < 2.40)signalLeptons.push_back(overlapMuons.at(iMu)); 
        } 
	       
        //Jets
        for (size_t iJet=0;iJet<overlapJets.size();iJet++) {
          if (overlapJets.at(iJet)->pT() > 25. && fabs(overlapJets.at(iJet)->eta()) < 2.40) {
	    signalJets.push_back(overlapJets.at(iJet));   
	    if (overlapJets.at(iJet)->btag())signalBJets.push_back(overlapJets.at(iJet));             
          }
	}

        //Variable definitions
        int nSignalLeptons = signalLeptons.size();
        int nBaselineLeptons = overlapElectrons.size() + overlapMuons.size();
        int nSignalJets = signalJets.size();
        int nSignalBJets = signalBJets.size();
	sort(signalJets.begin(), signalJets.end(), compareJetPt);
	sort(signalLeptons.begin(), signalLeptons.end(), compareParticlePt);

        //Preselection
        bool leadingBJets = isLeadingBJets(signalJets, signalBJets);

        bool preselection = 0; 
        if (nSignalLeptons == 1 && nBaselineLeptons == 1) {
	  if (nSignalBJets == 1 || nSignalBJets == 2) {
            if (leadingBJets) { 
              preselection = 1;
            }
          }
	}
	
        //Signal regions
	double mT=0; 
	if (nSignalLeptons) {
          mT = sqrt(2*signalLeptons.at(0)->pT()*met*(1-cos(signalLeptons.at(0)->phi()-event->missingmom().phi())));
        }
	
        double mCT=0;
        double mbb=0;
	if (nSignalJets>1) {
          mCT = sqrt(2*signalJets.at(0)->pT()*signalJets.at(1)->pT()*(1+cos(signalJets.at(0)->phi()-signalJets.at(1)->phi())));
          mbb = (signalJets.at(0)->mom()+signalJets.at(1)->mom()).m(); 
	}	

	bool SRA=false;
	bool SRB=false;
        if (nSignalBJets == 2 && (nSignalJets == 2 || nSignalJets == 3) && preselection) {
          if (met > 100. && mCT > 160. && mbb > 105. && mbb < 135.) {
            //SRA
            if (mT > 100. && mT < 130.) {
              _numSRA++;
	      SRA=true;
            }
            //SRB
            if (mT > 130.) {
              _numSRB++;
	      SRB=true;   
            }
          }
        }                      
	
	if (preselection &&  met>50. && mT>40. && mbb>40. && nSignalBJets==2 && met > 100. && mCT > 160. && mT>100. && mbb > 45. && mbb < 195.) {
	  vector<double> variables = {met, mCT, mbb, mT, signalJets.at(0)->pT(), signalLeptons.at(0)->pT()};
	  plots->fill(&variables);
	}

        cutFlowVector_str[0] = "No cuts ";
        cutFlowVector_str[1] = "1 lepton + 2 bjets";
        cutFlowVector_str[2] = "met > 100 GeV";
        cutFlowVector_str[3] = "mCT > 160 GeV";
        cutFlowVector_str[4] = "mT > 100 GeV";
        cutFlowVector_str[5] = "45 < mbb < 195 GeV";
        cutFlowVector_str[6] = "SRA";
        cutFlowVector_str[7] = "SRB";

        for (size_t j=0;j<NCUTS;j++){
          if(
             (j==0) ||
             
	     (j==1 && nSignalLeptons == 1 && nSignalBJets == 2) ||
             
             (j==2 && nSignalLeptons == 1 && nSignalBJets == 2 && met > 100.) ||

             (j==3 && nSignalLeptons == 1 && nSignalBJets == 2 && met > 100. && mCT > 160.) ||
 
             (j==4 && nSignalLeptons == 1 && nSignalBJets == 2 && met > 100. && mCT > 160. && mT>100.) ||
             
	     (j==5 && nSignalLeptons == 1 && nSignalBJets == 2 && met > 100. && mCT > 160. && mT>100. && mbb > 45. && mbb < 195.) ||  

             (j==6 && nSignalLeptons == 1 && nSignalBJets == 2 && met > 100. && mCT > 160. && mT>100. && mbb > 45. && mbb < 195. && SRA) ||  
            
	     (j==7 && nSignalLeptons == 1 && nSignalBJets == 2 && met > 100. && mCT > 160. && mT>100. && mbb > 45. && mbb < 195. && SRB) ){

            cutFlowVector[j]++;

          }
        }
      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        HEPUtilsAnalysis::add(other);

        Analysis_ATLAS_1LEPbb_20invfb* specificOther
                = dynamic_cast<Analysis_ATLAS_1LEPbb_20invfb*>(other);

        // Here we will add the subclass member variables:
        if (NCUTS != specificOther->NCUTS) NCUTS = specificOther->NCUTS;
        for (size_t j = 0; j < NCUTS; j++) {
          cutFlowVector[j] += specificOther->cutFlowVector[j];
          cutFlowVector_str[j] = specificOther->cutFlowVector_str[j];
        }
        _numSRA += specificOther->_numSRA;
        _numSRB += specificOther->_numSRB;
      }


      void collect_results() {

	string path = "ColliderBit/results/cutflow_";
	path.append(analysisRunName);
	path.append(".txt");
	cutflowFile.open(path.c_str());
	cutflowFile<<"XSEC_PER_EVENT: "<<xsec_per_event()<<endl;
	cutflowFile<<"XSEC: "<<xsec()<<endl;	
	cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
        cutflowFile << "CUT FLOW: ATLAS 1 lepton, 2 bjets paper "<<endl;
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------"<<endl;

        cutflowFile<< right << setw(60) << "CUT" << setw(20) << "RAW" << setw(20) << "SCALED" << setw(20) << " % " << endl;
        for (size_t j=0; j<NCUTS; j++) {
          cutflowFile << right << setw(60) << cutFlowVector_str[j].c_str() << setw(20) << cutFlowVector[j] << setw(20) << cutFlowVector[j]*xsec_per_event()*1000.*luminosity() << setw(20) << 100.*cutFlowVector[j]/cutFlowVector[0] << endl;
        }
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
	cutflowFile.close();

	plots->createFile(luminosity(),xsec_per_event());

        SignalRegionData results_SRA;
        results_SRA.analysis_name = "Analysis_ATLAS_1LEPbb_20invfb";
        results_SRA.sr_label = "SRA";
        results_SRA.n_observed = 4.;
        results_SRA.n_background = 5.69; 
        results_SRA.background_sys = 1.10;
        results_SRA.signal_sys = 0.; 
        results_SRA.n_signal = _numSRA;
        add_result(results_SRA);

        SignalRegionData results_SRB;
        results_SRB.analysis_name = "Analysis_ATLAS_1LEPbb_20invfb";
        results_SRB.sr_label = "SRB";
        results_SRB.n_observed = 3.;
        results_SRB.n_background = 2.67; 
        results_SRB.background_sys = 0.69;
        results_SRB.signal_sys = 0.; 
        results_SRB.n_signal = _numSRB;
        add_result(results_SRB);

      }

      bool isLeadingBJets(vector<HEPUtils::Jet*> jets, vector<HEPUtils::Jet*> bjets) {
        sort(jets.begin(), jets.end(), compareJetPt);
	sort(bjets.begin(), bjets.end(), compareJetPt);
	int nbjet = bjets.size();
	jets.resize(nbjet);
        return (jets == bjets);
      }

    };


    // Factory fn
    DEFINE_ANALYSIS_FACTORY(ATLAS_1LEPbb_20invfb)


  }
}
