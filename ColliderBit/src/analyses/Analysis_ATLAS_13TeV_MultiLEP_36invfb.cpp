///
///  \author Rose Kudzman-Blais
///  \date 2017 May
///
///  *********************************************


#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>
#include <algorithm>

#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "gambit/ColliderBit/mt2_bisect.h"
#include "gambit/ColliderBit/analyses/Perf_Plot.hpp"

using namespace std;

namespace Gambit {
  namespace ColliderBit {

    class Analysis_ATLAS_13TeV_MultiLEP_36invfb : public HEPUtilsAnalysis {
    private:

      // Numbers passing cuts
      double _numSR2_SF_loose, _numSR2_SF_tight, _numSR2_DF_100, _numSR2_DF_150, _numSR2_DF_200, _numSR2_DF_300, _numSR2_int, _numSR2_high, _numSR2_low, _numSR3_slep_a, _numSR3_slep_b, _numSR3_slep_c, _numSR3_slep_d, _numSR3_slep_e, _numSR3_WZ_0Ja, _numSR3_WZ_0Jb, _numSR3_WZ_0Jc, _numSR3_WZ_1Ja, _numSR3_WZ_1Jb, _numSR3_WZ_1Jc; 
      vector<int> cutFlowVector;
      vector<string> cutFlowVector_str;
      size_t NCUTS;

      Perf_Plot* plots;
      ofstream cutflowFile;
      string analysisRunName;


    public:

      Analysis_ATLAS_13TeV_MultiLEP_36invfb() {

	_numSR2_SF_loose=0;
	_numSR2_SF_tight=0;
	_numSR2_DF_100=0;
	_numSR2_DF_150=0;
	_numSR2_DF_200=0;
	_numSR2_DF_300=0;
	_numSR2_int=0;
	_numSR2_high=0;
	_numSR2_low=0;
	_numSR3_slep_a=0;
	_numSR3_slep_b=0;
	_numSR3_slep_c=0;
	_numSR3_slep_d=0;
	_numSR3_slep_e=0;
	_numSR3_WZ_0Ja=0;
	_numSR3_WZ_0Jb=0;
	_numSR3_WZ_0Jc=0;
	_numSR3_WZ_1Ja=0;
	_numSR3_WZ_1Jb=0;
	_numSR3_WZ_1Jc=0;


	NCUTS=9;
	set_luminosity(36.1);

        for (size_t i=0;i<NCUTS;i++){
          cutFlowVector.push_back(0);
          cutFlowVector_str.push_back("");
        }

        analysisRunName = "ATLAS_13TeV_MultiLEP_36invfb_test";
        vector<const char*> variables = {"met"};
        plots = new Perf_Plot(analysisRunName, &variables);

      }

      struct ptComparison {
        bool operator() (HEPUtils::Particle* i,HEPUtils::Particle* j) {return (i->pT()>j->pT());}
      } comparePt;
      
      struct ptJetComparison {
        bool operator() (HEPUtils::Jet* i,HEPUtils::Jet* j) {return (i->pT()>j->pT());}
      } compareJetPt;

      void analyze(const HEPUtils::Event* event) {
	HEPUtilsAnalysis::analyze(event);
        double met = event->met();

        // Baseline objects
        vector<HEPUtils::Particle*> baselineLeptons;
        vector<HEPUtils::Particle*> baselineElectrons;
        for (HEPUtils::Particle* electron : event->electrons()) {
	  if (electron->pT()>10. &&fabs(electron->eta())<2.47) {
            baselineElectrons.push_back(electron);
	    baselineLeptons.push_back(electron);
          }
	}

        vector<HEPUtils::Particle*> baselineMuons;
        for (HEPUtils::Particle* muon : event->muons()) {
	  if (muon->pT()>10. &&fabs(muon->eta())<2.4) {
            baselineMuons.push_back(muon);
	    baselineLeptons.push_back(muon);
          }
	}

        vector<HEPUtils::Jet*> baselineJets;
        for (HEPUtils::Jet* jet : event->jets()) {
          if (jet->pT()>20. &&fabs(jet->eta())<2.8) {
            baselineJets.push_back(jet);
	  }
        }

	//Overlap Removal + Signal Objects	
	vector<HEPUtils::Particle*> signalElectrons;
	vector<HEPUtils::Particle*> signalMuons;
	vector<HEPUtils::Particle*> signalLeptons;
	vector<HEPUtils::Jet*> signalJets;
	vector<HEPUtils::Jet*> signalBJets;

        const vector<double>  a = {0,10.};
        const vector<double>  b = {0,10000.};
        const vector<double> c = {0.77};
        HEPUtils::BinnedFn2D<double> _eff2d(a,b,c);
	vector<HEPUtils::Jet*> overlapJet;
        for (size_t iJet=0;iJet<baselineJets.size();iJet++) {
	  vector<HEPUtils::Particle*> overlapEl;
          bool hasTag=has_tag(_eff2d, baselineJets.at(iJet)->eta(), baselineJets.at(iJet)->pT());
          for (size_t iEl=0;iEl<baselineElectrons.size();iEl++) {
            if (fabs(baselineElectrons.at(iEl)->mom().deltaR_eta(baselineJets.at(iJet)->mom()))<0.2)overlapEl.push_back(baselineElectrons.at(iEl));
	  }
	  if (overlapEl.size()>0 && (baselineJets.at(iJet)->btag() && hasTag)) {
	    for (size_t iO=0;iO<overlapEl.size();iO++) {
	      baselineElectrons.erase(remove(baselineElectrons.begin(), baselineElectrons.end(), overlapEl.at(iO)), baselineElectrons.end());
	    }
	  }
	  if (overlapEl.size()>0 && !(baselineJets.at(iJet)->btag() && hasTag))overlapJet.push_back(baselineJets.at(iJet));
	}
	for (size_t iO=0;iO<overlapJet.size();iO++) {
	  baselineJets.erase(remove(baselineJets.begin(), baselineJets.end(), overlapJet.at(iO)), baselineJets.end());
	}

	for (size_t iEl=0;iEl<baselineElectrons.size();iEl++) {
	  bool overlap=false;
	  for (size_t iJet=0;iJet<baselineJets.size();iJet++) {
            if (fabs(baselineElectrons.at(iEl)->mom().deltaR_eta(baselineJets.at(iJet)->mom()))<0.4)overlap=true;
	    
	  }
	  if (!overlap) {
	    signalElectrons.push_back(baselineElectrons.at(iEl));
	    signalLeptons.push_back(baselineElectrons.at(iEl));
	  }
	}

	for (size_t iJet=0;iJet<baselineJets.size();iJet++) {
	  bool overlap=false;
	  for (size_t iMu=0;iMu<baselineMuons.size();iMu++) {
            if (fabs(baselineMuons.at(iMu)->mom().deltaR_eta(baselineJets.at(iJet)->mom()))<0.2 && baselineMuons.at(iMu)->pT()>0.7*baselineJets.at(iJet)->pT())overlap=true;
	    
	  }
	  if (!overlap) {
            bool hasTag=has_tag(_eff2d, baselineJets.at(iJet)->eta(), baselineJets.at(iJet)->pT());
	    signalJets.push_back(baselineJets.at(iJet));
	    if (baselineJets.at(iJet)->btag() && hasTag) {
	      if (fabs(baselineJets.at(iJet)->eta())<2.5 || (fabs(baselineJets.at(iJet)->eta())>2.4 && baselineJets.at(iJet)->pT()>50))signalBJets.push_back(baselineJets.at(iJet));
	    }
	  }
	}

	for (size_t iMu=0;iMu<baselineMuons.size();iMu++) {
	  bool overlap=false;
	  for (size_t iJet=0;iJet<signalJets.size();iJet++) {
            if (fabs(baselineMuons.at(iMu)->mom().deltaR_eta(signalJets.at(iJet)->mom()))<0.4)overlap=true;
	    
	  }
	  if (!overlap) {
	    signalMuons.push_back(baselineMuons.at(iMu));
	    signalLeptons.push_back(baselineMuons.at(iMu));
	  }
	}

	sort(signalJets.begin(),signalJets.end(),compareJetPt);	
        sort(signalLeptons.begin(),signalLeptons.end(),comparePt);
        size_t nSignalLeptons = signalLeptons.size();
        size_t nSignalJets = signalJets.size();
	size_t nSignalBJets = signalBJets.size();
	size_t nSignalNonBJets = nSignalJets-nSignalBJets;       
 
	//Preselection
        bool preselection=false; 
	if ((nSignalLeptons==2 || nSignalLeptons==3) && baselineLeptons.size()==nSignalLeptons) {
	  if (signalLeptons.at(0)->pT()>25 && signalLeptons.at(1)->pT()>20)preselection=true;
	}

        //Signal regions
	vector<vector<HEPUtils::Particle*>> OSSFpairs=getOSSFpair(signalLeptons);
	vector<vector<HEPUtils::Particle*>> OSpairs=getOSpair(signalLeptons);

	double mll=0;
	double mT2=0;
	double mjj=0;
	if (nSignalLeptons>1) {
	  mll=(signalLeptons.at(0)->mom()+signalLeptons.at(1)->mom()).m();

          double pLep1[3] = {signalLeptons.at(0)->mass(), signalLeptons.at(0)->mom().px(), signalLeptons.at(0)->mom().py()};
          double pLep2[3] = {signalLeptons.at(1)->mass(), signalLeptons.at(1)->mom().px(), signalLeptons.at(1)->mom().py()};
          double pMiss[3] = {0., event->missingmom().px(), event->missingmom().py() };
          double mn = 0.;

          mt2_bisect::mt2 mt2_calc;
          mt2_calc.set_momenta(pLep1,pLep2,pMiss);
          mt2_calc.set_mn(mn);
          mT2 = mt2_calc.get_mt2();
	}
	if (nSignalJets>1)mjj=(signalJets.at(0)->mom()+signalJets.at(1)->mom()).m();

	//Large jet veto
	bool large_jet_veto=true;
	for (size_t iJet=0;iJet<nSignalJets;iJet++) {
	  if (signalJets.at(iJet)->pT()>60 && fabs(signalJets.at(iJet)->eta())<2.4 && signalJets.at(iJet)->btag()==false)large_jet_veto=false;
	}

	//Bjet veto
	bool bjet_veto=true;
	for (size_t iJet=0;iJet<nSignalBJets;iJet++) {
	  if (signalBJets.at(iJet)->pT()>20 && fabs(signalBJets.at(iJet)->eta())<2.4)bjet_veto=false;
	}

	//2lep+0jet
	if (preselection && nSignalLeptons==2 && OSpairs.size()==1 && mll>40 && large_jet_veto && bjet_veto) {
	  if (OSSFpairs.size()==1) {
	    if (mT2>100 && mll>111)_numSR2_SF_loose++;
	    if (mT2>130 && mll>300)_numSR2_SF_tight++; 
	  }
	  if (OSSFpairs.size()==0) {
	    if (mT2>100)_numSR2_DF_100++;
	    if (mT2>150)_numSR2_DF_150++;
	    if (mT2>200)_numSR2_DF_200++;
	    if (mT2>300)_numSR2_DF_300++;
	  }
	}
	
	//2lep+jets
	if (preselection && nSignalLeptons==2 && OSpairs.size() && mll>40 && bjet_veto) {
	  if (signalLeptons.at(1)->pT()>25 && nSignalJets>1) {
	    if (signalJets.at(0)->pT()>30 && signalJets.at(1)->pT()>30) {
		HEPUtils::P4 Z=signalLeptons.at(0)->mom()+signalLeptons.at(1)->mom();
	     	double deltaR_ll=fabs(signalLeptons.at(0)->mom().deltaR_eta(signalLeptons.at(1)->mom()));
		double deltaR_jj=fabs(signalJets.at(0)->mom().deltaR_eta(signalJets.at(1)->mom()));
	      //SR2_int + SR2_high
	      if (nSignalNonBJets>1 && mll>81. && mll<101. && mjj>70. && mjj<100.) {
		HEPUtils::P4 W=signalJets.at(0)->mom()+signalJets.at(1)->mom();
	        double deltaPhi_met_W = fabs(W.phi()-event->missingmom().phi());
		if (Z.pT()>80. && W.pT()>100. && mT2>100. && deltaR_jj<1.5 && deltaR_ll<1.8 && deltaPhi_met_W>0.5 && deltaPhi_met_W<3.0) {
		  if (met>150)_numSR2_int++;
		  if (met>250)_numSR2_high++;
		}
	      }
	      //SR2_low_2J
	      if (nSignalNonBJets==2 && mll>81. && mll<101. && mjj>70. && mjj<90. && met>100. && Z.pT()>60.) {
		HEPUtils::P4 W=signalJets.at(0)->mom()+signalJets.at(1)->mom();
	        double deltaPhi_met_W=fabs(W.phi()-event->missingmom().phi());
	        double deltaPhi_met_Z=fabs(Z.phi()-event->missingmom().phi());
		if (deltaPhi_met_Z<0.8 && deltaPhi_met_W>1.5 && (met/Z.pT())>0.6 && (met/Z.pT())<1.6 && (met/W.pT())<0.8)_numSR2_low++;
	      }	
	      //SR2_low_3J
	      if (nSignalNonBJets>2 && nSignalNonBJets<6 && mll>86 && mll<96 && mjj>70 && mjj<90 && met>100 && Z.pT()>40 && deltaR_jj<2.2) {
		vector<HEPUtils::P4> W_ISR=get_W_ISR(signalJets,Z,event->missingmom());
		double deltaPhi_met_W=fabs(W_ISR.at(0).phi()-event->missingmom().phi());
		double deltaPhi_met_ISR=fabs(W_ISR.at(1).phi()-event->missingmom().phi());
		double deltaPhi_met_jet1=fabs(signalJets.at(0)->phi()-event->missingmom().phi());
		if (deltaPhi_met_W<2.2 && deltaPhi_met_ISR>2.4 && deltaPhi_met_jet1>2.6 && (met/W_ISR.at(1).pT())>0.4 && (met/W_ISR.at(1).pT())<0.8 && fabs(Z.eta())<1.6 && signalJets.at(2)->pT()>30.)_numSR2_low++;
	      }

	    }	
	  }
	}
	
	//3lep
	if (preselection && nSignalLeptons==3 && bjet_veto && OSSFpairs.size()) {
	  double mTmin=999;
	  double mSFOS=999;
	  for (size_t iPa=0;iPa<OSSFpairs.size();iPa++) {
	    for (size_t iLep=0;iLep<signalLeptons.size();iLep++) {
	      if (signalLeptons.at(iLep)!=OSSFpairs.at(iPa).at(0) && signalLeptons.at(iLep)!=OSSFpairs.at(iPa).at(1)) {
	        double mT = sqrt(2*signalLeptons.at(iLep)->pT()*met*(1-cos(signalLeptons.at(iLep)->phi()-event->missingmom().phi())));		
	        if (mT<mTmin) { 
		  mTmin=mT;
		  mSFOS=(OSSFpairs.at(iPa).at(0)->mom()+OSSFpairs.at(iPa).at(1)->mom()).m();
		}	
	      }
	    }
	  }
	  if (mSFOS<81.2 && met>130. && mTmin>110.) {
	    if (signalLeptons.at(2)->pT()>20. && signalLeptons.at(2)->pT()<30.)_numSR3_slep_a++;
	    if (signalLeptons.at(2)->pT()>30.)_numSR3_slep_b++;
	  }
	  if (mSFOS>101.2 && met>130. && mTmin>110.) {
	    if (signalLeptons.at(2)->pT()>20. && signalLeptons.at(2)->pT()<50.)_numSR3_slep_c++;
	    if (signalLeptons.at(2)->pT()>50. && signalLeptons.at(2)->pT()<80.)_numSR3_slep_d++;
	    if (signalLeptons.at(2)->pT()>80.)_numSR3_slep_e++;
	  }
	  if (mSFOS>81.2 && mSFOS<101.2 && nSignalNonBJets==0 && mTmin>110.) {
	    if (met>60. && met<120.)_numSR3_WZ_0Ja++;
	    if (met>120. && met<170.)_numSR3_WZ_0Jb++;
	    if (met>170.)_numSR3_WZ_0Jc++;
	  }
	  if (mSFOS>81.2 && mSFOS<101.2 && nSignalNonBJets>0) {
	    if (met>120. && met<200. && mTmin>110. && (signalLeptons.at(0)->mom()+signalLeptons.at(1)->mom()+signalLeptons.at(2)->mom()).pT()<120. && signalJets.at(0)->pT()>70.)_numSR3_WZ_1Ja++;
	    if (met>200. && mTmin>110. && mTmin<160.)_numSR3_WZ_1Jb++;
	    if (met>200. && signalLeptons.at(2)->pT()>35. && mTmin>160.)_numSR3_WZ_1Jc++;
	  }
	}

        if (preselection) {
          vector<double> variables = {met};
          plots->fill(&variables);
        }

	bool lepton_pT_veto=false;
	if (nSignalLeptons==2) {
	  if (signalLeptons.at(0)->pT()>20. && signalLeptons.at(1)->pT()>20.)lepton_pT_veto=true;
	}

        cutFlowVector_str[0] = "All events";
        cutFlowVector_str[1] = "2 baseline leptons";
        cutFlowVector_str[2] = "2 signal leptons";
        cutFlowVector_str[3] = "Opposite-sign";
        cutFlowVector_str[4] = "Lepton pT > 20 GeV";
        cutFlowVector_str[5] = "mll > 40 GeV";
        cutFlowVector_str[6] = "b-veto";
        cutFlowVector_str[7] = "met > 50 GeV";
        cutFlowVector_str[8] = "mT2 > 75 GeV";

        for (size_t j=0;j<NCUTS;j++){
          if(
             (j==0) ||

	     (j==1 && baselineLeptons.size()==2) ||

	     (j==2 && baselineLeptons.size()==2 && nSignalLeptons==2) ||

	     (j==3 && baselineLeptons.size()==2 && nSignalLeptons==2 && OSpairs.size()) ||

	     (j==4 && baselineLeptons.size()==2 && nSignalLeptons==2 && OSpairs.size() && lepton_pT_veto) ||

	     (j==5 && baselineLeptons.size()==2 && nSignalLeptons==2 && OSpairs.size() && lepton_pT_veto && mll>40.) ||

	     (j==6 && baselineLeptons.size()==2 && nSignalLeptons==2 && OSpairs.size() && lepton_pT_veto && mll>40. && bjet_veto)  ||

	     (j==7 && baselineLeptons.size()==2 && nSignalLeptons==2 && OSpairs.size() && lepton_pT_veto && mll>40. && bjet_veto && met>50.) ||            
 
             (j==8 && baselineLeptons.size()==2 && nSignalLeptons==2 && OSpairs.size() && lepton_pT_veto && mll>40. && bjet_veto && met>50. && mT2>75.) )

	  cutFlowVector[j]++;
	}

      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        
	HEPUtilsAnalysis::add(other);

        Analysis_ATLAS_13TeV_MultiLEP_36invfb* specificOther
                = dynamic_cast<Analysis_ATLAS_13TeV_MultiLEP_36invfb*>(other);

        // Here we will add the subclass member variables:
        if (NCUTS != specificOther->NCUTS) NCUTS = specificOther->NCUTS;
        for (size_t j = 0; j < NCUTS; j++) {
          cutFlowVector[j] += specificOther->cutFlowVector[j];
          cutFlowVector_str[j] = specificOther->cutFlowVector_str[j];
        }

	_numSR2_SF_loose+= specificOther->_numSR2_SF_loose;
	_numSR2_SF_tight+= specificOther->_numSR2_SF_tight;
	_numSR2_DF_100+= specificOther->_numSR2_DF_100;
	_numSR2_DF_150+= specificOther->_numSR2_DF_150;
	_numSR2_DF_200+= specificOther->_numSR2_DF_200;
	_numSR2_DF_300+= specificOther->_numSR2_DF_300;
	_numSR2_int+= specificOther->_numSR2_int;
	_numSR2_high+= specificOther->_numSR2_high;
	_numSR2_low+= specificOther->_numSR2_low;
	_numSR3_slep_a+= specificOther->_numSR3_slep_a;
	_numSR3_slep_b+= specificOther->_numSR3_slep_b;
	_numSR3_slep_c+= specificOther->_numSR3_slep_c;
	_numSR3_slep_d+= specificOther->_numSR3_slep_d;
	_numSR3_slep_e+= specificOther->_numSR3_slep_e;
	_numSR3_WZ_0Ja+= specificOther->_numSR3_WZ_0Ja;
	_numSR3_WZ_0Jb+= specificOther->_numSR3_WZ_0Jb;
	_numSR3_WZ_0Jc+= specificOther->_numSR3_WZ_0Jc;
	_numSR3_WZ_1Ja+= specificOther->_numSR3_WZ_1Ja;
	_numSR3_WZ_1Jb+= specificOther->_numSR3_WZ_1Jb;
	_numSR3_WZ_1Jc+= specificOther->_numSR3_WZ_1Jc;
      }


      void collect_results() {

        string path = "ColliderBit/results/cutflow_";
        path.append(analysisRunName);
        path.append(".txt");
        cutflowFile.open(path.c_str());
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
        cutflowFile << "CUT FLOW: ATLAS Multi-lepton paper "<<endl;
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------"<<endl;
        cutflowFile<< right << setw(60) << "CUT" << setw(20) << "RAW" << setw(20) << " % " << endl;
        for (size_t j=0; j<NCUTS; j++) {
          cutflowFile << right << setw(60) << cutFlowVector_str[j].c_str() << setw(20) << cutFlowVector[j] << setw(20) << 100.*cutFlowVector[j]/cutFlowVector[0] << endl;
        }
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
        cutflowFile.close();

        plots->createFile();

        //Now fill a results object with the results for each SR
        SignalRegionData results_SR2_SF_loose;
        results_SR2_SF_loose.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR2_SF_loose.sr_label = "SR2_SF_loose";
        results_SR2_SF_loose.n_observed = 133.;
        results_SR2_SF_loose.n_background = 133.; 
        results_SR2_SF_loose.background_sys = 22.;
        results_SR2_SF_loose.signal_sys = 0.; 
        results_SR2_SF_loose.n_signal = _numSR2_SF_loose;
	add_result(results_SR2_SF_loose);

        SignalRegionData results_SR2_SF_tight;
        results_SR2_SF_tight.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR2_SF_tight.sr_label = "SR2_SF_tight";
        results_SR2_SF_tight.n_observed = 9.;
        results_SR2_SF_tight.n_background = 9.8; 
        results_SR2_SF_tight.background_sys = 2.9;
        results_SR2_SF_tight.signal_sys = 0.; 
        results_SR2_SF_tight.n_signal = _numSR2_SF_tight;
	add_result(results_SR2_SF_tight);

        SignalRegionData results_SR2_DF_100;
        results_SR2_DF_100.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR2_DF_100.sr_label = "SR2_DF_100";
        results_SR2_DF_100.n_observed = 78.;
        results_SR2_DF_100.n_background = 68.; 
        results_SR2_DF_100.background_sys = 7.;
        results_SR2_DF_100.signal_sys = 0.; 
        results_SR2_DF_100.n_signal = _numSR2_DF_100;
	add_result(results_SR2_DF_100);

        SignalRegionData results_SR2_DF_150;
        results_SR2_DF_150.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR2_DF_150.sr_label = "SR2_DF_150";
        results_SR2_DF_150.n_observed = 11;
        results_SR2_DF_150.n_background = 11.5; 
        results_SR2_DF_150.background_sys = 3.1;
        results_SR2_DF_150.signal_sys = 0.; 
        results_SR2_DF_150.n_signal = _numSR2_DF_150;
	add_result(results_SR2_DF_150);

        SignalRegionData results_SR2_DF_200;
        results_SR2_DF_200.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR2_DF_200.sr_label = "SR2_DF_200";
        results_SR2_DF_200.n_observed = 6.;
        results_SR2_DF_200.n_background = 2.1; 
        results_SR2_DF_200.background_sys = 1.9;
        results_SR2_DF_200.signal_sys = 0.; 
        results_SR2_DF_200.n_signal = _numSR2_DF_200;
	add_result(results_SR2_DF_200);

        SignalRegionData results_SR2_DF_300;
        results_SR2_DF_300.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR2_DF_300.sr_label = "SR2_DF_300";
        results_SR2_DF_300.n_observed = 2.;
        results_SR2_DF_300.n_background = 0.6; 
        results_SR2_DF_300.background_sys = 0.6;
        results_SR2_DF_300.signal_sys = 0.; 
        results_SR2_DF_300.n_signal = _numSR2_DF_300;
	add_result(results_SR2_DF_300);

        SignalRegionData results_SR2_int;
        results_SR2_int.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR2_int.sr_label = "SR2_int";
        results_SR2_int.n_observed = 2.;
        results_SR2_int.n_background = 4.1; 
        results_SR2_int.background_sys = 2.6;
        results_SR2_int.signal_sys = 0.; 
        results_SR2_int.n_signal = _numSR2_int;
	add_result(results_SR2_int);

        SignalRegionData results_SR2_high;
        results_SR2_high.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR2_high.sr_label = "SR2_high";
        results_SR2_high.n_observed = 0.;
        results_SR2_high.n_background = 1.6; 
        results_SR2_high.background_sys = 1.6;
        results_SR2_high.signal_sys = 0.; 
        results_SR2_high.n_signal = _numSR2_high;
	add_result(results_SR2_high);

        SignalRegionData results_SR2_low;
        results_SR2_low.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR2_low.sr_label = "SR2_low";
        results_SR2_low.n_observed = 11.;
        results_SR2_low.n_background = 4.2; 
        results_SR2_low.background_sys = 3.8;
        results_SR2_low.signal_sys = 0.; 
        results_SR2_low.n_signal = _numSR2_low;
	add_result(results_SR2_low);

        SignalRegionData results_SR3_slep_a;
        results_SR3_slep_a.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR3_slep_a.sr_label = "SR3_slep_a";
        results_SR3_slep_a.n_observed = 4.;
        results_SR3_slep_a.n_background = 2.23; 
        results_SR3_slep_a.background_sys = 0.79;
        results_SR3_slep_a.signal_sys = 0.; 
        results_SR3_slep_a.n_signal = _numSR3_slep_a;
	add_result(results_SR3_slep_a);

        SignalRegionData results_SR3_slep_b;
        results_SR3_slep_b.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR3_slep_b.sr_label = "SR3_slep_b";
        results_SR3_slep_b.n_observed = 3.;
        results_SR3_slep_b.n_background = 2.79; 
        results_SR3_slep_b.background_sys = 0.43;
        results_SR3_slep_b.signal_sys = 0.; 
        results_SR3_slep_b.n_signal = _numSR3_slep_b;
	add_result(results_SR3_slep_b);

        SignalRegionData results_SR3_slep_c;
        results_SR3_slep_c.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR3_slep_c.sr_label = "SR3_slep_c";
        results_SR3_slep_c.n_observed = 9.;
        results_SR3_slep_c.n_background = 5.42; 
        results_SR3_slep_c.background_sys = 0.93;
        results_SR3_slep_c.signal_sys = 0.; 
        results_SR3_slep_c.n_signal = _numSR3_slep_c;
	add_result(results_SR3_slep_c);

        SignalRegionData results_SR3_slep_d;
        results_SR3_slep_d.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR3_slep_d.sr_label = "SR3_slep_d";
        results_SR3_slep_d.n_observed = 0.;
        results_SR3_slep_d.n_background = 1.42; 
        results_SR3_slep_d.background_sys = 0.38;
        results_SR3_slep_d.signal_sys = 0.; 
        results_SR3_slep_d.n_signal = _numSR3_slep_d;
	add_result(results_SR3_slep_d);

        SignalRegionData results_SR3_slep_e;
        results_SR3_slep_e.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR3_slep_e.sr_label = "SR3_slep_e";
        results_SR3_slep_e.n_observed = 0.;
        results_SR3_slep_e.n_background = 1.14; 
        results_SR3_slep_e.background_sys = 0.23;
        results_SR3_slep_e.signal_sys = 0.; 
        results_SR3_slep_e.n_signal = _numSR3_slep_e;
	add_result(results_SR3_slep_e);

        SignalRegionData results_SR3_WZ_0Ja;
        results_SR3_WZ_0Ja.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR3_WZ_0Ja.sr_label = "SR3_WZ_0Ja";
        results_SR3_WZ_0Ja.n_observed = 21.;
        results_SR3_WZ_0Ja.n_background = 21.74; 
        results_SR3_WZ_0Ja.background_sys = 2.85;
        results_SR3_WZ_0Ja.signal_sys = 0.; 
        results_SR3_WZ_0Ja.n_signal = _numSR3_WZ_0Ja;
	add_result(results_SR3_WZ_0Ja);

        SignalRegionData results_SR3_WZ_0Jb;
        results_SR3_WZ_0Jb.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR3_WZ_0Jb.sr_label = "SR3_WZ_0Jb";
        results_SR3_WZ_0Jb.n_observed = 1.;
        results_SR3_WZ_0Jb.n_background = 2.68; 
        results_SR3_WZ_0Jb.background_sys = 0.46;
        results_SR3_WZ_0Jb.signal_sys = 0.; 
        results_SR3_WZ_0Jb.n_signal = _numSR3_WZ_0Jb;
	add_result(results_SR3_WZ_0Jb);

        SignalRegionData results_SR3_WZ_0Jc;
        results_SR3_WZ_0Jc.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR3_WZ_0Jc.sr_label = "SR3_WZ_0Jc";
        results_SR3_WZ_0Jc.n_observed = 2.;
        results_SR3_WZ_0Jc.n_background = 1.56; 
        results_SR3_WZ_0Jc.background_sys = 0.33;
        results_SR3_WZ_0Jc.signal_sys = 0.; 
        results_SR3_WZ_0Jc.n_signal = _numSR3_WZ_0Jc;
	add_result(results_SR3_WZ_0Jc);

        SignalRegionData results_SR3_WZ_1Ja;
        results_SR3_WZ_1Ja.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR3_WZ_1Ja.sr_label = "SR3_WZ_1Ja";
        results_SR3_WZ_1Ja.n_observed = 1.;
        results_SR3_WZ_1Ja.n_background = 2.21; 
        results_SR3_WZ_1Ja.background_sys = 0.53;
        results_SR3_WZ_1Ja.signal_sys = 0.; 
        results_SR3_WZ_1Ja.n_signal = _numSR3_WZ_1Ja;
	add_result(results_SR3_WZ_1Ja);

        SignalRegionData results_SR3_WZ_1Jb;
        results_SR3_WZ_1Jb.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR3_WZ_1Jb.sr_label = "SR3_WZ_1Jb";
        results_SR3_WZ_1Jb.n_observed = 3.;
        results_SR3_WZ_1Jb.n_background = 1.82; 
        results_SR3_WZ_1Jb.background_sys = 0.26;
        results_SR3_WZ_1Jb.signal_sys = 0.; 
        results_SR3_WZ_1Jb.n_signal = _numSR3_WZ_1Jb;
	add_result(results_SR3_WZ_1Jb);

        SignalRegionData results_SR3_WZ_1Jc;
        results_SR3_WZ_1Jb.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR3_WZ_1Jb.sr_label = "SR3_WZ_1Jc";
        results_SR3_WZ_1Jb.n_observed = 4.;
        results_SR3_WZ_1Jb.n_background = 1.26; 
        results_SR3_WZ_1Jb.background_sys = 0.34;
        results_SR3_WZ_1Jb.signal_sys = 0.; 
        results_SR3_WZ_1Jb.n_signal = _numSR3_WZ_1Jc;
	add_result(results_SR3_WZ_1Jc);

      }

      vector<vector<HEPUtils::Particle*>> getOSSFpair(vector<HEPUtils::Particle*> leptons) {
        vector<vector<HEPUtils::Particle*>> OSSFpair_container;
        for (size_t iLe1=0;iLe1<leptons.size();iLe1++) {
          for (size_t iLe2=0;iLe2<leptons.size();iLe2++) {
            if (leptons.at(iLe1)->abspid()==leptons.at(iLe2)->abspid() && leptons.at(iLe1)->pid()!=leptons.at(iLe2)->pid()) {
              vector<HEPUtils::Particle*> OSSFpair;
              OSSFpair.push_back(leptons.at(iLe1));
              OSSFpair.push_back(leptons.at(iLe2));
              OSSFpair_container.push_back(OSSFpair);
            }
          }
        }
        return OSSFpair_container;
      }

      vector<vector<HEPUtils::Particle*>> getOSpair(vector<HEPUtils::Particle*> leptons) {
        vector<vector<HEPUtils::Particle*>> OSpair_container;
        for (size_t iLe1=0;iLe1<leptons.size();iLe1++) {
          for (size_t iLe2=0;iLe2<leptons.size();iLe2++) {
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

      vector<HEPUtils::P4> get_W_ISR(vector<HEPUtils::Jet*> jets, HEPUtils::P4 Z, HEPUtils::P4 met) {
	HEPUtils::P4 Z_met_sys=Z+met;
	double deltaR_min=999;
	size_t Wjets_id1;
	size_t Wjets_id2;
	for (size_t i=0;i<jets.size();i++) {
	  for (size_t j=0;j<jets.size();j++) {
	    if (i!=j) {
	      HEPUtils::P4 jj_sys=jets.at(i)->mom()+jets.at(j)->mom();
	      double deltaR=fabs(jj_sys.deltaR_eta(Z_met_sys));
	      if (deltaR<deltaR_min) {
		deltaR_min=deltaR;
		Wjets_id1=i;
		Wjets_id2=j;
	      }
	    }
	  }
	}
	HEPUtils::P4 W=jets.at(Wjets_id2)->mom()+jets.at(Wjets_id1)->mom();
	HEPUtils::P4 ISR;
	for (size_t k=0;k<jets.size();k++) {
	  if ((k=!Wjets_id1) && (k!=Wjets_id2))ISR+=(jets.at(k)->mom());
	}
	vector<HEPUtils::P4> W_ISR;
	W_ISR.push_back(W);
	W_ISR.push_back(ISR);
	return W_ISR;
      }
		 

    };


    // Factory fn
    DEFINE_ANALYSIS_FACTORY(ATLAS_13TeV_MultiLEP_36invfb)


  }
}
