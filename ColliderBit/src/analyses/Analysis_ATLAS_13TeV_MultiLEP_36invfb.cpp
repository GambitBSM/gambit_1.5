///
///  \author Rose Kudzman-Blais
///  \date 2017 May
///
///  *********************************************


#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>

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
      double _numSR2_SF_loose, _numSR2_SF_tight, _numSR2_DF_100, _numSR2_DF_150, _numSR2_DF_200, _numSR2_DF_300, _numSR2_int, _numSR2_high, _numSR2_low_2J, _num_low_3J, _numSR3_slep_a, _numSR3_slep_b, _numSR3_slep_c, _numSR3_slep_d, _numSR3_slep_e, _numSR3_WZ_0Ja, _numSR3_WZ_0Jb, _numSR3_WZ_0Jc, _numSR3_WZ_1Ja, _numSR3_WZ_1Jb, _numSR3_WZ_1Jc; 
      vector<int> cutFlowVector1, cutFlowVector2, cutFlowVector3, cutFlowVector4;
      vector<string> cutFlowVector_str1, cutFlowVector_str2, cutFlowVector_str3, cutFlowVector_str4;
      size_t NCUTS1, NCUTS2, NCUTS3, NCUTS4;

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
	_numSR2_low_2J=0;
	_num_low_3J=0;
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


	NCUTS1=10;
	NCUTS2=7;
	NCUTS3=7;
	NCUTS4=5;
	set_luminosity(35.9);

        for (size_t i=0;i<NCUTS1;i++){
          cutFlowVector1.push_back(0);
          cutFlowVector_str1.push_back("");
        }
        for (size_t i=0;i<NCUTS2;i++){
          cutFlowVector2.push_back(0);
          cutFlowVector_str2.push_back("");
        }
        for (size_t i=0;i<NCUTS3;i++){
          cutFlowVector3.push_back(0);
          cutFlowVector_str3.push_back("");
        }
        for (size_t i=0;i<NCUTS4;i++){
          cutFlowVector4.push_back(0);
          cutFlowVector_str4.push_back("");
        }

        analysisRunName = "ATLAS_13TeV_MultiLEP_36invfb_test";
        vector<const char*> variables = {"met"};
        plots = new Perf_Plot(analysisRunName, &variables);

      }


      void analyze(const HEPUtils::Event* event) {
	HEPUtilsAnalysis::analyze(event);

        // Missing energy
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
	  if (muon->pT()>10. &&fabs(muon->eta())<2.5) {
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

        vector<HEPUtils::Jet*> baselinePhotons;
        for (HEPUtils::Jet* photons : event->photons()) {
          if (photon->pT()>25. && fabs(photon->eta())<2.37 && (fabs(photon->eta())<1.37 || fabs(photon->eta())>1.52)) {
            baselineJets.push_back(jet);
          }
        }

	//Overlap Removal + Signal Objects	
	vector<HEPUtils::Particle*> signalElectrons;
	vector<HEPUtils::Particle*> signalMuons;
	vector<HEPUtils::Particle*> signalLeptons;
	vector<HEPUtils::Jet*> signalJets;
	vector<HEPUtils::Jet*> signalBJets;
	vector<HEPUtils::Jet*> signalNonBJets;

	vector<size_t> overlapJet;
        for (size_t iJet=0;iJet<baselineJets.size();iJet++) {
	  vector<size_t> overlapEl;
          for (size_t iEl=0;iEl<baselineElectrons.size();iEl++) {
            if (fabs(baselineElectrons.at(iEl)->mom().deltaR_eta(baselineJets.at(iJet)->mom()))<0.2)overlapEl.push_back(iEl);
	  }
	  if (overlapEl.size()>0 && baselineJets.at(iJet)->btag()) {
	    overlapJet.push_back(iJet)
	    for (size_t iO=0;iO<overlapEl.size();iO++) {
	      baselineElectrons.erase(baselineElectrons.begin()+overlapEl.at(iO));  
	    }
	  }
	 if (overlap.size()==0)overlapJet.push_back(iJet);
	}
	for (size_t iO=0;overlapJet.size();i0++) {
	  baselineJets.erase(baselineJets.begin()+overlapJet.at(iO));
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
	    signalJets.push_back(baselineJets.at(iJet));
	    if (baselineJets.at(iJet)->btag())signalBJets.push_back(baselineJets.at(iJet));
	    if (!baselineJets.at(iJet)->btag())signalNonBJets.push_back(baselineJets.at(iJet));
	  }
	}

	for (size_t iMu=0;iMu<baselineMuons.size();iMu++) {
	  bool overlap=false;
	  for (size_t iJet=0;iJet<signalJets.size();iJet++) {
            if (fabs(baselineMuons.at(iMu)->mom().deltaR_eta(signalJets.at(iJet)->mom()))<0.4)overlap=true;
	    
	  }
	  if (!overlap) {
	    signalMuons.push_back(baselineMuons.at(iEl));
	    signalLeptons.push_back(baselineMuons.at(iEl));
	  }
	}

	sort(signalJets.begin(),signalJets.end(),compareJetPt);	
        sort(signalLeptons.begin(),signalLeptons.end(),comparePt);
        int nSignalLeptons = signalLeptons.size();
        int nSignalJets = signalJets.size();
	int nSignalNonBJets = signalNonBJets.size();       
 
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
	for (size_t iJet;iJet<nSignalJets;iJet++) {
	  if (signalJets.at(iJet)->pT()>60 && fabs(signalJets.at(iJet)->eta())<2.4 && signalJets.at(iJet)->btag()==false)large_jet_veto=false;
	}

	//Bjet veto
	bool bjet_veto=true;
	for (size_t iJet;iJet<nSignalBJets;iJet++);
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
		HEPUtils::P4 Z=signalLeptons.at(0)->mom()+signalLeptons.at(1)->mom()
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
	      if (nSignalNonBJets==2 && mll>81. && mll<101. && mjj>70. && mjj<90. && met>100. && Z_pT>60.) {
		HEPUtils::P4 W=signalJets.at(0)->mom()+signalJets.at(1)->mom();
	        double deltaPhi_met_W=fabs(W.phi()-event->missingmom().phi());
	        double deltaPhi_met_Z=fabs(Z.phi()-event->missingmom().phi());
		if (deltaPhi_met_Z<0.8 && deltaPhi_met_W>1.5 && (met/Z.pT())>0.6 && (met/Z.pT())<1.6 && (met/W.pT())<0.8)_numSR2_low_2J++;
	      }	
	      //SR2_low_3J
	      if (nSignalNonBJets>2 && nSignalNonBJets<6 && mll>86 && mll<96 && mjj>70 && mjj<90 && met>100 && Z.pT()>40 && deltaR_jj<2.2) {
		vector<HEPUtils::P4> W_ISR=get_W_ISR(signalJets,Z,event->missingmom());
		double deltaPhi_met_W=fabs(W_ISR.at(0).phi()-event->missingmom().phi());
		double deltaPhi_met_ISR=fabs(W_ISR.at(1).phi()-event->missingmom().phi());
		double deltaPhi_met_jet1=fabs(signalJets.at(0).phi()-event->missingmom().phi());
		if (deltaPhi_met_W<2.2 && delta_met_ISR>2.4 && deltaPhi_met_jet1>2.6 && (met/W_ISR.at(1).pT())>0.4 && (met/W_ISR.at(1).pT())<0.8 && fabs(Z.eta())<1.6 && signalJets.at(2)->pT()>30.)_numSR2_low_3J++;
	    }

	  }	
	}
	
	//3lep
	if (preselection && nSignalLeptons==3 && bjet_veto && OSSFpairs.size()) {
	  double mTmin=999;
	  double mSFOS=999;
	  for (size_t iPa=0;iPa<OSSFpairs.size();iPa++) {
	    for (size_t iLep;iLep<signalLeptons.size();iLep++) {
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
	    if (met>170.)_numSR3_WZ_0Jc;
	  }
	  if (mSFOS>81.2 && mSFOS<101.2 && nSignalNonBJets>0) {
	    if (met>120. && met<200. && mTmin>110. && (signalLeptons.at(0)->mom()+signalLeptons.at(1)->mom()+signalLeptons.at(2)->mom()).pT()<120. && signalJets.at(0)>70.)_numSR3_WZ_1Ja++;
	    if (met>200. && mTmin>110. && mTmin<160.)_numSR3_WZ_1Jb++;
	    if (met>200. && signalLeptons.at(2)->pT()>35. && mTmin>160.)_numSR3_WZ_1Jc+;
	  }
	}

		
        if (preselection) {
          vector<double> variables = {met};
          plots->fill(&variables);
        }

        cutFlowVector_str1[0] = "All events";
        cutFlowVector_str1[1] = "2 light leptons";
        cutFlowVector_str1[2] = "Same-sign";
        cutFlowVector_str1[3] = "3rd lepton veto";
        cutFlowVector_str1[4] = "Low mass veto";
        cutFlowVector_str1[5] = "Bjet veto";
        cutFlowVector_str1[6] = "met > 60 GeV";
        cutFlowVector_str1[7] = "0 or 1 ISR jet";
        cutFlowVector_str1[8] = "mT < 100 GeV";
        cutFlowVector_str1[9] = "pT_ll > 100 GeV";

        cutFlowVector_str2[0] = "All events";
        cutFlowVector_str2[1] = "3 leptons";
        cutFlowVector_str2[2] = "Low mass & conversions veto";
        cutFlowVector_str2[3] = "Bjet veto";
        cutFlowVector_str2[4] = "met > 50 GeV";
        cutFlowVector_str2[5] = "mT > 100 GeV";
        cutFlowVector_str2[6] = "mll > 75 GeV";

        cutFlowVector_str3[0] = "All events";
        cutFlowVector_str3[1] = "3 leptons";
        cutFlowVector_str3[2] = "Low mass & conversions veto";
        cutFlowVector_str3[3] = "Bjet veto";
        cutFlowVector_str3[4] = "met > 50 GeV";
        cutFlowVector_str3[5] = "mT2 < 100 GeV";
        cutFlowVector_str3[6] = "mll < 75 GeV";

        cutFlowVector_str4[0] = "All events";
        cutFlowVector_str4[1] = "4 leptons";
        cutFlowVector_str4[2] = "Low mass veto";
        cutFlowVector_str4[3] = "Bjet veto";
        cutFlowVector_str4[4] = "met > 100 GeV";

        for (size_t j=0;j<NCUTS1;j++){
          if(
             (j==0) ||

	     (j==1 && nSignalLightLeptons==2) ||

	     (j==2 && nSignalLightLeptons==2) ||

	     (j==3 && nSignalLightLeptons==2 && signalLightLeptons.at(0)->pid()*signalLightLeptons.at(1)->pid()>0) ||

	     (j==4 && nSignalLightLeptons==2 && signalLightLeptons.at(0)->pid()*signalLightLeptons.at(1)->pid()>0 && nSignalLeptons==2) ||

	     (j==5 && nSignalLightLeptons==2 && signalLightLeptons.at(0)->pid()*signalLightLeptons.at(1)->pid()>0 && nSignalLeptons==2 && low_mass_veto) ||

	     (j==6 && nSignalLightLeptons==2 && signalLightLeptons.at(0)->pid()*signalLightLeptons.at(1)->pid()>0 && nSignalLeptons==2 && low_mass_veto && bjet_veto)  ||

	     (j==7 && nSignalLightLeptons==2 && signalLightLeptons.at(0)->pid()*signalLightLeptons.at(1)->pid()>0 && nSignalLeptons==2 && low_mass_veto && bjet_veto && met>60 && ISRjet) ||            
 
             (j==8 && nSignalLightLeptons==2 && signalLightLeptons.at(0)->pid()*signalLightLeptons.at(1)->pid()>0 && nSignalLeptons==2 && low_mass_veto && bjet_veto && met>60 && ISRjet && mT<100) ||
             
             (j==9 && nSignalLightLeptons==2 && signalLightLeptons.at(0)->pid()*signalLightLeptons.at(1)->pid()>0 && nSignalLeptons==2 && low_mass_veto && bjet_veto && met>60 && ISRjet && mT<100 && pT_ll>100) )

	  cutFlowVector1[j]++;
	}

        for (size_t j=0;j<NCUTS2;j++){
          if(
             (j==0) ||

             (j==1 && nSignalLeptons==3) ||

             (j==2 && nSignalLeptons==3 && low_mass_veto && conversion_veto) ||

             (j==3 && nSignalLeptons==3 && low_mass_veto && conversion_veto && bjet_veto) ||

             (j==4 && nSignalLeptons==3 && low_mass_veto && conversion_veto && bjet_veto && met>50.) ||

             (j==5 && nSignalLeptons==3 && low_mass_veto && conversion_veto && bjet_veto && met>50. && mT>100.) ||

             (j==6 && nSignalLeptons==3 && low_mass_veto && conversion_veto && bjet_veto && met>50. && mT>100. && mll>75.) )

          cutFlowVector2[j]++;
        }

        for (size_t j=0;j<NCUTS3;j++){
          if(
             (j==0) ||

             (j==1 && nSignalLeptons==3) ||

             (j==2 && nSignalLeptons==3 && low_mass_veto && conversion_veto) ||

             (j==3 && nSignalLeptons==3 && low_mass_veto && conversion_veto && bjet_veto) ||

             (j==4 && nSignalLeptons==3 && low_mass_veto && conversion_veto && bjet_veto && met>50.) ||

             (j==5 && nSignalLeptons==3 && low_mass_veto && conversion_veto && bjet_veto && met>50. && mT<100.) ||

             (j==6 && nSignalLeptons==3 && low_mass_veto && conversion_veto && bjet_veto && met>50. && mT<100. && mll<75.) )

          cutFlowVector3[j]++;
        }

        for (size_t j=0;j<NCUTS4;j++){
          if(
             (j==0) ||
        
             (j==1 && nSignalLeptons==4) ||
        
             (j==2 && nSignalLeptons==4 && low_mass_veto) ||
        
             (j==3 && nSignalLeptons==4 && low_mass_veto && bjet_veto) ||
      
             (j==4 && nSignalLeptons==4 && low_mass_veto && bjet_veto && met>100.) )

          cutFlowVector4[j]++;
        }
      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        
	HEPUtilsAnalysis::add(other);

        Analysis_ATLAS_13TeV_MultiLEP_36invfb* specificOther
                = dynamic_cast<Analysis_ATLAS_13TeV_MultiLEP_36invfb*>(other);

        // Here we will add the subclass member variables:
        if (NCUTS1 != specificOther->NCUTS1) NCUTS1 = specificOther->NCUTS1;
	if (NCUTS2 != specificOther->NCUTS2) NCUTS2 = specificOther->NCUTS2;
	if (NCUTS3 != specificOther->NCUTS3) NCUTS3 = specificOther->NCUTS3;
	if (NCUTS4 != specificOther->NCUTS4) NCUTS4 = specificOther->NCUTS4;
        for (size_t j = 0; j < NCUTS1; j++) {
          cutFlowVector1[j] += specificOther->cutFlowVector1[j];
          cutFlowVector_str1[j] = specificOther->cutFlowVector_str1[j];
        }
        for (size_t j = 0; j < NCUTS2; j++) {
          cutFlowVector2[j] += specificOther->cutFlowVector2[j];
          cutFlowVector_str2[j] = specificOther->cutFlowVector_str2[j];
        }
        for (size_t j = 0; j < NCUTS3; j++) {
          cutFlowVector3[j] += specificOther->cutFlowVector3[j];
          cutFlowVector_str3[j] = specificOther->cutFlowVector_str3[j];
        }
        for (size_t j = 0; j < NCUTS4; j++) {
          cutFlowVector4[j] += specificOther->cutFlowVector4[j];
          cutFlowVector_str4[j] = specificOther->cutFlowVector_str4[j];
        }
        _numSR1 += specificOther->_numSR1;
        _numSR2 += specificOther->_numSR2;
        _numSR3 += specificOther->_numSR3;
        _numSR4 += specificOther->_numSR4;
        _numSR5 += specificOther->_numSR5;
        _numSR6 += specificOther->_numSR6;
        _numSR7 += specificOther->_numSR7;
        _numSR8 += specificOther->_numSR8;
      }


      void collect_results() {

        string path = "ColliderBit/results/cutflow_";
        path.append(analysisRunName);
        path.append(".txt");
        cutflowFile.open(path.c_str());
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
        cutflowFile << "CUT FLOW: CMS Multi-lepton paper "<<endl;
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------"<<endl;
        cutflowFile<< right << setw(60) << "CUT" << setw(20) << "RAW" << setw(20) << " % " << endl;
        for (size_t j=0; j<NCUTS1; j++) {
          cutflowFile << right << setw(60) << cutFlowVector_str1[j].c_str() << setw(20) << cutFlowVector1[j] << setw(20) << 100.*cutFlowVector1[j]/cutFlowVector1[0] << endl;
        }
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
        cutflowFile<< right << setw(60) << "CUT" << setw(20) << "RAW" << setw(20) << " % " << endl;
        for (size_t j=0; j<NCUTS2; j++) {
          cutflowFile << right << setw(60) << cutFlowVector_str2[j].c_str() << setw(20) << cutFlowVector2[j] << setw(20) << 100.*cutFlowVector2[j]/cutFlowVector2[0] << endl;
        }
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
        cutflowFile<< right << setw(60) << "CUT" << setw(20) << "RAW" << setw(20) << " % " << endl;
        for (size_t j=0; j<NCUTS3; j++) {
          cutflowFile << right << setw(60) << cutFlowVector_str3[j].c_str() << setw(20) << cutFlowVector3[j] << setw(20) << 100.*cutFlowVector3[j]/cutFlowVector1[0] << endl;
        }
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
        cutflowFile<< right << setw(60) << "CUT" << setw(20) << "RAW" << setw(20) << " % " << endl;
        for (size_t j=0; j<NCUTS4; j++) {
          cutflowFile << right << setw(60) << cutFlowVector_str4[j].c_str() << setw(20) << cutFlowVector4[j] << setw(20) << 100.*cutFlowVector4[j]/cutFlowVector4[0] << endl;
        }
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
        cutflowFile.close();

        plots->createFile();

        //Now fill a results object with the results for each SR
        SignalRegionData results_SR1;
        results_SR1.analysis_name = "Analysis_ATLAS_13TeV_MultiLEP_36invfb";
        results_SR1.sr_label = "SR1";
        results_SR1.n_observed = 13.;
        results_SR1.n_background = 12.; 
        results_SR1.background_sys = 3.;
        results_SR1.signal_sys = 0.; 
        results_SR1.n_signal = _numSR1;
	add_result(results_SR1);

      }

      struct ptComparison {
        bool operator() (HEPUtils::Particle* i,HEPUtils::Particle* j) {return (i->pT()>j->pT());}
      } comparePt;
      
      struct ptJetComparison {
        bool operator() (HEPUtils::Jet* i,HEPUtils::Jet* j) {return (i->pT()>j->pT());}
      } compareJetPt;

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
	for (size_t i;i<jets.size();i++) {
	  for (size_t j;j<jets.size();j++) {
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
	  if (k=!Wjets_id1 && k!=Wjets_id2)ISR+=jets.at(k)->mom();
	}
	vector<HEPUtils::P4> W_ISR=W+ISR;
	return W_ISR;
      }
		 


    };


    // Factory fn
    DEFINE_ANALYSIS_FACTORY(ATLAS_13TeV_MultiLEP_36invfb)


  }
}
