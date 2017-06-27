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

    class Analysis_CMS_13TeV_2LEPsoft_36invfb : public HEPUtilsAnalysis {
    private:

      // Numbers passing cuts
      double _numSR1, _numSR2, _numSR3, _numSR4, _numSR5, _numSR6, _numSR7, _numSR8, _numSR9, _numSR10, _numSR11, _numSR12; 
      vector<int> cutFlowVector1;
      vector<string> cutFlowVector_str1;
      size_t NCUTS1;

      Perf_Plot* plots;
      ofstream cutflowFile;
      string analysisRunName;

    public:

      struct ptComparison {
        bool operator() (HEPUtils::Particle* i,HEPUtils::Particle* j) {return (i->pT()>j->pT());}
      } comparePt;

      Analysis_CMS_13TeV_2LEPsoft_36invfb() {

        _numSR1=0;
	_numSR2=0;
	_numSR3=0;
	_numSR4=0;
	_numSR5=0; 
	_numSR6=0;
	_numSR7=0; 
	_numSR8=0;
	_numSR9=0;
	_numSR10=0; 
	_numSR11=0;
	_numSR12=0; 

	NCUTS1=1;
	set_luminosity(35.9);

        for (size_t i=0;i<NCUTS1;i++){
          cutFlowVector1.push_back(0);
          cutFlowVector_str1.push_back("");
        }

        analysisRunName = "CMS_13TeV_2LEPsoft_36invfb_test";
        vector<const char*> variables = {"met"};
        plots = new Perf_Plot(analysisRunName, &variables);

      }


      void analyze(const HEPUtils::Event* event) {
	HEPUtilsAnalysis::analyze(event);

        // Missing energy
        double met = event->met();

        // Signal objects
        vector<HEPUtils::Particle*> signalLeptons;
        vector<HEPUtils::Particle*> signalElectrons;
        for (HEPUtils::Particle* electron : event->electrons()) {
	  if (electron->pT()>5. && electron->pT()<30.  &&fabs(electron->eta())<2.5) {
            signalElectrons.push_back(electron);
	    signalLeptons.push_back(electron);
          }
	}

        vector<HEPUtils::Particle*> signalMuons;
        for (HEPUtils::Particle* muon : event->muons()) {
	  if (muon->pT()>5. && muon->pT()<30. && fabs(muon->eta())<2.4) {
            signalMuons.push_back(muon);
	    signalLeptons.push_back(muon);
          }
	}

        size_t nSignalBJets;
	vector<HEPUtils::Jet*> signalJets;
	const vector<double>  a = {0,10.};
        const vector<double>  b = {0,10000.};
        const vector<double> c = {0.7};
        HEPUtils::BinnedFn2D<double> _eff2d(a,b,c);
        for (HEPUtils::Jet* jet : event->jets()) {
          bool hasTag=has_tag(_eff2d, jet->eta(), jet->pT());
          if (jet->pT()>25. && fabs(jet->eta())<2.4) {
	   signalJets.push_back(jet); 
	   if (jet->btag() && hasTag)nSignalBJets++;
	  }
        }

        sort(signalLeptons.begin(),signalLeptons.end(),comparePt);
        size_t nSignalLeptons = signalLeptons.size();
	size_t nSignalMuons = signalMuons.size();
	size_t nSignalJets = signalJets.size();

	//Variable definitions
	double m_ll=0;
	double pT_ll=0;
	double hT=0;
	double mTauTau=0;
	vector<double> mT;
	
	if (nSignalLeptons>=2) {
	  m_ll=(signalLeptons.at(0)->mom()+signalLeptons.at(1)->mom()).m();
	  pT_ll=(signalLeptons.at(0)->mom()+signalLeptons.at(1)->mom()).pT();

	  double determinant = signalLeptons.at(0)->mom().px()*signalLeptons.at(1)->mom().py()-signalLeptons.at(0)->mom().py()*signalLeptons.at(1)->mom().px();
    	  float xi_1 = (event->missingmom().px()*signalLeptons.at(1)->mom().py()-signalLeptons.at(1)->mom().px()*event->missingmom().py())/determinant;
    	  float xi_2 = (event->missingmom().py()*signalLeptons.at(0)->mom().px()-signalLeptons.at(0)->mom().py()*event->missingmom().px())/determinant;
  	  mTauTau = (1.+xi_1)*(1.+xi_2)*2*signalLeptons.at(0)->mom().dot(signalLeptons.at(1)->mom());
	}

	for (size_t iJet=0;iJet<nSignalJets;iJet++)hT+=signalJets.at(iJet)->pT();

	for (size_t iLep=0;iLep<nSignalLeptons;iLep++)mT.push_back(sqrt(2*signalLeptons.at(iLep)->pT()*met*(1-cos(signalLeptons.at(iLep)->phi()-event->missingmom().phi()))));
        

	//Preselection
        bool preselection=false; 

	if (nSignalLeptons==2 && nSignalBJets==0) {
	  if (signalLeptons.at(0)->pid()*signalLeptons.at(1)->pid()<0. && signalLeptons.at(0)->abspid()==signalLeptons.at(1)->abspid()) {
  	    if (m_ll<50. && pT_ll>3. && met>125. && met/hT<1.4 && met/hT>0.6 && hT>100. && m_ll>4. && (m_ll<9. || m_ll>0.5) && (mTauTau<0. || mTauTau>160.) && mT.at(0)<70. && mT.at(1)<70.) {
	      preselection=true;
	    }
	  }
	}	

	//SR
	if (preselection && met>125. && met<200. && nSignalMuons>0) {
	  if (m_ll>4. && m_ll<10.)_numSR1++;
	  if (m_ll>10. && m_ll<20.)_numSR2++;
	  if (m_ll>20. && m_ll<30.)_numSR3++;
	  if (m_ll>30. && m_ll<50.)_numSR4++;
	}
	if (preselection && met>200. && met<250.) {
	  if (m_ll>4. && m_ll<10.)_numSR5++;
	  if (m_ll>10. && m_ll<20.)_numSR6++;
	  if (m_ll>20. && m_ll<30.)_numSR7++;
	  if (m_ll>30. && m_ll<50.)_numSR8++;
	}
	if (preselection && met>250.) {
	  if (m_ll>4. && m_ll<10.)_numSR9++;
	  if (m_ll>10. && m_ll<20.)_numSR10++;
	  if (m_ll>20. && m_ll<30.)_numSR11++;
	  if (m_ll>30. && m_ll<50.)_numSR12++;
	}

        if (preselection) {
          vector<double> variables = {met};
          plots->fill(&variables);
        }

        cutFlowVector_str1[0] = "All events";

        for (size_t j=0;j<NCUTS1;j++){
          if(
             (j==0) )

          cutFlowVector1[j]++;
        }
      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        
	HEPUtilsAnalysis::add(other);

        Analysis_CMS_13TeV_2LEPsoft_36invfb* specificOther
                = dynamic_cast<Analysis_CMS_13TeV_2LEPsoft_36invfb*>(other);

        // Here we will add the subclass member variables:
        if (NCUTS1 != specificOther->NCUTS1) NCUTS1 = specificOther->NCUTS1;
        for (size_t j = 0; j < NCUTS1; j++) {
          cutFlowVector1[j] += specificOther->cutFlowVector1[j];
          cutFlowVector_str1[j] = specificOther->cutFlowVector_str1[j];
        }
        _numSR1 += specificOther->_numSR1;
        _numSR2 += specificOther->_numSR2;
        _numSR3 += specificOther->_numSR3;
        _numSR4 += specificOther->_numSR4;
        _numSR5 += specificOther->_numSR5;
        _numSR6 += specificOther->_numSR6;
        _numSR7 += specificOther->_numSR7;
        _numSR8 += specificOther->_numSR8;
        _numSR9 += specificOther->_numSR5;
        _numSR10 += specificOther->_numSR6;
        _numSR11 += specificOther->_numSR7;
        _numSR12 += specificOther->_numSR8;
      }


      void collect_results() {

        string path = "ColliderBit/results/cutflow_";
        path.append(analysisRunName);
        path.append(".txt");
        cutflowFile.open(path.c_str());
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
        cutflowFile << "CUT FLOW: CMS 2 soft leptons paper "<<endl;
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------"<<endl;
        cutflowFile<< right << setw(60) << "CUT" << setw(20) << "RAW" << setw(20) << " % " << endl;
        for (size_t j=0; j<NCUTS1; j++) {
          cutflowFile << right << setw(60) << cutFlowVector_str1[j].c_str() << setw(20) << cutFlowVector1[j] << setw(20) << 100.*cutFlowVector1[j]/cutFlowVector1[0] << endl;
        }
        cutflowFile << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
        cutflowFile.close();

        plots->createFile();

        //Now fill a results object with the results for each SR
        SignalRegionData results_SR1;
        results_SR1.analysis_name = "Analysis_CMS_13TeV_2LEPsoft_36invfb";
        results_SR1.sr_label = "SR1";
        results_SR1.n_observed = 2.;
        results_SR1.n_background = 3.5; 
        results_SR1.background_sys = 1.;
        results_SR1.signal_sys = 0.; 
        results_SR1.n_signal = _numSR1;
	add_result(results_SR1);

	SignalRegionData results_SR2;
        results_SR2.analysis_name = "Analysis_CMS_13TeV_2LEPsoft_36invfb";
        results_SR2.sr_label = "SR1";
        results_SR2.n_observed = 15.;
        results_SR2.n_background = 12.;
        results_SR2.background_sys = 2.3;
        results_SR2.signal_sys = 0.;
        results_SR2.n_signal = _numSR2;
        add_result(results_SR2);

        SignalRegionData results_SR3;
        results_SR3.analysis_name = "Analysis_CMS_13TeV_2LEPsoft_36invfb";
        results_SR3.sr_label = "SR3";
        results_SR3.n_observed = 19.;
        results_SR3.n_background = 17.;
        results_SR3.background_sys = 2.4;
        results_SR3.signal_sys = 0.;
        results_SR3.n_signal = _numSR3;
        add_result(results_SR3);

        SignalRegionData results_SR4;
        results_SR4.analysis_name = "Analysis_CMS_13TeV_2LEPsoft_36invfb";
        results_SR4.sr_label = "SR4";
        results_SR4.n_observed = 18.;
        results_SR4.n_background = 11.;
        results_SR4.background_sys = 2.;
        results_SR4.signal_sys = 0.;
        results_SR4.n_signal = _numSR4;
        add_result(results_SR4);

        SignalRegionData results_SR5;
        results_SR5.analysis_name = "Analysis_CMS_13TeV_2LEPsoft_36invfb";
        results_SR5.sr_label = "SR5";
        results_SR5.n_observed = 1.;
        results_SR5.n_background = 1.6;
        results_SR5.background_sys = 0.7;
        results_SR5.signal_sys = 0.;
        results_SR5.n_signal = _numSR5;
        add_result(results_SR5);

        SignalRegionData results_SR6;
        results_SR6.analysis_name = "Analysis_CMS_13TeV_2LEPsoft_36invfb";
        results_SR6.sr_label = "SR6";
        results_SR6.n_observed = 0.;
        results_SR6.n_background = 3.5;
        results_SR6.background_sys = 0.9;
        results_SR6.signal_sys = 0.;
        results_SR6.n_signal = _numSR6;
        add_result(results_SR6);

        SignalRegionData results_SR7;
        results_SR7.analysis_name = "Analysis_CMS_13TeV_2LEPsoft_36invfb";
        results_SR7.sr_label = "SR7";
        results_SR7.n_observed = 3.;
        results_SR7.n_background = 2.;
        results_SR7.background_sys = 0.7;
        results_SR7.signal_sys = 0.;
        results_SR7.n_signal = _numSR7;
        add_result(results_SR7);

        SignalRegionData results_SR8;
        results_SR8.analysis_name = "Analysis_CMS_13TeV_2LEPsoft_36invfb";
        results_SR8.sr_label = "SR8";
        results_SR8.n_observed = 1.;
        results_SR8.n_background = 0.51;
        results_SR8.background_sys = 0.51;
        results_SR8.signal_sys = 0.;
        results_SR8.n_signal = _numSR8;
        add_result(results_SR8);

        SignalRegionData results_SR9;
        results_SR9.analysis_name = "Analysis_CMS_13TeV_2LEPsoft_36invfb";
        results_SR9.sr_label = "SR9";
        results_SR9.n_observed = 2.;
        results_SR9.n_background = 1.4;
        results_SR9.background_sys = 0.7;
        results_SR9.signal_sys = 0.;
        results_SR9.n_signal = _numSR9;
        add_result(results_SR9);

        SignalRegionData results_SR10;
        results_SR10.analysis_name = "Analysis_CMS_13TeV_2LEPsoft_36invfb";
        results_SR10.sr_label = "SR10";
        results_SR10.n_observed = 1.;
        results_SR10.n_background = 1.5;
        results_SR10.background_sys = 0.6;
        results_SR10.signal_sys = 0.;
        results_SR10.n_signal = _numSR10;
        add_result(results_SR10);

        SignalRegionData results_SR11;
        results_SR11.analysis_name = "Analysis_CMS_13TeV_2LEPsoft_36invfb";
        results_SR11.sr_label = "SR11";
        results_SR11.n_observed = 2.;
        results_SR11.n_background = 1.5;
        results_SR11.background_sys = 0.8;
        results_SR11.signal_sys = 0.;
        results_SR11.n_signal = _numSR11;
        add_result(results_SR11);

        SignalRegionData results_SR12;
        results_SR12.analysis_name = "Analysis_CMS_13TeV_2LEPsoft_36invfb";
        results_SR12.sr_label = "SR12";
        results_SR12.n_observed = 0.;
        results_SR12.n_background = 1.2;
        results_SR12.background_sys = 0.6;
        results_SR12.signal_sys = 0.;
        results_SR12.n_signal = _numSR12;
        add_result(results_SR12);
      }

    };


    // Factory fn
    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_2LEPsoft_36invfb)


  }
}
