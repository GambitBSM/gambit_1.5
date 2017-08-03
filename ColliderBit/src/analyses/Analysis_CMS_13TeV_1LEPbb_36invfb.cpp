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
#include "gambit/ColliderBit/analyses/Perf_Plot.hpp"

using namespace std;

namespace Gambit {
  namespace ColliderBit {

    class Analysis_CMS_13TeV_1LEPbb_36invfb : public HEPUtilsAnalysis {
    private:

      double _numSRA, _numSRB; 

      vector<int> cutFlowVector;
      vector<string> cutFlowVector_str;
      vector<double> cutFlowVectorCMS_225_75;
      double xsecCMS_225_75;
      size_t NCUTS;

      Perf_Plot* plots;
      ofstream cutflowFile;
      string analysisRunName;

    public:

      Analysis_CMS_13TeV_1LEPbb_36invfb() {

        _numSRA=0;
        _numSRB=0;

        NCUTS=10;
        set_luminosity(35.9);

	xsecCMS_225_75=1165;
        for (size_t i=0;i<NCUTS;i++){
          cutFlowVector.push_back(0);
	  cutFlowVectorCMS_225_75.push_back(0);
          cutFlowVector_str.push_back("");
        }

        analysisRunName = "CMS_13TeV_1LEPbb_36invfb_225_75";
      }


      void analyze(const HEPUtils::Event* event) {
	HEPUtilsAnalysis::analyze(event);
        double met = event->met();

        // Baseline objects
        const vector<double> a={0,10.};
        const vector<double> b={0,10000.};
        const vector<double> cEl={0.83};
        HEPUtils::BinnedFn2D<double> _eff2dEl(a,b,cEl);
        vector<HEPUtils::Particle*> baselineElectrons;
        for (HEPUtils::Particle* electron : event->electrons()) {
          bool hasTrig=has_tag(_eff2dEl, electron->eta(), electron->pT());
          if (electron->pT()>5. && electron->abseta()<2.5 && hasTrig)baselineElectrons.push_back(electron);
        }

        const vector<double> cMu={0.89};
        HEPUtils::BinnedFn2D<double> _eff2dMu(a,b,cMu);
        vector<HEPUtils::Particle*> baselineMuons;
        for (HEPUtils::Particle* muon : event->muons()) {
          bool hasTrig=has_tag(_eff2dMu, muon->eta(), muon->pT());
          if (muon->pT()>5. && muon->abseta()<2.4 && hasTrig)baselineMuons.push_back(muon);
        }

        vector<HEPUtils::Particle*> baselineTaus;
        for (HEPUtils::Particle* tau : event->taus()) {
          if (tau->pT()>20. && tau->abseta()<2.3)baselineTaus.push_back(tau);
        }

        vector<HEPUtils::Jet*> baselineJets;
        for (HEPUtils::Jet* jet : event->jets()) {
          if (jet->pT()>25. &&fabs(jet->eta())<2.4)baselineJets.push_back(jet);
        }

        // Signal objects
        vector<HEPUtils::Particle*> signalLeptons;
        vector<HEPUtils::Particle*> signalElectrons;
        vector<HEPUtils::Particle*> signalMuons;
        vector<HEPUtils::Jet*> signalJets;   
        vector<HEPUtils::Jet*> signalBJets;

	for (size_t iEl=0;iEl<baselineElectrons.size();iEl++) {
	  if (baselineElectrons.at(iEl)->pT()>30. && baselineElectrons.at(iEl)->abseta()<1.44)signalElectrons.push_back(baselineElectrons.at(iEl));
	} 
       
        for (size_t iMu=0;iMu<baselineMuons.size();iMu++) {
          if (baselineMuons.at(iMu)->pT()>25. && baselineMuons.at(iMu)->abseta()<2.1)signalMuons.push_back(baselineMuons.at(iMu));
        }

        const vector<double> c={0.65};
        HEPUtils::BinnedFn2D<double> _eff2d(a,b,c);
	for (size_t iJet=0;iJet<baselineJets.size();iJet++) {
           if (baselineJets.at(iJet)->pT()>30.) {
	    signalJets.push_back(baselineJets.at(iJet));                
            bool hasTag=has_tag(_eff2d, baselineJets.at(iJet)->eta(), baselineJets.at(iJet)->pT()); 
	    if (baselineJets.at(iJet)->btag() && hasTag)signalBJets.push_back(baselineJets.at(iJet));
	  }
        }

	signalLeptons=signalElectrons;
	signalLeptons.insert(signalLeptons.end(),signalMuons.begin(),signalMuons.end());
        int nSignalLeptons=signalLeptons.size();
        int nSignalJets=signalJets.size();
        int nSignalBJets=signalBJets.size();

	//Variables
	bool preselection=false;	
        bool lepton2_veto=true;
	bool tau_veto=true;
        double mCT=0;
        double mbb=0;
	double mT=0;

	if ((baselineMuons.size()+baselineElectrons.size())>1)lepton2_veto=false;
	if (baselineTaus.size()>0)tau_veto=false;
        if (nSignalLeptons>0 && met>50. && lepton2_veto && tau_veto && nSignalJets==2 && nSignalBJets==2)preselection=true;
	
	if (nSignalBJets>1) {
	  mCT=sqrt(2*signalBJets.at(0)->pT()*signalBJets.at(1)->pT()*(1+cos(signalBJets.at(0)->mom().deltaPhi(signalBJets.at(1)->mom()))));
          mbb=(signalBJets.at(0)->mom()+signalBJets.at(1)->mom()).m();
	}
	if (signalLeptons.size()>0)mT=sqrt(2*signalLeptons.at(0)->pT()*met*(1-cos(signalLeptons.at(0)->mom().deltaPhi(event->missingmom()))));

	//Signal Regions
        if (preselection && mbb>90 && mbb<150 && mCT>170. && met>125. && mT>150.) {
          //SRA
          if (met>125. && met<200.)_numSRA++;
          //SRB
          if (met>200.)_numSRB++;   
	}

        cutFlowVector_str[0] = "All events";
        cutFlowVector_str[1] = "$\\geq$ 1 signal lepton; $E_{T}^{miss} > 50 GeV$";
        cutFlowVector_str[2] = "2nd lepton veto";
        cutFlowVector_str[3] = "Tau veto";
        cutFlowVector_str[4] = "2 jets";
        cutFlowVector_str[5] = "2 bjets";
        cutFlowVector_str[6] = "$90 < m_{bb} < 150 GeV$";
        cutFlowVector_str[7] = "$m_{CT} > 170 GeV$";
        cutFlowVector_str[8] = "$E_{T}^{miss} > 125 GeV$";
        cutFlowVector_str[9] = "$m_{T} > 150 GeV$";

	cutFlowVectorCMS_225_75[0]=7297.6;	
	cutFlowVectorCMS_225_75[1]=1320.5;	
	cutFlowVectorCMS_225_75[2]=1265.3;	
	cutFlowVectorCMS_225_75[3]=1259.0;	
	cutFlowVectorCMS_225_75[4]=680.8;	
	cutFlowVectorCMS_225_75[5]=299.0;	
	cutFlowVectorCMS_225_75[6]=258.4;	
	cutFlowVectorCMS_225_75[7]=50.9;	
	cutFlowVectorCMS_225_75[8]=38.4;	
	cutFlowVectorCMS_225_75[9]=4.7;	

        for (size_t j=0;j<NCUTS;j++){
          if(
             (j==0) ||

	     (j==1 && nSignalLeptons>=1 && met>50) ||

	     (j==2 && nSignalLeptons>=1 && met>50 && lepton2_veto) ||

	     (j==3 && nSignalLeptons>=1 && met>50 && lepton2_veto && tau_veto) ||

	     (j==4 && nSignalLeptons>=1 && met>50 && lepton2_veto && tau_veto && nSignalJets==2) ||

	     (j==5 && nSignalLeptons>=1 && met>50 && lepton2_veto && tau_veto && nSignalJets==2 && nSignalBJets==2) ||

	     (j==6 && nSignalLeptons>=1 && met>50 && lepton2_veto && tau_veto && nSignalJets==2 && nSignalBJets==2 && mbb>90 && mbb<150) ||

	     (j==7 && nSignalLeptons>=1 && met>50 && lepton2_veto && tau_veto && nSignalJets==2 && nSignalBJets==2 && mbb>90 && mbb<150 && mCT>170.) ||            
 
             (j==8 && nSignalLeptons>=1 && met>50 && lepton2_veto && tau_veto && nSignalJets==2 && nSignalBJets==2 && mbb>90 && mbb<150 && mCT>170. && met>125.) ||
             
             (j==9 && nSignalLeptons>=1 && met>50 && lepton2_veto && tau_veto && nSignalJets==2 && nSignalBJets==2 && mbb>90 && mbb<150 && mCT>170. && met>125. && mT>150.) )

            cutFlowVector[j]++;
	}

      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        HEPUtilsAnalysis::add(other);

        Analysis_CMS_13TeV_1LEPbb_36invfb* specificOther
                = dynamic_cast<Analysis_CMS_13TeV_1LEPbb_36invfb*>(other);

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

        if (analysisRunName.find("225_75") != string::npos) {
          cutflowFile<<"\\begin{table}[H] \n\\caption{$\\tilde{\\chi}_{1}^{\\pm}\\tilde{\\chi}_{2}^{0}$ decay via $W/h, [\\tilde{\\chi}_{2}^{0}\\tilde{\\chi}_{1}^{\\pm},\\tilde{\\chi}_{1}^{0}]: [225,75] [GeV]$} \n\\makebox[\\linewidth]{ \n\\renewcommand{\\arraystretch}{0.4} \n\\begin{tabular}{c c c c c} \n\\hline"<<endl;
          cutflowFile<<"& ATLAS & GAMBIT & GAMBIT/ATLAS & $\\sigma$-corrected GAMBIT/ATLAS \\\\ \\hline"<<endl;
          cutflowFile<<"$\\sigma (pp\\to \\tilde{\\chi}_{1}^{\\pm}, \\tilde{\\chi}_{2}^{0})$ &"<<setprecision(4)<<xsecCMS_225_75<<" $fb$ &"<<setprecision(4)<<xsec()<<"$fb$ &"<<setprecision(4)<< xsec()/xsecCMS_225_75<<" & 1\\\\ \\hline"<<endl;
          cutflowFile<<"\\multicolumn{5}{c}{Expected events at 35.9 $fb^{-1}$} \\\\ \\hline"<<endl;
          for (size_t i=0; i<NCUTS; i++) {
            cutflowFile<<cutFlowVector_str[i]<<"&"<<setprecision(4)<<cutFlowVectorCMS_225_75[i]<<"&"<<setprecision(4)<<cutFlowVector[i]*xsec_per_event()*luminosity()<<"&"<<setprecision(4)<<cutFlowVector[i]*xsec_per_event()*luminosity()/cutFlowVectorCMS_225_75[i]<<"&"<<setprecision(4)<<(xsecCMS_225_75/xsec())*cutFlowVector[i]*xsec_per_event()*luminosity()/cutFlowVectorCMS_225_75[i]<<"\\\\"<< endl;
          }
          cutflowFile<<"\\hline \\multicolumn{5}{c}{Percentage (\\%)} \\\\ \\hline"<<endl;
          for (size_t i=0; i<NCUTS; i++) {
            cutflowFile<<cutFlowVector_str[i]<<"&"<<setprecision(4)<<cutFlowVectorCMS_225_75[i]*100./cutFlowVectorCMS_225_75[1]<<"&"<<setprecision(4)<<cutFlowVector[i]*100./cutFlowVector[1]<<"& - & -\\\\"<< endl;
          }
          cutflowFile<<"\\end{tabular} \n} \n\\end{table}"<<endl;
	}
        cutflowFile.close();

        SignalRegionData results_SRA;
        results_SRA.analysis_name = "Analysis_CMS_13TeV_1LEPbb_36invfb";
        results_SRA.sr_label = "SRA";
        results_SRA.n_observed = 11.;
        results_SRA.n_background = 7.5; 
        results_SRA.background_sys = 2.5;
        results_SRA.signal_sys = 0.; 
        results_SRA.n_signal = _numSRA;
        add_result(results_SRA);

        SignalRegionData results_SRB;
        results_SRB.analysis_name = "Analysis_CMS_13TeV_1LEPbb_36invfb";
        results_SRB.sr_label = "SRB";
        results_SRB.n_observed = 7.;
        results_SRB.n_background = 8.7; 
        results_SRB.background_sys = 2.2;
        results_SRB.signal_sys = 0.; 
        results_SRB.n_signal = _numSRB;
        add_result(results_SRB);

      }

    };


    // Factory fn
    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_1LEPbb_36invfb)


  }
}
