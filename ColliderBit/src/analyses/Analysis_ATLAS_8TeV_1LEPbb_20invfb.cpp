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

    class Analysis_ATLAS_8TeV_1LEPbb_20invfb : public HEPUtilsAnalysis {
    private:

      // Numbers passing cuts
      double _numSRA, _numSRB; 
      vector<int> cutFlowVector;
      vector<double> cutFlowVectorATLAS_130_0;
      vector<double> cutFlowVectorATLAS_250_0;
      double xsecATLAS_130_0;
      double xsecATLAS_250_0;
      vector<string> cutFlowVector_str;
      size_t NCUTS;

      Perf_Plot* plots_2bjets;	
      Perf_Plot* plots_mbb;	
      Perf_Plot* plots_HEPmct;	
      Perf_Plot* plots_HEPmt;	
      Perf_Plot* plots_HEPnbj;	
      Perf_Plot* plots_HEPmbb;	
      ofstream cutflowFile;
      string analysisRunName;

    public:

      struct particleComparison {
  	bool operator() (HEPUtils::Particle* i,HEPUtils::Particle* j) {return (i->pT()>j->pT());}
      }	compareParticlePt;

      struct jetComparison {
  	bool operator() (HEPUtils::Jet* i,HEPUtils::Jet* j) {return (i->pT()>j->pT());}
      }	compareJetPt;

      Analysis_ATLAS_8TeV_1LEPbb_20invfb() {

        _numSRA=0;
        _numSRB=0;

        NCUTS=8;
        set_luminosity(20.3);
        xsecATLAS_130_0=4240.;
        xsecATLAS_250_0=320.;

        for (size_t i=0;i<NCUTS;i++){
          cutFlowVector.push_back(0);
          cutFlowVectorATLAS_130_0.push_back(0);
          cutFlowVectorATLAS_250_0.push_back(0);
          cutFlowVector_str.push_back("");
        }

	analysisRunName = "ATLAS_8TeV_1LEPbb_20invfb_130_0";

	vector<const char*> variablesNames = {"met","mct","mbb","mt","j0pt","lpt","nbj","j1pt","j0eta","j1eta","jjdeltaR"};
	plots_2bjets = new Perf_Plot(analysisRunName+"_2bjets", &variablesNames);
	plots_mbb = new Perf_Plot(analysisRunName+"_mbb", &variablesNames);
	plots_HEPmct = new Perf_Plot(analysisRunName+"_HEPmct", &variablesNames);
	plots_HEPmt = new Perf_Plot(analysisRunName+"_HEPmt", &variablesNames);
	plots_HEPnbj = new Perf_Plot(analysisRunName+"_HEPnbj", &variablesNames);
	plots_HEPmbb = new Perf_Plot(analysisRunName+"_HEPmbb", &variablesNames);

      }


      void analyze(const HEPUtils::Event* event) {
	  
	HEPUtilsAnalysis::analyze(event);
        double met = event->met();

        // Baseline objects
        vector<HEPUtils::Particle*> baselineElectrons;
        for (HEPUtils::Particle* electron : event->electrons()) {
          if (electron->pT()>10. && electron->abseta()<2.47)baselineElectrons.push_back(electron);
        }
	ATLAS::applyMediumIDElectronSelection(baselineElectrons);

        vector<HEPUtils::Particle*> baselineMuons;
        for (HEPUtils::Particle* muon : event->muons()) {
          if (muon->pT()>10. && muon->abseta()<2.4)baselineMuons.push_back(muon);
        }

        vector<HEPUtils::Jet*> baselineJets;
	// Matthias jet smearing implemented roughly from 
	// https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2015-017/
        //const vector<double>  b3 = {0,50.,70.,100.,150.,200.,1000.,10000.};
        //const vector<double> cJets = {0.145,0.115,0.095,0.075,0.06,0.04,0.01};
	//static HEPUtils::BinnedFn2D<double> _resJets2D(a,b3,cJets);                                                                                                                     
        for (HEPUtils::Jet* jet : event->jets()) {
	  //const double _resJets = _resJets2D.get_at(jet->abseta(), jet->pT());
	  //std::random_device rd;
	  //std::mt19937 gen(rd());
	  //std::normal_distribution<> d(1., _resJets);
	  // Smear by a Gaussian centered on 1 with width given by the (fractional) resolution                                                                                    
	  //double smear_factor = d(gen);
	  /// @todo Is this the best way to smear? Should we preserve the mean jet energy, or pT, or direction?                                                                   
	  //jet->set_mom(HEPUtils::P4::mkXYZM(jet->mom().px()*smear_factor, jet->mom().py()*smear_factor, jet->mom().pz()*smear_factor, jet->mass()));
          //if (jet->pT()>20. && fabs(jet->eta())<4.5) baselineJets.push_back(jet);
          if (jet->pT()>20. && jet->abseta()<4.5)baselineJets.push_back(jet);
        }

        //Overlap Removal
        vector<HEPUtils::Particle*> overlapElectrons;
	vector<HEPUtils::Particle*> overlapMuons;
	vector<HEPUtils::Jet*> overlapJets;

        vector<size_t> overlapEl;
        for (size_t iEl1=0;iEl1<baselineElectrons.size();iEl1++) {
          for (size_t iEl2=0;iEl2<baselineElectrons.size();iEl2++) {
            if (baselineElectrons.at(iEl1)->mom().deltaR_eta(baselineElectrons.at(iEl2)->mom())<0.1) {
                if (baselineElectrons.at(iEl1)->pT()>baselineElectrons.at(iEl2)->pT())overlapEl.push_back(iEl2);
                if (baselineElectrons.at(iEl1)->pT()<baselineElectrons.at(iEl2)->pT())overlapEl.push_back(iEl1);
            }
          }
        }
        for (size_t iO=0;iO<overlapEl.size();iO++)baselineElectrons.erase(baselineElectrons.begin()+overlapEl.at(iO));

        for (size_t iJet=0;iJet<baselineJets.size();iJet++) {
          bool overlap=false;
          for (size_t iEl=0;iEl<baselineElectrons.size();iEl++) {
            if (fabs(baselineElectrons.at(iEl)->mom().deltaR_eta(baselineJets.at(iJet)->mom()))<0.2)overlap=true;
          }
          if (!overlap)overlapJets.push_back(baselineJets.at(iJet));
        }

        for (size_t iEl=0;iEl<baselineElectrons.size();iEl++) {
          bool overlap=false;
          for (size_t iJet=0;iJet<overlapJets.size();iJet++) {
            if (fabs(baselineElectrons.at(iEl)->mom().deltaR_eta(overlapJets.at(iJet)->mom()))<0.4)overlap=true;
          }
          if (!overlap)overlapElectrons.push_back(baselineElectrons.at(iEl));
        }

        for (size_t iMu=0;iMu<baselineMuons.size();iMu++) {
          bool overlap=false;
          for (size_t iJet=0;iJet<overlapJets.size();iJet++) {
            if (fabs(baselineMuons.at(iMu)->mom().deltaR_eta(overlapJets.at(iJet)->mom()))<0.4)overlap=true;
          }
          if (!overlap)overlapMuons.push_back(baselineMuons.at(iMu));
        }

	//Signal Objects
        vector<HEPUtils::Particle*> signalLeptons;
        vector<HEPUtils::Particle*> signalElectrons;
        vector<HEPUtils::Particle*> signalMuons;
        vector<HEPUtils::Jet*> signalJets;   
	vector<HEPUtils::Jet*> signalBJets;

	for (size_t iEl=0;iEl<overlapElectrons.size();iEl++) {
	  if (overlapElectrons.at(iEl)->pT()>25.)signalElectrons.push_back(overlapElectrons.at(iEl));
        }
	ATLAS::applyTightIDElectronSelection(signalElectrons);
        
        for (size_t iMu=0;iMu<overlapMuons.size();iMu++) {
          if (overlapMuons.at(iMu)->pT()>25.)signalMuons.push_back(overlapMuons.at(iMu)); 
        } 
	       
	const vector<double> aBJet = {0,2.1,10.};
	const vector<double> bBJet = {0,30.,40.,50.,70.,10000.};
        const vector<double> cBJet={0.54,0.63,0.67,0.7,0.75,0.35,0.42,0.44,0.46,0.49};
        HEPUtils::BinnedFn2D<double> _eff2dBJet(aBJet,bBJet,cBJet);
        for (size_t iJet=0;iJet<overlapJets.size();iJet++) {
          if (overlapJets.at(iJet)->pT()>25. && overlapJets.at(iJet)->abseta()<2.40) {
	    signalJets.push_back(overlapJets.at(iJet));  
            bool hasTag=has_tag(_eff2dBJet, overlapJets.at(iJet)->eta(), overlapJets.at(iJet)->pT());
	    if (overlapJets.at(iJet)->btag() && hasTag)signalBJets.push_back(overlapJets.at(iJet));             
          }
	}

        signalLeptons=signalElectrons;
        signalLeptons.insert(signalLeptons.end(),signalMuons.begin(),signalMuons.end());
        int nSignalLeptons=signalLeptons.size();
        int nBaselineLeptons=overlapElectrons.size()+overlapMuons.size();
	int nSignalElectrons=signalElectrons.size();
	int nSignalMuons=signalMuons.size();
        int nSignalJets=signalJets.size();
        int nSignalBJets=signalBJets.size();
	sort(signalJets.begin(), signalJets.end(), compareJetPt);
	sort(signalLeptons.begin(), signalLeptons.end(), compareParticlePt);

	//Variables	
	double mT=0; 
        double mCT=0;
        double mbb=0;
        bool leadingBJets=isLeadingBJets(signalJets, signalBJets);
	bool lepton_overlap=true;
        bool preselection=false; 

        const vector<double>  aLep = {0,10.};
        const vector<double>  bLep = {0,10000.};
        const vector<double> cEl = {0.95};
        const vector<double> cMu1 = {0.7};
        const vector<double> cMu2 = {0.85};
        HEPUtils::BinnedFn2D<double> _eff2dMu1(aLep,bLep,cMu1);
        HEPUtils::BinnedFn2D<double> _eff2dMu2(aLep,bLep,cMu2);
        HEPUtils::BinnedFn2D<double> _eff2dEl(aLep,bLep,cEl);

	for (size_t iEl=0;iEl<overlapElectrons.size();iEl++) {
	  for (size_t iMu=0;iMu<overlapMuons.size();iMu++) {
	    if(fabs(overlapElectrons.at(iEl)->mom().deltaR_eta(overlapMuons.at(iMu)->mom()))<0.1)lepton_overlap=false;
	  }
	}
        for (size_t iMu1=0;iMu1<overlapMuons.size();iMu1++) {
          for (size_t iMu2=0;iMu2<overlapMuons.size();iMu2++) {
            if(fabs(overlapMuons.at(iMu1)->mom().deltaR_eta(overlapMuons.at(iMu2)->mom()))<0.05)lepton_overlap=false;
          }
        }
        if (lepton_overlap && nSignalLeptons==1 && nBaselineLeptons==1 && (nSignalJets==2 || nSignalJets==3) && leadingBJets) {
	  if (nSignalMuons==1) {	
            bool hasTrig1=has_tag(_eff2dMu1,signalMuons.at(0)->eta(),signalMuons.at(0)->pT());
            bool hasTrig2=has_tag(_eff2dMu2,signalMuons.at(0)->eta(),signalMuons.at(0)->pT());
	    if (signalMuons.at(0)->abseta()<1.05 && hasTrig1)preselection=true;
	    if (signalMuons.at(0)->abseta()>1.05 && hasTrig2)preselection=true;
	  }
	  if (nSignalElectrons==1) {	
            bool hasTrig=has_tag(_eff2dEl,signalElectrons.at(0)->eta(),signalElectrons.at(0)->pT());
	    if (hasTrig)preselection=true;
	  }
	}

	if (nSignalLeptons)mT=sqrt(2*signalLeptons.at(0)->pT()*met*(1-cos(signalLeptons.at(0)->phi()-event->missingmom().phi())));      
	if (nSignalJets>1) {
          mCT=sqrt(2*signalJets.at(0)->pT()*signalJets.at(1)->pT()*(1+cos(signalJets.at(0)->phi()-signalJets.at(1)->phi())));
          mbb=(signalJets.at(0)->mom()+signalJets.at(1)->mom()).m(); 
	}	

	bool SRA=false;
	bool SRB=false;
        if (preselection && nSignalBJets==2 && met>100. && mCT>160. && mbb>105. && mbb<135.) {
          if (mT>100. && mT<130.) {
            _numSRA++;
            SRA=true;
          }
          if (mT>130.) {
            _numSRB++;
	    SRB=true;   
          }
        }                      

	if (preselection) {
	  vector<double> variables={met, mCT, mbb, mT, signalJets.at(0)->pT(), signalLeptons.at(0)->pT(), nSignalBJets, signalJets.at(1)->pT(),signalJets.at(0)->eta(), signalJets.at(1)->eta(), signalJets.at(0)->mom().deltaR_eta(signalJets.at(1)->mom())};
	  if (met>50. && mT>40. && mbb>40. && nSignalBJets==2)plots_2bjets->fill(&variables);
	  if (met>50. && mT>40. && mbb>40. && nSignalBJets==2 && met>100. && mCT>160. && mT>100. && mbb>45. && mbb<195.)plots_mbb->fill(&variables);
          if (nSignalBJets==2 && met>100. && mT>100. && mbb>45. && mbb<195. && (mbb<105. || mbb>135.))plots_HEPmct->fill(&variables);
          if (nSignalBJets==2 && met>100. && mCT>160. && mbb>45. && mbb<195. && (mbb<105. || mbb>135.))plots_HEPmt->fill(&variables);
          if (nSignalBJets==2 && met>100 && mCT>160. && mT>100)plots_HEPmbb->fill(&variables);
          if (met>100. && mCT>160. && mT>100. && mbb>105. && mbb<135.)plots_HEPnbj->fill(&variables);
	}

	cutFlowVector_str[1] = "Lepton + 2 b-jets";
        cutFlowVector_str[2] = "$E_{T}^{miss} > 100 GeV$";
        cutFlowVector_str[3] = "$m_{CT} > 160 GeV$";
        cutFlowVector_str[4] = "$m_{T} > 100 GeV$";
        cutFlowVector_str[5] = "$45 GeV < m_{bb} < 195 GeV$";
        cutFlowVector_str[6] = "SRA";
        cutFlowVector_str[7] = "SRB";

	cutFlowVectorATLAS_130_0[0] = 100000;
	cutFlowVectorATLAS_130_0[1] = 531.1;
        cutFlowVectorATLAS_130_0[2] = 163.7;
        cutFlowVectorATLAS_130_0[3] = 70.4;
        cutFlowVectorATLAS_130_0[4] = 9.7;
        cutFlowVectorATLAS_130_0[5] = 9.6;
        cutFlowVectorATLAS_130_0[6] = 7.2;
        cutFlowVectorATLAS_130_0[7] = 0.3;

	cutFlowVectorATLAS_250_0[0] = 99000;
	cutFlowVectorATLAS_250_0[1] = 71.3;
        cutFlowVectorATLAS_250_0[2] = 45.2;
        cutFlowVectorATLAS_250_0[3] = 15.0;
        cutFlowVectorATLAS_250_0[4] = 8.1;
        cutFlowVectorATLAS_250_0[5] = 8.0;
        cutFlowVectorATLAS_250_0[6] = 1.3;
        cutFlowVectorATLAS_250_0[7] = 4.4;

        for (size_t j=0;j<NCUTS;j++){
          if(
             (j==0) ||
             
	     (j==1 && preselection && nSignalBJets==2) ||
             
             (j==2 && preselection && nSignalBJets==2 && met>100.) ||

             (j==3 && preselection && nSignalBJets==2 && met>100. && mCT>160.) ||
 
             (j==4 && preselection && nSignalBJets==2 && met>100. && mCT>160. && mT>100.) ||
             
	     (j==5 && preselection && nSignalBJets==2 && met>100. && mCT>160. && mT>100. && mbb>45. && mbb<195.) ||  

             (j==6 && SRA) ||  
            
	     (j==7 && SRB) ){

            cutFlowVector[j]++;

          }
        }
      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        HEPUtilsAnalysis::add(other);

        Analysis_ATLAS_8TeV_1LEPbb_20invfb* specificOther
                = dynamic_cast<Analysis_ATLAS_8TeV_1LEPbb_20invfb*>(other);

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

        if (analysisRunName.find("250_0") != string::npos) {
          cutflowFile<<"\\begin{table}[H] \n\\caption{$\\tilde{\\chi}_{1}^{\\pm}\\tilde{\\chi}_{2}^{0}$ decay via $W/h$, $[\\tilde{\\chi}_{1}^{\\pm}\\tilde{\\chi}_{2}^{0},\\tilde{\\chi}_{1}^{0}]: [250,0] [GeV]$} \n\\makebox[\\linewidth]{ \n\\renewcommand{\\arraystretch}{0.4} \n\\begin{tabular}{c c c c c} \n\\hline"<<endl;
	  cutflowFile<<"& ATLAS & GAMBIT & GAMBIT/ATLAS & $\\sigma$-corrected GAMBIT/ATLAS \\\\ \\hline"<<endl;
	  cutflowFile<<"$\\sigma (pp\\to \\tilde{\\chi}_{1}^{\\pm}, \\tilde{\\chi}_{2}^{0})$ &"<<setprecision(4)<<xsecATLAS_250_0<<" $fb$ &"<<setprecision(4)<<xsec()<<"$fb$ &"<<setprecision(4)<<xsec()/xsecATLAS_250_0<<" & 1\\\\"<<endl;
	  cutflowFile<<"Generated Events &"<<setprecision(4)<<cutFlowVectorATLAS_250_0[0]<<"&"<<setprecision(4)<<cutFlowVector[0]<<"& - & -\\\\ \\hline"<<endl;
	  cutflowFile<<"\\multicolumn{5}{c}{Expected events at 20.3 $fb^{-1}$} \\\\ \\hline"<<endl;
          for (size_t i=1; i<NCUTS; i++) {
            cutflowFile<<cutFlowVector_str[i]<<"&"<<setprecision(4)<<cutFlowVectorATLAS_250_0[i]<<"&"<<setprecision(4)<<cutFlowVector[i]*xsec_per_event()*luminosity()<<"&"<<setprecision(4)<<cutFlowVector[i]*xsec_per_event()*luminosity()/cutFlowVectorATLAS_250_0[i]<<"&"<<setprecision(4)<<(xsecATLAS_250_0/xsec())*cutFlowVector[i]*xsec_per_event()*luminosity()/cutFlowVectorATLAS_250_0[i]<<"\\\\"<< endl;
          }
	  cutflowFile<<"\\hline \\multicolumn{5}{c}{Percentage (\\%)} \\\\ \\hline"<<endl;
          for (size_t i=1; i<NCUTS; i++) {
            cutflowFile<<cutFlowVector_str[i]<<"&"<<setprecision(4)<<cutFlowVectorATLAS_250_0[i]*100./cutFlowVectorATLAS_250_0[1]<<"&"<<setprecision(4)<<cutFlowVector[i]*100./cutFlowVector[1]<<"& - & -\\\\"<< endl;
          }
          cutflowFile<<"\\end{tabular} \n} \n\\end{table}"<<endl;
	}
        if (analysisRunName.find("130_0") != string::npos) {
          cutflowFile<<"\\begin{table}[H] \n\\caption{$\\tilde{\\chi}_{1}^{\\pm}\\tilde{\\chi}_{2}^{0}$ decay via $W/h$, $[\\tilde{\\chi}_{1}^{\\pm}\\tilde{\\chi}_{2}^{0},\\tilde{\\chi}_{1}^{0}]: [130,0] [GeV]$} \n\\makebox[\\linewidth]{ \n\\renewcommand{\\arraystretch}{0.4} \n\\begin{tabular}{c c c c c} \n\\hline"<<endl;
	  cutflowFile<<"& ATLAS & GAMBIT & GAMBIT/ATLAS & $\\sigma$-corrected GAMBIT/ATLAS \\\\ \\hline"<<endl;
	  cutflowFile<<"$\\sigma (pp\\to \\tilde{\\chi}_{1}^{\\pm}, \\tilde{\\chi}_{2}^{0})$ &"<<setprecision(4)<<xsecATLAS_130_0<<" $fb$ &"<<setprecision(4)<<xsec()<<"$fb$ &"<<setprecision(4)<<xsec()/xsecATLAS_130_0<<" & 1\\\\"<<endl;
	  cutflowFile<<"Generated Events &"<<setprecision(4)<<cutFlowVectorATLAS_130_0[0]<<"&"<<setprecision(4)<<cutFlowVector[0]<<"& - & -\\\\ \\hline"<<endl;
	  cutflowFile<<"\\multicolumn{5}{c}{Expected events at 20.3 $fb^{-1}$} \\\\ \\hline"<<endl;
          for (size_t i=1; i<NCUTS; i++) {
            cutflowFile<<cutFlowVector_str[i]<<"&"<<setprecision(4)<<cutFlowVectorATLAS_130_0[i]<<"&"<<setprecision(4)<<cutFlowVector[i]*xsec_per_event()*luminosity()<<"&"<<setprecision(4)<<cutFlowVector[i]*xsec_per_event()*luminosity()/cutFlowVectorATLAS_130_0[i]<<"&"<<setprecision(4)<<(xsecATLAS_130_0/xsec())*cutFlowVector[i]*xsec_per_event()*luminosity()/cutFlowVectorATLAS_130_0[i]<<"\\\\"<< endl;
          }
	  cutflowFile<<"\\hline \\multicolumn{5}{c}{Percentage (\\%)} \\\\ \\hline"<<endl;
          for (size_t i=1; i<NCUTS; i++) {
            cutflowFile<<cutFlowVector_str[i]<<"&"<<setprecision(4)<<cutFlowVectorATLAS_130_0[i]*100./cutFlowVectorATLAS_130_0[1]<<"&"<<setprecision(4)<<cutFlowVector[i]*100./cutFlowVector[1]<<"& - & -\\\\"<< endl;
          }
          cutflowFile<<"\\end{tabular} \n} \n\\end{table}"<<endl;
	}
	cutflowFile.close();

	plots_2bjets->createFile(luminosity(),xsec_per_event());
	plots_mbb->createFile(luminosity(),xsec_per_event());
	plots_HEPmct->createFile(luminosity(),xsec_per_event());
	plots_HEPmt->createFile(luminosity(),xsec_per_event());
	plots_HEPnbj->createFile(luminosity(),xsec_per_event());
	plots_HEPmbb->createFile(luminosity(),xsec_per_event());

	SignalRegionData results_SRA;
        results_SRA.analysis_name = "Analysis_ATLAS_8TeV_1LEPbb_20invfb";
        results_SRA.sr_label = "SRA";
        results_SRA.n_observed = 4.;
        results_SRA.n_background = 5.69; 
        results_SRA.background_sys = 1.10;
        results_SRA.signal_sys = 0.; 
        results_SRA.n_signal = _numSRA;
        add_result(results_SRA);

        SignalRegionData results_SRB;
        results_SRB.analysis_name = "Analysis_ATLAS_8TeV_1LEPbb_20invfb";
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
    DEFINE_ANALYSIS_FACTORY(ATLAS_8TeV_1LEPbb_20invfb)


  }
}
