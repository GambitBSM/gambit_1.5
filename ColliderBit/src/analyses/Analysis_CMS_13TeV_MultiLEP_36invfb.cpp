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

    class Analysis_CMS_13TeV_MultiLEP_36invfb : public HEPUtilsAnalysis {
    private:

      // Numbers passing cuts
      double _numSR1, _numSR2, _numSR3, _numSR4, _numSR5, _numSR6, _numSR7, _numSR8; 
      vector<int> cutFlowVector1, cutFlowVector2, cutFlowVector3, cutFlowVector4;
      vector<string> cutFlowVector_str1, cutFlowVector_str2, cutFlowVector_str3, cutFlowVector_str4;
      size_t NCUTS1, NCUTS2, NCUTS3, NCUTS4;

      Perf_Plot* plots;
      ofstream cutflowFile;
      string analysisRunName;

    public:

      Analysis_CMS_13TeV_MultiLEP_36invfb() {

        _numSR1=0;
	_numSR2=0;
	_numSR3=0;
	_numSR4=0;
	_numSR5=0; 
	_numSR6=0;
	_numSR7=0; 
	_numSR8=0;

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

        analysisRunName = "CMS_13TeV_MultiLEP_36invfb_test";
        vector<const char*> variables = {"met","mT","mT2"};
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
	  if (electron->pT()>=10. &&fabs(electron->eta())<2.5) {
            baselineElectrons.push_back(electron);
	    baselineLeptons.push_back(electron);
          }
	}

        vector<HEPUtils::Particle*> baselineMuons;
        for (HEPUtils::Particle* muon : event->muons()) {
	  if (muon->pT()>=10. &&fabs(muon->eta())<2.4) {
            baselineMuons.push_back(muon);
	    baselineLeptons.push_back(muon);
          }
	}

        vector<HEPUtils::Particle*> baselineTaus;
        for (HEPUtils::Particle* tau : event->taus()) {
          if (tau->pT()>=20. &&fabs(tau->eta())<2.3) {
            baselineTaus.push_back(tau);
            baselineLeptons.push_back(tau);
          }
        }

        vector<HEPUtils::Jet*> baselineJets;
        for (HEPUtils::Jet* jet : event->jets()) {
          if (jet->pT()>=25. &&fabs(jet->eta())<2.4) {
            baselineJets.push_back(jet);
	  }
        }

        // Signal objects
        vector<HEPUtils::Particle*> signalLeptons=baselineLeptons;
        vector<HEPUtils::Particle*> signalElectrons=baselineElectrons;
        vector<HEPUtils::Particle*> signalMuons=baselineMuons;
        vector<HEPUtils::Particle*> signalTaus=baselineTaus;
        vector<HEPUtils::Jet*> signalJets;   
        vector<HEPUtils::Particle*> signalLightLeptons=signalElectrons;
	signalLightLeptons.insert(signalLightLeptons.end(),signalMuons.begin(),signalMuons.end());
        sort(signalLightLeptons.begin(),signalLightLeptons.end(),comparePt);
        sort(signalLeptons.begin(),signalLeptons.end(),comparePt);

	//Jets (overlap removal) 
        for (size_t iJet=0;iJet<baselineJets.size();iJet++) {
	  bool overlap=false;            
	  for (size_t iLe=0;iLe<baselineLeptons.size();iLe++) {
	    if (fabs(baselineLeptons.at(iLe)->mom().deltaR_eta(baselineJets.at(iJet)->mom()))<0.4)overlap=true;
	  }
	  if (!overlap)signalJets.push_back(baselineJets.at(iJet));
	}

        //Variable definitions
        int nSignalLeptons = signalLeptons.size();
        int nSignalJets = signalJets.size();
	int nSignalElectrons = signalElectrons.size();
	int nSignalMuons = signalMuons.size();
	int nSignalTaus = signalTaus.size(); 
	int nSignalLightLeptons = signalLightLeptons.size();
        
	//Preselection
        bool preselection=false; 
 
	//Bjet veto
	bool bjet_veto=true;
	for (size_t iJet=0;iJet<signalJets.size();iJet++) {
	  if (signalJets.at(iJet)->btag())bjet_veto=false;
	}

	//Low-mass veto
	bool low_mass_veto=true;
	bool conversion_veto=true;
	
	vector<vector<HEPUtils::Particle*>> OSSFpair_cont = getOSSFpair(signalLeptons);
	if (OSSFpair_cont.size()>0) {
	  for (size_t iPa=0;iPa<OSSFpair_cont.size();iPa++) {
	    if (OSSFpair_cont.at(iPa).at(0)->mass()+OSSFpair_cont.at(iPa).at(1)->mass()<12)low_mass_veto=false;
	    if (nSignalLeptons==3) {
	      if (OSSFpair_cont.at(iPa).at(0)->abspid()!=15 && abs(signalLeptons.at(0)->mass()+signalLeptons.at(1)->mass()+signalLeptons.at(2)->mass()-91.2)<15)conversion_veto=false;
	    }	
	    if (abs(OSSFpair_cont.at(iPa).at(0)->mass()+OSSFpair_cont.at(iPa).at(1)->mass()-91.2)<15)conversion_veto=false;
	  }
	}
        
	if (bjet_veto && low_mass_veto)preselection=true;

        //Signal regions

	//2 same-sign leptons
	bool two_ss_leptons=false;	
	if (preselection && nSignalLeptons==2 && nSignalTaus==0) {
	  if (signalLeptons.at(0)->pid()*signalLeptons.at(1)->pid()>0 && met>60 && conversion_veto) {
	    if (signalLeptons.at(0)->abspid()==11 && signalLeptons.at(0)->pT()>25) {
	      if (signalLeptons.at(1)->abspid()==11 && signalLeptons.at(1)->pT()>15)two_ss_leptons=true;
	      if (signalLeptons.at(1)->abspid()==13 && signalLeptons.at(1)->pT()>10)two_ss_leptons=true;
  	    }  	
            if (signalLeptons.at(0)->abspid()==13 && signalLeptons.at(0)->pT()>20) {
              if (signalLeptons.at(1)->abspid()==11 && signalLeptons.at(1)->pT()>15)two_ss_leptons=true;
              if (signalLeptons.at(1)->abspid()==13 && signalLeptons.at(1)->pT()>10)two_ss_leptons=true; 
	    }
	  }
	}

        bool three_more_leptons=false;
        if (preselection && met>50 && nSignalLeptons>=3 && nSignalTaus<=2 && conversion_veto) {
          
	  if (nSignalTaus==2) {
            //2 taus + 1 electron
            if (nSignalElectrons==1) {
	      if (signalElectrons.at(0)->pT()>30.)three_more_leptons=true;
	    }
	    //2 taus + 1 muon
            if (nSignalMuons==1) { 
	      if (signalMuons.at(0)->pT()>25.)three_more_leptons=true;
	    }
	    //2 taus + electron(s)/muon(s)
	    if (nSignalTaus<2) {
              if (signalLightLeptons.at(0)->abspid()==11 && signalLightLeptons.at(0)->pT()>30) {
                if (signalLightLeptons.at(1)->abspid()==11 && signalLightLeptons.at(1)->pT()>15)three_more_leptons=true;
                if (signalLightLeptons.at(1)->abspid()==13 && signalLightLeptons.at(1)->pT()>10)three_more_leptons=true;
              }
              if (signalLightLeptons.at(0)->abspid()==13 && signalLightLeptons.at(0)->pT()>25) {
                if (signalLightLeptons.at(1)->abspid()==11 && signalLightLeptons.at(1)->pT()>15)three_more_leptons=true;
                if (signalLightLeptons.at(1)->abspid()==13 && signalLightLeptons.at(1)->pT()>10)three_more_leptons=true;
	      }
	    }
            for (size_t iLe=0;iLe<signalLeptons.size();iLe++) {
              if (fabs(signalLeptons.at(iLe)->eta())>=2.1)three_more_leptons=false;
            }
          }

	  else {
	    //0 or 1 tau + 1 (leading) muon + electron(s)
	    if (signalLeptons.at(0)->abspid()==13 && nSignalMuons == 1) {
              if (signalLightLeptons.at(0)->pT()>25 && signalLightLeptons.at(1)->pT()>15)three_more_leptons=true;
            }
	    //Every other possibility
	    else {
	      if (signalLightLeptons.at(0)->abspid()==11 && signalLightLeptons.at(0)->pT()>25) {
                if (signalLightLeptons.at(1)->abspid()==11 && signalLightLeptons.at(1)->pT()>15)three_more_leptons=true;
                if (signalLightLeptons.at(1)->abspid()==13 && signalLightLeptons.at(1)->pT()>10)three_more_leptons=true;
              }
              if (signalLightLeptons.at(0)->abspid()==13 && signalLightLeptons.at(0)->pT()>20) {
                if (signalLightLeptons.at(1)->abspid()==11 && signalLightLeptons.at(1)->pT()>15)three_more_leptons=true;
                if (signalLightLeptons.at(1)->abspid()==13 && signalLightLeptons.at(1)->pT()>10)three_more_leptons=true;
              }
            }
	  }
	
	}

	double pT_ll=0;
	double mT=0;
	double mT2=0;
	if (two_ss_leptons) {
         //should I use leading lepton here or not for mT???
	  pT_ll = signalLeptons.at(0)->pT()+signalLeptons.at(1)->pT();
	  mT = sqrt(2*signalLeptons.at(1)->pT()*met*(1-cos(signalLeptons.at(1)->phi()-event->missingmom().phi()))); 
	  if (nSignalJets==0 && met>140 && mT>100)_numSR1++;
	  if (nSignalJets==1 && met>200 && mT<100 && pT_ll<100)_numSR2++;
	}
	
	if (three_more_leptons) {
	  //0 taus + 3 electron(s)/muon(s)
	  if (nSignalLightLeptons==3 && nSignalTaus==0) {
	    //should I use leading lepton here or not for mT???
	    mT = sqrt(2*signalLeptons.at(2)->pT()*met*(1-cos(signalLeptons.at(2)->phi()-event->missingmom().phi())));
	    if (mT>120 && met>200)_numSR3++;
	    if (met>250)_numSR4++;
	  }
	  //1 tau + 2 electron(s)/muon(s)
	  if (nSignalLightLeptons==2 && nSignalTaus==1) {
              
	    //MT2 for leading lepton + tau
	    double pLep1[3] = {signalLightLeptons.at(0)->mass(), signalLightLeptons.at(0)->mom().px(), signalLightLeptons.at(0)->mom().py()};
            double pTau[3] = {signalTaus.at(0)->mass(), signalTaus.at(0)->mom().px(), signalTaus.at(0)->mom().py()};
            double pMiss[3] = {0., event->missingmom().px(), event->missingmom().py() };
            double mn = 0.;

            mt2_bisect::mt2 mt2_calc;
            mt2_calc.set_momenta(pLep1,pTau,pMiss);
            mt2_calc.set_mn(mn);
            mT2 = mt2_calc.get_mt2();

            if (mT2>50 && met>200)_numSR5++;
	  }
	  //2 taus + 1 electron/muon
	  if (nSignalLightLeptons==1 && nSignalTaus==2) {
	    
	    //MT2 for lepton + leading tau
	    sort(signalTaus.begin(),signalTaus.end(),comparePt);
	    double pLep1[3] = {signalLightLeptons.at(0)->mass(), signalLightLeptons.at(0)->mom().px(), signalLightLeptons.at(0)->mom().py()};
            double pTau[3] = {signalTaus.at(0)->mass(), signalTaus.at(0)->mom().px(), signalTaus.at(0)->mom().py()};
            double pMiss[3] = {0., event->missingmom().px(), event->missingmom().py()};
            double mn = 0.;

            mt2_bisect::mt2 mt2_calc;
            mt2_calc.set_momenta(pLep1,pTau,pMiss);
            mt2_calc.set_mn(mn);
            mT2 = mt2_calc.get_mt2();

	    if (mT2>50 && met>200)_numSR6++;
	    if (met>75)_numSR7++;
	  }
	  //>3 leptons
	  if (nSignalLeptons>3 && met>200)_numSR8++;
	}

        if (preselection) {
          vector<double> variables = {met, mT, mT2};
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

	bool ISRjet=false;
	vector<HEPUtils::Jet*> ISRjets;
	for (size_t iJet=0;iJet<signalJets.size();iJet++) {
	  if (signalJets.at(iJet)->pT()>40)ISRjets.push_back(signalJets.at(iJet));
	}
	if (ISRjets.size()==0 || ISRjets.size()==1)ISRjet=true;

	double mll=0;
	if (OSSFpair_cont.size()!=0) {
	  mll = get_mll(OSSFpair_cont);
	}
	else {
	  vector<vector<HEPUtils::Particle*>> OSpair_cont = getOSpair(signalLeptons);
	  mll = get_mll(OSpair_cont);
	}
   
        for (size_t j=0;j<NCUTS1;j++){
          if(
             (j==0) ||

	     (j==1 && nSignalLightLeptons==2) ||

	     (j==2 && nSignalLightLeptons==2 && signalLightLeptons.at(0)->pid()*signalLightLeptons.at(1)->pid()>0) ||

	     (j==3 && nSignalLightLeptons==2 && signalLightLeptons.at(0)->pid()*signalLightLeptons.at(1)->pid()>0 && nSignalLeptons==2) ||

	     (j==4 && nSignalLightLeptons==2 && signalLightLeptons.at(0)->pid()*signalLightLeptons.at(1)->pid()>0 && nSignalLeptons==2 && low_mass_veto) ||

	     (j==5 && nSignalLightLeptons==2 && signalLightLeptons.at(0)->pid()*signalLightLeptons.at(1)->pid()>0 && nSignalLeptons==2 && low_mass_veto && bjet_veto)  ||
	     
	     (j==6 && nSignalLightLeptons==2 && signalLightLeptons.at(0)->pid()*signalLightLeptons.at(1)->pid()>0 && nSignalLeptons==2 && low_mass_veto && bjet_veto && met>60) ||            

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

        Analysis_CMS_13TeV_MultiLEP_36invfb* specificOther
                = dynamic_cast<Analysis_CMS_13TeV_MultiLEP_36invfb*>(other);

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
        results_SR1.analysis_name = "Analysis_CMS_13TeV_MultiLEP_36invfb";
        results_SR1.sr_label = "SR1";
        results_SR1.n_observed = 13.;
        results_SR1.n_background = 12.; 
        results_SR1.background_sys = 3.;
        results_SR1.signal_sys = 0.; 
        results_SR1.n_signal = _numSR1;
	add_result(results_SR1);

	SignalRegionData results_SR2;
        results_SR2.analysis_name = "Analysis_CMS_13TeV_MultiLEP_36invfb";
        results_SR2.sr_label = "SR1";
        results_SR2.n_observed = 18.;
        results_SR2.n_background = 18.;
        results_SR2.background_sys = 4.;
        results_SR2.signal_sys = 0.;
        results_SR2.n_signal = _numSR2;
        add_result(results_SR2);

        SignalRegionData results_SR3;
        results_SR3.analysis_name = "Analysis_CMS_13TeV_MultiLEP_36invfb";
        results_SR3.sr_label = "SR3";
        results_SR3.n_observed = 19.;
        results_SR3.n_background = 19.;
        results_SR3.background_sys = 4.;
        results_SR3.signal_sys = 0.;
        results_SR3.n_signal = _numSR3;
        add_result(results_SR3);

        SignalRegionData results_SR4;
        results_SR4.analysis_name = "Analysis_CMS_13TeV_MultiLEP_36invfb";
        results_SR4.sr_label = "SR4";
        results_SR4.n_observed = 128.;
        results_SR4.n_background = 142.;
        results_SR4.background_sys = 34.;
        results_SR4.signal_sys = 0.;
        results_SR4.n_signal = _numSR4;
        add_result(results_SR4);

        SignalRegionData results_SR5;
        results_SR5.analysis_name = "Analysis_CMS_13TeV_MultiLEP_36invfb";
        results_SR5.sr_label = "SR5";
        results_SR5.n_observed = 18.;
        results_SR5.n_background = 22.;
        results_SR5.background_sys = 5.;
        results_SR5.signal_sys = 0.;
        results_SR5.n_signal = _numSR5;
        add_result(results_SR5);

        SignalRegionData results_SR6;
        results_SR6.analysis_name = "Analysis_CMS_13TeV_MultiLEP_36invfb";
        results_SR6.sr_label = "SR6";
        results_SR6.n_observed = 2;
        results_SR6.n_background = 1.2;
        results_SR6.background_sys = 0.6;
        results_SR6.signal_sys = 0.;
        results_SR6.n_signal = _numSR6;
        add_result(results_SR6);

        SignalRegionData results_SR7;
        results_SR7.analysis_name = "Analysis_CMS_13TeV_MultiLEP_36invfb";
        results_SR7.sr_label = "SR7";
        results_SR7.n_observed = 82.;
        results_SR7.n_background = 109.;
        results_SR7.background_sys = 28.;
        results_SR7.signal_sys = 0.;
        results_SR7.n_signal = _numSR7;
        add_result(results_SR7);

        SignalRegionData results_SR8;
        results_SR8.analysis_name = "Analysis_CMS_13TeV_MultiLEP_36invfb";
        results_SR8.sr_label = "SR8";
        results_SR8.n_observed = 166.;
        results_SR8.n_background = 197.;
        results_SR8.background_sys = 42.;
        results_SR8.signal_sys = 0.;
        results_SR8.n_signal = _numSR8;
        add_result(results_SR4);

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

      double get_mll(vector<vector<HEPUtils::Particle*>> pair_cont) { 
	double mll;
	vector<double> mll_container;
	for (size_t iPa=0;iPa<pair_cont.size();iPa++) {
	  mll_container.push_back(pair_cont.at(iPa).at(0)->mass()+pair_cont.at(iPa).at(1)->mass());
	}	  
	sort(mll_container.begin(),mll_container.end());
	if (mll_container.size()==1)mll=mll_container.at(0);
	if (mll_container.size()>1) {
 	  auto const it_lower = lower_bound(mll_container.begin(), mll_container.end(), 91.2);
 	  auto const it_upper = upper_bound(mll_container.begin(), mll_container.end(), 91.2);
	  if (it_lower==mll_container.end() || it_upper==mll_container.end())mll=*mll_container.end();
	  else {
	    double mll_lower_delta=fabs(91.2-*it_lower);
	    double mll_upper_delta=fabs(91.2-*it_upper);
	    if (mll_lower_delta<mll_upper_delta)mll=*it_lower;
	    if (mll_lower_delta>=mll_upper_delta)mll=*it_upper;
	  }
	}
	return mll;
      }

      struct ptComparison {
        bool operator() (HEPUtils::Particle* i,HEPUtils::Particle* j) {return (i->pT()>j->pT());}
      } comparePt;

    };


    // Factory fn
    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_MultiLEP_36invfb)


  }
}
