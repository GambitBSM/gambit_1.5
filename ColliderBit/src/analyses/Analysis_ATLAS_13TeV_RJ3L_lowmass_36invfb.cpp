#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>

#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "gambit/ColliderBit/mt2_bisect.h"

#include "RestFrames/RestFrames.hh"
#include "TLorentzVector.h"

using namespace std;

/* The ATLAS 13 TeV 3 lepton low mass recursive jigsaw search

   Based on code kindly supplied by Abhishek Sharma

   Note that use of ROOT is compulsory for the RestFrames package

   Based on: arXiv link N/A at present
  
   Code adapted by Martin White

   KNOWN ISSUES

   1) Need to check overlap removal step when the paper comes out. For now, have assumed it is the same as the stop analysis.

*/

namespace Gambit {
  namespace ColliderBit {

    bool sortByPT_RJ3L(HEPUtils::Jet* jet1, HEPUtils::Jet* jet2) { return (jet1->pT() > jet2->pT()); }
    //bool sortByMass(HEPUtils::Jet* jet1, HEPUtils::Jet* jet2) { return (jet1->mass() > jet2->mass()); }

    bool SortLeptons(const pair<TLorentzVector,int> lv1, const pair<TLorentzVector,int> lv2)
    //bool VariableConstruction::SortLeptons(const lep lv1, const lep lv2) 
    {
      return lv1.first.Pt() > lv2.first.Pt();
    }
    
    // Class to randomly select m elements from an n-d vector
    template<class BidiIter >
    BidiIter random_unique(BidiIter begin, BidiIter end, size_t num_random) {
      size_t left = std::distance(begin, end);
      while (num_random--) {
        BidiIter r = begin;
        std::advance(r, rand()%left);
        std::swap(*begin, *r);
        ++begin;
        --left;
      }
      return begin;
    }
    
    
    class Analysis_ATLAS_13TeV_RJ3L_lowmass_36invfb : public HEPUtilsAnalysis {
    private:

      // Numbers passing cuts
      int _numSR;

      vector<int> cutFlowVector;
      vector<string> cutFlowVector_str;
      int NCUTS; //=16;

      	// Recursive jigsaw objects (using RestFrames)

	RestFrames::LabRecoFrame*       LAB_3L;
	RestFrames::DecayRecoFrame*     C1N2_3L;
	RestFrames::DecayRecoFrame*     C1a_3L;
	RestFrames::DecayRecoFrame*     N2b_3L;
	
	RestFrames::DecayRecoFrame*     Wa_3L;
	RestFrames::DecayRecoFrame*     Zb_3L;
	
	RestFrames::VisibleRecoFrame*   L1a_3L;
	RestFrames::VisibleRecoFrame*   L1b_3L;
	RestFrames::VisibleRecoFrame*   L2b_3L;
	
	RestFrames::InvisibleRecoFrame* X1a_3L;
	RestFrames::InvisibleRecoFrame* X1b_3L;
	
	RestFrames::InvisibleGroup*    INV_3L;
	
	RestFrames::SetMassInvJigsaw*     X1_mass_3L;
	RestFrames::SetRapidityInvJigsaw* X1_eta_3L;
	
	RestFrames::ContraBoostInvJigsaw* X1X1_contra_3L;

      
      // Debug histos

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

      void LeptonJetOverlapRemoval(vector<HEPUtils::Particle*> &lepvec, vector<HEPUtils::Jet*> &jetvec, double DeltaRMax) {
        //Routine to do lepton-jet check
        //Discards leptons if they are within DeltaRMax of a jet

        vector<HEPUtils::Particle*> Survivors;

        for(unsigned int itlep = 0; itlep < lepvec.size(); itlep++) {
          bool overlap = false;
          HEPUtils::P4 lepmom=lepvec.at(itlep)->mom();
          for(unsigned int itjet= 0; itjet < jetvec.size(); itjet++) {
            HEPUtils::P4 jetmom=jetvec.at(itjet)->mom();
            double dR;

            dR=jetmom.deltaR_eta(lepmom);

            if(fabs(dR) <= DeltaRMax) overlap=true;
          }
          if(overlap) continue;
          Survivors.push_back(lepvec.at(itlep));
        }
        lepvec=Survivors;

        return;
      }


    public:

      Analysis_ATLAS_13TeV_RJ3L_lowmass_36invfb() {

	_numSR=0;
	
        NCUTS=10;
        set_luminosity(36.);

        for(int i=0;i<NCUTS;i++){
          cutFlowVector.push_back(0);
          cutFlowVector_str.push_back("");
        }

	// Recursive jigsaw stuff

	
	LAB_3L     = new RestFrames::LabRecoFrame("LAB_3L","lab");
	C1N2_3L    = new RestFrames::DecayRecoFrame("C1N2_3L","#tilde{#chi}^{ #pm}_{1} #tilde{#chi}^{ 0}_{2}");
	C1a_3L     = new RestFrames::DecayRecoFrame("C1a_3L","#tilde{#chi}^{ #pm}_{1}");
	N2b_3L     = new RestFrames::DecayRecoFrame("N2b_3L","#tilde{#chi}^{ 0}_{2}");
	
	L1a_3L      = new RestFrames::VisibleRecoFrame("L1a_3L","#it{l}_{1a}");
	L1b_3L      = new RestFrames::VisibleRecoFrame("L1b_3L","#it{l}_{1b}");
	L2b_3L      = new RestFrames::VisibleRecoFrame("L2b_3L","#it{l}_{2b}");
	
	X1a_3L      = new RestFrames::InvisibleRecoFrame("X1a_3L","#tilde{#chi}^{ 0}_{1 a} + #nu_{a}");
	X1b_3L      = new RestFrames::InvisibleRecoFrame("X1b_3L","#tilde{#chi}^{ 0}_{1 b}");
	
	
	LAB_3L->SetChildFrame(*C1N2_3L);
	
	C1N2_3L->AddChildFrame(*C1a_3L);
	C1N2_3L->AddChildFrame(*N2b_3L);
	
	C1a_3L->AddChildFrame(*L1a_3L);
	C1a_3L->AddChildFrame(*X1a_3L);
	
	N2b_3L->AddChildFrame(*L1b_3L);
	N2b_3L->AddChildFrame(*L2b_3L);
	N2b_3L->AddChildFrame(*X1b_3L);
	
	
	if(LAB_3L->InitializeTree())
	  std::cout << "...Contructor::3L Successfully initialized reconstruction trees" << std::endl;
	else
	  std::cout << "...Constructor::3L Failed initializing reconstruction trees" << std::endl;
	
	//setting the invisible components
	INV_3L = new RestFrames::InvisibleGroup("INV_3L","Invisible system LSP mass Jigsaw");
	INV_3L->AddFrame(*X1a_3L); 
	INV_3L->AddFrame(*X1b_3L);
	
	
	// Set di-LSP mass to minimum Lorentz-invariant expression
	X1_mass_3L = new RestFrames::SetMassInvJigsaw("X1_mass_3L", "Set M_{#tilde{#chi}_{1}^{ 0} #tilde{#chi}_{1}^{ 0}} to minimum");
	INV_3L->AddJigsaw(*X1_mass_3L);
	
	// Set di-LSP rapidity to that of visible particles and neutrino
	X1_eta_3L = new RestFrames::SetRapidityInvJigsaw("X1_eta_3L", "#eta_{#tilde{#chi}_{1}^{ 0} #tilde{#chi}_{1}^{ 0}} = #eta_{3#it{l}}");
	INV_3L->AddJigsaw(*X1_eta_3L);
	X1_eta_3L->AddVisibleFrames(C1N2_3L->GetListVisibleFrames());
	
	
	X1X1_contra_3L = new RestFrames::ContraBoostInvJigsaw("X1X1_contra_3L","Contraboost invariant Jigsaw");
	INV_3L->AddJigsaw(*X1X1_contra_3L);
	X1X1_contra_3L->AddVisibleFrames(C1a_3L->GetListVisibleFrames(),0);
	X1X1_contra_3L->AddVisibleFrames(N2b_3L->GetListVisibleFrames(),1);
	X1X1_contra_3L->AddInvisibleFrames(C1a_3L->GetListInvisibleFrames(),0);
	X1X1_contra_3L->AddInvisibleFrames(N2b_3L->GetListInvisibleFrames(),1);
	
	LAB_3L->InitializeAnalysis();
    	
      }



      void analyze(const HEPUtils::Event* event) {
        HEPUtilsAnalysis::analyze(event);

        // Missing energy
        HEPUtils::P4 ptot = event->missingmom();
    	TVector3 ETMiss;
	ETMiss.SetXYZ(ptot.px(),ptot.py(),0.0);

        // Baseline lepton objects
        vector<HEPUtils::Particle*> baselineElectrons, baselineMuons, baselineTaus;
	
        for (HEPUtils::Particle* electron : event->electrons()) {
          if (electron->pT() > 10. && electron->abseta() < 2.47) baselineElectrons.push_back(electron);
        }
        for (HEPUtils::Particle* muon : event->muons()) {
          if (muon->pT() > 10. && muon->abseta() < 2.4) baselineMuons.push_back(muon);
        }
	
	// Photons
	vector<HEPUtils::Particle*> signalPhotons;
	for (HEPUtils::Particle* photon : event->photons()) {
	  signalPhotons.push_back(photon);
        }
	
	
	// No taus used in 13 TeV analysis?
	//for (HEPUtils::Particle* tau : event->taus()) {
	//if (tau->pT() > 10. && tau->abseta() < 2.47) baselineTaus.push_back(tau);
        //}
        //ATLAS::applyTauEfficiencyR1(baselineTaus);


        // Jets
        vector<HEPUtils::Jet*> bJets;
        vector<HEPUtils::Jet*> nonBJets;
        vector<HEPUtils::Jet*> trueBJets; //for debugging

        // Get b jets
        /// @note We assume that b jets have previously been 100% tagged
        const std::vector<double>  a = {0,10.};
        const std::vector<double>  b = {0,10000.};
        const std::vector<double> c = {0.77}; // set b-tag efficiency to 77%
        HEPUtils::BinnedFn2D<double> _eff2d(a,b,c);
        for (HEPUtils::Jet* jet : event->jets()) {
          bool hasTag=has_tag(_eff2d, jet->eta(), jet->pT());
          if (jet->pT() > 20. && fabs(jet->eta()) < 4.5) {
            if(jet->btag() && hasTag && fabs(jet->eta()) < 2.4 && jet->pT() > 20.){
              bJets.push_back(jet);
            } else {
              nonBJets.push_back(jet);
            }
          }
        }


        // Overlap removal
        vector<HEPUtils::Particle*> signalElectrons;
        vector<HEPUtils::Particle*> signalMuons;
	vector<HEPUtils::Particle*> signalLeptons;
        vector<HEPUtils::Particle*> electronsForVeto;
        vector<HEPUtils::Particle*> muonsForVeto;

        vector<HEPUtils::Jet*> signalJets;
        vector<HEPUtils::Jet*> signalBJets;
        vector<HEPUtils::Jet*> signalNonBJets;

	// Overlap removal is the same as the 8 TeV analysis
        JetLeptonOverlapRemoval(nonBJets,baselineElectrons,0.2);
        LeptonJetOverlapRemoval(baselineElectrons,nonBJets,0.4);
        LeptonJetOverlapRemoval(baselineElectrons,bJets,0.4);
        LeptonJetOverlapRemoval(baselineMuons,nonBJets,0.4);
        LeptonJetOverlapRemoval(baselineMuons,bJets,0.4);

	
	// Also we have already sorted jets by their b tag properties, so reset the b tag variable for each jet to the right category
	// i.e. this was previously 100% true for true b jets then the efficiency map was applied above
        for (HEPUtils::Jet* jet : bJets) {
	  jet->set_btag(true);
	  signalJets.push_back(jet);
	  signalBJets.push_back(jet);
        }
	
        for (HEPUtils::Jet* jet : nonBJets) {
	  if(jet->pT() > 20. && fabs(jet->eta()) < 2.4) {
	    jet->set_btag(false);
	    signalJets.push_back(jet);
	    signalNonBJets.push_back(jet);
	  }
	}

        //Put signal jets in pT order
        std::sort(signalJets.begin(), signalJets.end(), sortByPT_RJ3L);
        std::sort(signalBJets.begin(), signalBJets.end(), sortByPT_RJ3L);
        std::sort(signalNonBJets.begin(), signalNonBJets.end(), sortByPT_RJ3L);

	for (HEPUtils::Particle* electron : baselineElectrons) {
          signalElectrons.push_back(electron);
	  signalLeptons.push_back(electron);
        }

        for (HEPUtils::Particle* muon : baselineMuons) {
          signalMuons.push_back(muon);
	  signalLeptons.push_back(muon);
        }
	
        // We now have the signal electrons, muons, jets and b jets- move on to the analysis
	
	bool m_is3Lep=false;
	if(signalLeptons.size()==3)m_is3Lep=true;
	
	// Define the various variables that will be used later
	
	bool m_pass3L_presel = false;
	bool m_foundSFOS = false;

	double m_lept1Pt  = -999;
	double m_lept1sign=-999;
	
	double m_lept2Pt =-999;
	double m_lept2sign =-999;
	
	double m_lept3Pt =-999;
	double m_lept3sign =-999;

	//Di-Lepton System: Calculated for OS Pairs
	double m_mll=-999;

	//Tri-Lepton System:
	double m_mTW=-999;

	
	// Some lab frame angles and stuff

	double m_H4PP = -999;
	double m_HT4PP = -999;
	double m_RPT_HT4PP = -999;

	if(m_is3Lep){
	  
	  TLorentzVector metLV;
	  //TLorentzVector bigFatJet;
	  metLV.SetPxPyPzE(ptot.px(),ptot.py(),0.,sqrt(ptot.px()*ptot.px()+ptot.py()*ptot.py()));
	  
	  //Put the Leptons in a more useful form
	  vector<pair<TLorentzVector,int> > myLeptons;
	  //vector<lep> myLeptons;
	  for(unsigned int ilep=0; ilep<signalLeptons.size(); ilep++)
	    {
	      pair<TLorentzVector,int> temp;
	      TLorentzVector tlv_temp;
	      
	      tlv_temp.SetPtEtaPhiM(signalLeptons[ilep]->pT(),signalLeptons[ilep]->eta(),signalLeptons[ilep]->phi(),0.0);
	      temp.first = tlv_temp;
	      int lepton_charge=0;
	      if(signalLeptons[ilep]->pid()<0)lepton_charge=-1;
	      if(signalLeptons[ilep]->pid()>0)lepton_charge=1;
	      temp.second = lepton_charge;
	      //temp.third = lepton_origin->at(lep_signal_index[ilep]);
	      //temp.fourth = lepton_type->at(lep_signal_index[ilep]);
	      //temp = make_tuple(tlv_temp,lepton_charge->at(lep_signal_index[ilep]),lepton_origin->at(lep_signal_index_[ilep]),lepton_type->at(lepton_signal_index[ilep]));
	      myLeptons.push_back(temp);
	    }

	  sort(myLeptons.begin(), myLeptons.end(), SortLeptons);

	  
	  
	  if(myLeptons[0].first.Pt()<25.0 || myLeptons[1].first.Pt()<25.0 || myLeptons[2].first.Pt()<20.0)
	    {
	      m_pass3L_presel=false;
	    }
	  else{
	    m_pass3L_presel=true;
	  }

	  //if(!m_pass3L_presel)return;
	  
	  //Tri-Lepton System
	  //Here we choose leptons based on where they "come from"
	  //lept1 and lept2 are the lepton pair with invariant mass closest to the Z-Mass
	  //lept3 is the remaining lepton
	  //This is meant to emulate lept1 and lept2 being produced by the Z, while lept3 is produced by the W
	  
	  double diff = 10000000000.0;
	  int Zlep1 = -99;
	  int Zlep2 = -99;
	  double Zmass = -999.0;
	  bool foundSFOS = false;
	  
	  for(unsigned int i=0; i<myLeptons.size(); i++)
	    {
	      for(unsigned int j=i+1; j<myLeptons.size(); j++)
		{
		  //Opposite-Sign
                if(myLeptons[i].second*myLeptons[j].second<0)
		  {
                    //Same-Flavor
                    if(abs(myLeptons[i].second)==abs(myLeptons[j].second))
		      {
                        double mass = (myLeptons[i].first+myLeptons[j].first).M();
			double massdiff = fabs(mass-91.1876);
                        if(massdiff<diff)
			  {
                            diff=massdiff;
                            Zmass=mass;
                            Zlep1 = i;
                            Zlep2 = j;
                            foundSFOS = true;
			  }
		      }
		  }
		}
	    }

	  if(!foundSFOS)
	    {
	      m_foundSFOS=false;
	    }
	  else {
	    m_foundSFOS=true;
	  }
	  
	  if(m_foundSFOS){

	    int Wlep1 = -999;
	    if( (Zlep1==0 && Zlep2==1) || (Zlep1==1 && Zlep2==0) ) Wlep1=2;
	    else if( (Zlep1==0 && Zlep2==2) || (Zlep1==2 && Zlep2==0) ) Wlep1=1;
	    else if((Zlep1==1 && Zlep2==2) || (Zlep1==2 && Zlep2==1) ) Wlep1=0;

	    //Knowing the indices, we perform assignments
	    m_lept1Pt   = myLeptons[0].first.Pt();
	    m_lept1sign = myLeptons[0].second;
	    
	    m_lept2Pt   = myLeptons[1].first.Pt();
	    m_lept2sign = myLeptons[1].second;
	    
	    m_lept3Pt   = myLeptons[2].first.Pt();
	    m_lept3sign = myLeptons[2].second;
	    
	    
	    m_mll  = Zmass; //based on mass minimization

	    vector<TLorentzVector> Leptons;
	    Leptons.push_back(myLeptons[Wlep1].first);
	    Leptons.push_back(myLeptons[Zlep1].first);
	    Leptons.push_back(myLeptons[Zlep2].first);
	    
	    
	    TLorentzVector metLV;
	    metLV.SetPxPyPzE(ptot.px(),ptot.py(),0.,sqrt(ptot.px()*ptot.px()+ptot.py()*ptot.py()));
	    
	    double wlepMetphi = myLeptons[Wlep1].first.DeltaPhi(metLV);
	    
	    m_mTW = sqrt(2*myLeptons[Wlep1].first.Pt()*metLV.Pt()*(1-cos(wlepMetphi)));  
	    	    
	    INV_3L->SetLabFrameThreeVector(ETMiss); //set the MET in the event
	    
	    
	    L1a_3L->SetLabFrameFourVector(Leptons[0]); // Set lepton from W
	    L1b_3L->SetLabFrameFourVector(Leptons[1]); // Set lepton1 from Z
	    L2b_3L->SetLabFrameFourVector(Leptons[2]); // Set lepton2 from Z

	    if(!LAB_3L->AnalyzeEvent()) cout << "FillEvent:: Something went wrong..." << endl;
	    
	    TLorentzVector l1;
	    TLorentzVector l2;
	    TLorentzVector l3 = L1a_3L->GetFourVector(*LAB_3L);
	  
	    //if(DEBUG) cout << "WlepPt: " << m_WlepPt << " Wlepsign: " << m_Wlepsign << endl;
	    if (L1b_3L->GetFourVector(*LAB_3L).Pt() > L2b_3L->GetFourVector(*LAB_3L).Pt()){
	      l1 = L1b_3L->GetFourVector(*LAB_3L);
	      l2 = L2b_3L->GetFourVector(*LAB_3L);
	      
	    }
	    else {
	      
	   
	      l2 = L1b_3L->GetFourVector(*LAB_3L);
	      l1 = L2b_3L->GetFourVector(*LAB_3L);
	    }


	    // More lab frame stuff
 
	    //if(DEBUG)  cout << "Zlep1: " << m_Zlep1Pt << " " << m_Zlep1sign << " Zlep2Pt: " << m_Zlep2Pt << " " << m_Zlep2sign << endl;
	    TLorentzVector vP_V1aPP  = L1a_3L->GetFourVector(*C1N2_3L);
	    TLorentzVector vP_V1bPP  = L1b_3L->GetFourVector(*C1N2_3L);
	    TLorentzVector vP_V2bPP  = L2b_3L->GetFourVector(*C1N2_3L);
	    TLorentzVector vP_I1aPP  = X1a_3L->GetFourVector(*C1N2_3L);
	    TLorentzVector vP_I1bPP  = X1b_3L->GetFourVector(*C1N2_3L);
	    
	    TLorentzVector vP_V1aPa  = L1a_3L->GetFourVector(*C1a_3L);
	    TLorentzVector vP_I1aPa  = X1a_3L->GetFourVector(*C1a_3L);
	    
	    TLorentzVector vP_V1bPb = L1b_3L->GetFourVector(*N2b_3L);
	    TLorentzVector vP_V2bPb = L2b_3L->GetFourVector(*N2b_3L);
	    TLorentzVector vP_I1bPb = X1b_3L->GetFourVector(*N2b_3L);
	    
	    

	    //Variables w/ 4 objects 
	    
	    /// Defined in the PP-frame
	    //Four vector sum of all visible objets + four vector sum of inv objects
	  	    
	    //Scalar sum of all visible objects + vector sum of invisible momenta 
	    m_H4PP = vP_V1aPP.P() + vP_V1bPP.P() + vP_V2bPP.P() + (vP_I1aPP + vP_I1bPP).P();//H(3,1)PP
	    m_HT4PP = vP_V1aPP.Pt() + vP_V1bPP.Pt() + vP_V2bPP.Pt() + (vP_I1aPP + vP_I1bPP).Pt();//HT(3,1)PP
	    
	    // Invisible components again
	    TLorentzVector vP_IaLAB = X1a_3L->GetFourVector(*LAB_3L);
	    TLorentzVector vP_IbLAB = X1b_3L->GetFourVector(*LAB_3L);
	  
	    
	    // Testing for low mass 3L
	    TLorentzVector p_Ia_Lab = X1a_3L->GetFourVector(*LAB_3L);
	    TLorentzVector p_Ib_Lab = X1b_3L->GetFourVector(*LAB_3L);
	    TVector3 lab_to_pp = C1N2_3L->GetBoostInParentFrame();
	    
	    /// Defined in the P-frame    
	    
	    ////Calculation of dRll_I_PP;
	    //m_dRll_I_PP = (vP_V1bPP+vP_V1bPP).DeltaR(vP_I1bPP);
	    //m_R_Ib_Ia = (vP_V1bPP + vP_V2bPP + vP_I1bPP).P()/(vP_V1aPP+vP_I1aPP).P();
	    
	    // signal variables
	    TLorentzVector vP_Va = C1a_3L->GetVisibleFourVector(*C1a_3L);
	    TLorentzVector vP_Vb = N2b_3L->GetVisibleFourVector(*N2b_3L);
	    
	    TVector3 vP_PP = C1N2_3L->GetFourVector(*LAB_3L).Vect();
	    double Pt_PP = vP_PP.Pt();
	    m_RPT_HT4PP = Pt_PP / (Pt_PP + m_HT4PP);
	    
	    	    
	    // mt_min here
	    
	    /*double min0 = -999;
	    double min1 = -999;
	    double min2 = -999;
	    double lepmetphi0 = myLeptons[0].first.DeltaPhi(metLV);
	    double lepmetphi1 = myLeptons[1].first.DeltaPhi(metLV);
	    double lepmetphi2 = myLeptons[2].first.DeltaPhi(metLV);

	    if (myLeptons[0].second == -myLeptons[1].second) min0 =  sqrt(2*myLeptons[2].first.Pt()*metLV.Pt()*(1-cos(lepmetphi2)));
	    if (myLeptons[0].second == -myLeptons[2].second) min1 =  sqrt(2*myLeptons[1].first.Pt()*metLV.Pt()*(1-cos(lepmetphi1)));
	    if (myLeptons[1].second == -myLeptons[2].second) min2 =  sqrt(2*myLeptons[0].first.Pt()*metLV.Pt()*(1-cos(lepmetphi0)));
	    
	    if (min0 > 0 && min1 > 0) m_min_mt = min(min0,min1);
	    else if (min0 > 0 && min2 > 0) m_min_mt = min(min0,min2);
	    else if (min1 > 0 && min2 > 0) m_min_mt = min(min1,min2);
	    else if (min0 > 0 && min1 < 0 && min2 < 0) m_min_mt = min0;
	    else if (min1 > 0 && min0 < 0 && min2 < 0) m_min_mt = min1;
	    else if (min2 > 0 && min0 < 0 && min1 <0) m_min_mt = min2;*/
	  	    
	    
	  } // end of if(m_foundSFOS)
	  
	} // end of if(m_is3Lep)


	// Cutflow check

	cutFlowVector_str[0] = "No cuts ";
        cutFlowVector_str[1] = "Preselection ";
        cutFlowVector_str[2] = "75 GeV < mll < 105 GeV ";
        cutFlowVector_str[3] = "mTW > 100 GeV ";
        cutFlowVector_str[4] = "m_HT4PP/m_H4PP > 0.9 ";
        cutFlowVector_str[5] = "m_H4PP > 250 GeV ";
        cutFlowVector_str[6] = "pT_PP/(pT_PP + HT_PP(3,1)) ";
	
	for(int j=0;j<NCUTS;j++){
	  
	  if( (j==0) ||
	      
	      (j==1 && m_is3Lep && (((m_lept1sign*m_lept2sign<0 && abs(m_lept1sign)==abs(m_lept2sign)) || (m_lept1sign*m_lept3sign<0 && abs(m_lept1sign)==abs(m_lept3sign)) || (m_lept2sign*m_lept3sign<0 && abs(m_lept2sign)==abs(m_lept3sign)))) && m_lept1Pt>60. && m_lept2Pt>40. && m_lept3Pt>30. && signalBJets.size()==0 && signalJets.size()==0) ||
	      
	      (j==2 && m_is3Lep && (((m_lept1sign*m_lept2sign<0 && abs(m_lept1sign)==abs(m_lept2sign)) || (m_lept1sign*m_lept3sign<0 && abs(m_lept1sign)==abs(m_lept3sign)) || (m_lept2sign*m_lept3sign<0 && abs(m_lept2sign)==abs(m_lept3sign)))) && m_lept1Pt>60. && m_lept2Pt>40. && m_lept3Pt>30. && signalBJets.size()==0 && signalJets.size()==0 && m_mll>75. && m_mll<105.) ||
	      
	      (j==3 && m_is3Lep && (((m_lept1sign*m_lept2sign<0 && abs(m_lept1sign)==abs(m_lept2sign)) || (m_lept1sign*m_lept3sign<0 && abs(m_lept1sign)==abs(m_lept3sign)) || (m_lept2sign*m_lept3sign<0 && abs(m_lept2sign)==abs(m_lept3sign)))) && m_lept1Pt>60. && m_lept2Pt>40. && m_lept3Pt>30. && signalBJets.size()==0 && signalJets.size()==0 && m_mll>75. && m_mll<105. && m_mTW>100.) ||
	      
	      (j==4 && m_is3Lep && (((m_lept1sign*m_lept2sign<0 && abs(m_lept1sign)==abs(m_lept2sign)) || (m_lept1sign*m_lept3sign<0 && abs(m_lept1sign)==abs(m_lept3sign)) || (m_lept2sign*m_lept3sign<0 && abs(m_lept2sign)==abs(m_lept3sign)))) && m_lept1Pt>60. && m_lept2Pt>40. && m_lept3Pt>30. && signalBJets.size()==0 && signalJets.size()==0 && m_mll>75. && m_mll<105. && m_mTW>100. && m_HT4PP/m_H4PP > 0.9) ||
	      
	      (j==5 && m_is3Lep && (((m_lept1sign*m_lept2sign<0 && abs(m_lept1sign)==abs(m_lept2sign)) || (m_lept1sign*m_lept3sign<0 && abs(m_lept1sign)==abs(m_lept3sign)) || (m_lept2sign*m_lept3sign<0 && abs(m_lept2sign)==abs(m_lept3sign)))) && m_lept1Pt>60. && m_lept2Pt>40. && m_lept3Pt>30. && signalBJets.size()==0 && signalJets.size()==0 && m_mll>75. && m_mll<105. && m_mTW>100. && m_HT4PP/m_H4PP > 0.9 && m_H4PP > 250.) ||
	      
	      (j==6 && m_is3Lep && (((m_lept1sign*m_lept2sign<0 && abs(m_lept1sign)==abs(m_lept2sign)) || (m_lept1sign*m_lept3sign<0 && abs(m_lept1sign)==abs(m_lept3sign)) || (m_lept2sign*m_lept3sign<0 && abs(m_lept2sign)==abs(m_lept3sign)))) && m_lept1Pt>60. && m_lept2Pt>40. && m_lept3Pt>30. && signalBJets.size()==0 && signalJets.size()==0 && m_mll>75. && m_mll<105. && m_mTW>100. && m_HT4PP/m_H4PP > 0.9 && m_H4PP > 250. && m_RPT_HT4PP < 0.05)
	      
	      )cutFlowVector[j]++;
	}
	
	// Now apply the signal region cuts

	if(m_is3Lep && (((m_lept1sign*m_lept2sign<0 && abs(m_lept1sign)==abs(m_lept2sign)) || (m_lept1sign*m_lept3sign<0 && abs(m_lept1sign)==abs(m_lept3sign)) || (m_lept2sign*m_lept3sign<0 && abs(m_lept2sign)==abs(m_lept3sign)))) && m_lept1Pt>60. && m_lept2Pt>40. && m_lept3Pt>30. && signalBJets.size()==0 && signalJets.size()==0 && m_mll>75. && m_mll<105. && m_mTW>100. && m_HT4PP/m_H4PP > 0.9 && m_H4PP > 250. && m_RPT_HT4PP < 0.05)_numSR++;

	
	return;
	
      }
      

      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        HEPUtilsAnalysis::add(other);

        Analysis_ATLAS_13TeV_RJ3L_lowmass_36invfb* specificOther
                = dynamic_cast<Analysis_ATLAS_13TeV_RJ3L_lowmass_36invfb*>(other);

        // Here we will add the subclass member variables:
        if (NCUTS != specificOther->NCUTS) NCUTS = specificOther->NCUTS;
        for (int j=0; j<NCUTS; j++) {
          cutFlowVector[j] += specificOther->cutFlowVector[j];
          cutFlowVector_str[j] = specificOther->cutFlowVector_str[j];
        }

	_numSR+= specificOther->_numSR;

      }


      void collect_results() {

	double scale_by=1.;
	cout << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
	cout << "CUT FLOW: ATLAS 13 TeV 3 lep low mass RJ signal region "<<endl;
	cout << "------------------------------------------------------------------------------------------------------------------------------"<<endl;
	cout<< right << setw(40) << "CUT" << setw(20) << "RAW" << setw(20) << "SCALED"
	    << setw(20) << "%" << setw(20) << "clean adj RAW"<< setw(20) << "clean adj %" << endl;
	for (int j=0; j<NCUTS; j++) {
	  cout << right << setw(40) << cutFlowVector_str[j].c_str() << setw(20)
	       << cutFlowVector[j] << setw(20) << cutFlowVector[j]*scale_by << setw(20)
	       << 100.*cutFlowVector[j]/cutFlowVector[0] << "%" << setw(20)
	       << cutFlowVector[j]*scale_by << setw(20) << 100.*cutFlowVector[j]/cutFlowVector[0]<< "%" << endl;
	}
	cout << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;

	/// Register results objects with the results for each SR; obs & bkg numbers from the paper

	static const string ANAME = "Analysis_ATLAS_13TeV_RJ3L_lowmass_36invfb";

	/*int _numSRA_TT, _numSRA_TW, _numSRA_T0;
	int _numSRB_TT, _numSRB_TW, _numSRB_T0;
	int _numSRC1, _numSRC2, _numSRC3, _numSRC4, _numSRC5;
	int _numSRD_low, _numSRD_high, _numSRE;*/
	
        add_result(SignalRegionData(ANAME, "SR", 20, {_numSR,  0.}, {10.31, 1.96}));     
	
        /*SignalRegionData results_SRA_TT;
        results_SRA_TT.analysis_name = "Analysis_ATLAS_13TeV_RJ3L_lowmass_36invfb";
        results_SRA_TT.sr_label = "SRA_TT";
        results_SRA_TT.n_observed = 11.;
        results_SRA_TT.n_background = 15.8;
        results_SRA_TT.background_sys = 1.9;
        results_SRA_TT.signal_sys = 0.;
        results_SRA_TT.n_signal = _numSRA1;

        add_result(results_SRA_TT);*/

	 //deletion of RJ 3-lepton pointers
	delete LAB_3L;
	delete C1N2_3L;
	delete C1a_3L;
	delete N2b_3L;
	delete L1a_3L;
	delete L1b_3L;
	delete L2b_3L;
	delete X1a_3L;
	delete X1b_3L;
	delete INV_3L;
	delete X1_mass_3L;
	delete X1_eta_3L;
	delete X1X1_contra_3L;
	
	
        return;
      }

    };


    DEFINE_ANALYSIS_FACTORY(ATLAS_13TeV_RJ3L_lowmass_36invfb)


  }
}
