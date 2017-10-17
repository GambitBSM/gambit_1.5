#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>

#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "gambit/ColliderBit/mt2_bisect.h"

using namespace std;

/* The ATLAS 0 lepton direct stop analysis

   Based on: https://arxiv.org/abs/1709.04183
  
   Code by Martin White

   KNOWN ISSUES

   a) Should apply Very Loose selection to electron candidates
   b) Cannot apply requirement that the ETmiss calculated from tracking information is aligned in phi with that calculated from the calo system.
   c) We do not apply the tau veto. Could approximate by removing events with tagged taus in?

*/

namespace Gambit {
  namespace ColliderBit {

    bool sortByPT13(HEPUtils::Jet* jet1, HEPUtils::Jet* jet2) { return (jet1->pT() > jet2->pT()); }
    bool sortByMass(HEPUtils::Jet* jet1, HEPUtils::Jet* jet2) { return (jet1->mass() > jet2->mass()); }

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
    
    
    class Analysis_ATLAS_13TeV_0LEPStop_20invfb : public HEPUtilsAnalysis {
    private:

      // Numbers passing cuts
      int _numSRA_TT, _numSRA_TW, _numSRA_T0;
      int _numSRB_TT, _numSRB_TW, _numSRB_T0;
      int _numSRC1, _numSRC2, _numSRC3, _numSRC4, _numSRC5;
      int _numSRD_low, _numSRD_high, _numSRE;

      vector<int> cutFlowVector;
      vector<string> cutFlowVector_str;
      int NCUTS; //=16;

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

      Analysis_ATLAS_13TeV_0LEPStop_20invfb() {

	_numSRA_TT=0;  _numSRA_TW=0; _numSRA_T0=0;
	_numSRB_TT=0; _numSRB_TW=0; _numSRB_T0=0;
	_numSRC1=0;  _numSRC2=0;  _numSRC3=0; _numSRC4=0; _numSRC5=0;
	_numSRD_low=0; _numSRD_high=0; _numSRE=0;
	
        NCUTS=120;
        set_luminosity(36.);

        for(int i=0;i<NCUTS;i++){
          cutFlowVector.push_back(0);
          cutFlowVector_str.push_back("");
        }

      }



      void analyze(const HEPUtils::Event* event) {
        HEPUtilsAnalysis::analyze(event);

        // Missing energy
        HEPUtils::P4 ptot = event->missingmom();
        double met = event->met();


        // Baseline lepton objects
        vector<HEPUtils::Particle*> baselineElectrons, baselineMuons, baselineTaus;
	
        for (HEPUtils::Particle* electron : event->electrons()) {
          if (electron->pT() > 7. && electron->abseta() < 2.47) baselineElectrons.push_back(electron);
        }
        for (HEPUtils::Particle* muon : event->muons()) {
          if (muon->pT() > 6. && muon->abseta() < 2.7) baselineMuons.push_back(muon);
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
          if (jet->pT() > 20. && fabs(jet->eta()) < 2.8) {
            if(jet->btag() && hasTag && fabs(jet->eta()) < 2.5 && jet->pT() > 20.){
              bJets.push_back(jet);
            } else {
              nonBJets.push_back(jet);
            }
          }
        }


        // Overlap removal
        vector<HEPUtils::Particle*> signalElectrons;
        vector<HEPUtils::Particle*> signalMuons;
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

	
	// It seems that there are no extra signal jet requirements (unlike 8 TeV analysis)
	// Also we have already sorted jets by their b tag properties, so reset the b tag variable for each jet to the right category
	// i.e. this was previously 100% true for true b jets then the efficiency map was applied above
        for (HEPUtils::Jet* jet : bJets) {
	  jet->set_btag(true);
	  signalJets.push_back(jet);
	  signalBJets.push_back(jet);
        }
	
        for (HEPUtils::Jet* jet : nonBJets) {
	  jet->set_btag(false);
	  signalJets.push_back(jet);
	  signalNonBJets.push_back(jet);
	}

        //Put signal jets in pT order
        std::sort(signalJets.begin(), signalJets.end(), sortByPT13);
        std::sort(signalBJets.begin(), signalBJets.end(), sortByPT13);
        std::sort(signalNonBJets.begin(), signalNonBJets.end(), sortByPT13);

	for (HEPUtils::Particle* electron : baselineElectrons) {
          signalElectrons.push_back(electron);
        }

        for (HEPUtils::Particle* muon : baselineMuons) {
          signalMuons.push_back(muon);
        }

        // We now have the signal electrons, muons, jets and b jets- move on to the analysis

	// Common preselection for all signal regions

	// At least four jets, one of which is b-tagged
	bool passPreSelJetCut=false;
	int nJets = signalJets.size();
	if(nJets>=4){
          if(signalJets[0]->pT() > 80.
             && signalJets[1]->pT() > 80.
             && signalJets[2]->pT() > 40.
             && signalJets[3]->pT() > 40.){
	    if(signalJets[0]->btag() ||
	       signalJets[1]->btag() ||
	       signalJets[2]->btag() ||
	       signalJets[3]->btag())passPreSelJetCut=true;
	  }
        }
	
        // Lepton veto
        int nElectrons = signalElectrons.size();
        int nMuons = signalMuons.size();
        bool cut_LeptonVeto=true;
        if((nElectrons + nMuons)>0.)cut_LeptonVeto=false;

	//MET > 250 GeV
        bool cut_METGt250=false;
        if(met>250.)cut_METGt250=true;

	//MET > 400 GeV
        bool cut_METGt400=false;
        if(met>400.)cut_METGt400=true;

	//MET > 500 GeV
	bool cut_METGt500=false;
        if(met>500.)cut_METGt500=true;

	//MET > 550 GeV
	bool cut_METGt550=false;
        if(met>550.)cut_METGt550=true;

        //Calculate dphi(jet,met) for the two leading jets
        bool cut_dPhiJetsPresel=false;
        bool cut_dPhiJet2=false;
        bool cut_dPhiJet1=false;
        double dphi_jetmet1=9999;
        if(nJets>0)dphi_jetmet1=std::acos(std::cos(signalJets.at(0)->phi()-ptot.phi()));
        double dphi_jetmet2=9999;
        if(nJets>1)dphi_jetmet2=std::acos(std::cos(signalJets.at(1)->phi()-ptot.phi()));
	if(dphi_jetmet2>0.4)cut_dPhiJet2=true;
        if(dphi_jetmet1>0.4)cut_dPhiJet1=true;
        if(cut_dPhiJet1 && cut_dPhiJet2)cut_dPhiJetsPresel=true;

	bool passPresel=false;
	if(passPreSelJetCut && cut_LeptonVeto && cut_METGt250 && cut_dPhiJetsPresel)passPresel=true;

	
	/*********************************************************/
	/*                                                       */
	/* SIGNAL REGIONS A & B                                  */
	/*                                                       */
	/*********************************************************/
	
        //Number of b jets >=2
        bool passBJetCut_AB=false;
        if(signalBJets.size()>=2)passBJetCut_AB=true;

	// Have an additional requirement on dPhiJetMet for 3rd jet
	bool cut_dPhiJet3=false;
	bool cut_dPhiJets_AB=false;
	double dphi_jetmet3=9999;
	if(nJets>2)dphi_jetmet3=std::acos(std::cos(signalJets.at(2)->phi()-ptot.phi()));
	if(dphi_jetmet3>0.4)cut_dPhiJet3=true;
	if(cut_dPhiJetsPresel && cut_dPhiJet3)cut_dPhiJets_AB=true;

	// Need to recluster jets at this point (R=0.8 and R=1.2)
	vector<HEPUtils::Jet*> signalJets_0_8=get_jets(signalJets,0.8);
	vector<HEPUtils::Jet*> signalJets_1_2=get_jets(signalJets,1.2);

	// Require two R=1.2 jets
	int numJets_1_2=signalJets_1_2.size();
	//bool passJet_1_2_Cut=false;
	//if(numJets_1_2==2)passJet_1_2_Cut=true;

	// Classify events based on R=1.2 jet content
	//Put signal jets in decreasing mass order
        std::sort(signalJets_1_2.begin(), signalJets_1_2.end(), sortByMass);
	/*bool isTT=false;
	bool isTW=false;
	bool isT0=false;
	if(passJet_1_2_Cut){
	  if(signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()>120.)isTT=true;
	  if(signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()>60. && signalJets_1_2[1]->mass()<120.)isTW=true;
	  if(signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<60.)isT0=true;
	  }*/

	//Calculate dphi(b,met) for the closest b-jet in phi to MET
	// Put phi in range -pi to pi
	
        double dphi_bjetmet_min=9999.;
        double minphi =9999.;
        int whichb=0;
        for(size_t j=0; j<signalBJets.size(); j++) {
          if(fabs(std::acos(std::cos(signalBJets.at(j)->phi()-ptot.phi())))<minphi) {
            minphi = fabs(std::acos(std::cos(signalBJets.at(j)->phi()-ptot.phi())));
            dphi_bjetmet_min = minphi;
            whichb=j;
          }
        }

	double mT_bjetmet_min = 0;
        if(passBJetCut_AB) mT_bjetmet_min = sqrt(2*signalBJets.at(whichb)->pT()*met*(1-std::cos(dphi_bjetmet_min)));
	bool cut_mT_bmin_Gt200=false;
        if(mT_bjetmet_min>200.)cut_mT_bmin_Gt200=true;

	bool cut_mT_bmin_Gt250=false;
        if(mT_bjetmet_min>250.)cut_mT_bmin_Gt250=true;
	
	bool cut_mT_bmin_Gt350=false;
        if(mT_bjetmet_min>350.)cut_mT_bmin_Gt350=true;


	
	//Calculate dphi(b,met) for the furthest b-jet in phi to MET
        double dphi_bjetmet_max=0.;
        double maxphi = 0.;
        whichb=0;
        for(size_t j=0; j<signalBJets.size(); j++) {
          if(fabs(std::acos(std::cos(signalBJets.at(j)->phi()-ptot.phi())))>maxphi) {
            maxphi = fabs(std::acos(std::cos(signalBJets.at(j)->phi()-ptot.phi())));
            dphi_bjetmet_max = maxphi;
            whichb=j;
          }
        }

	double mT_bjetmet_max = 0;
        if(passBJetCut_AB) mT_bjetmet_max = sqrt(2*signalBJets.at(whichb)->pT()*met*(1-std::cos(dphi_bjetmet_max)));
	bool cut_mT_bmax_Gt200=false;
        if(mT_bjetmet_max>200.)cut_mT_bmax_Gt200=true;

	bool cut_mT_bmax_Gt300=false;
        if(mT_bjetmet_max>300.)cut_mT_bmax_Gt300=true;
	
	bool cut_mT_bmax_Gt450=false;
        if(mT_bjetmet_max>450.)cut_mT_bmax_Gt450=true;
	
	//Sort 0.8 jets by pT
        std::sort(signalJets_0_8.begin(), signalJets_0_8.end(), sortByPT13);
	//bool cut_leading08jetMGt60=false;
	//if(signalJets_0_8.size()>0 && signalJets_0_8[0]->mass()>60.)cut_leading08jetMGt60=true;

	// Apply requirement on stransverse mass

	// First need to find jet combinations that minimise the chi-squared
	// Paper claims that "Initially, single or pairs of R=0.4 jets form W candidates
	// which are then combined with additional b-tagged jets in the event to construct top candidates"
	// So we need to pick 2 b-jets (at random?) if there are more than 2 in the event.
	// We need to try and make top candidates using the combinations (b+nonb,b+nonb), (b+nonb+nonb,b+nonb) and (b+nonb+nonb,b+nonb+nonb)

	// First make two R=0.4 jet collections:
	// topreco_bjets: contains 2 b jets picked at random from the event (t=top)
	// topreco_ojets: contains all other jets in the event (potentially including some b jets, o=other)
	double topMass = 173.2;
	double wMass = 80.4;
      
	vector<int> b_jet_indices;
	for(unsigned int i=0;i<signalBJets.size();i++){
	  b_jet_indices.push_back(i);
	}

	double mt2=0;
	
	// Note that the whole calculation assumes that there are two b jets
	// The ATLAS cutflow applies this cut last, which I don't understand
	// Perhaps they take the jets with the highest b weight for this calculation (and the other cuts)
	// The SRs will not be affected, but the cutflow will be
	// We will instead choose the two b jets at random
	
	int bjet0=-1;
	int bjet1=-1;
	
	if(signalBJets.size()>1){
	  
	  random_unique(b_jet_indices.begin(), b_jet_indices.end(), 2);
	  bjet0=b_jet_indices[0];
	  bjet1=b_jet_indices[1];
	  
	  vector<HEPUtils::Jet*> topreco_bjets;
	  vector<HEPUtils::Jet*> topreco_ojets;
	  for (unsigned int i=0;i<signalBJets.size();i++) {
	    if(i==bjet0 || i ==bjet1){
	    topreco_bjets.push_back(signalBJets[i]);
	    }
	    else {
	      topreco_ojets.push_back(signalBJets[i]);
	    }
	  }
	  
	  for (unsigned int i=0;i<signalNonBJets.size();i++) {
	    topreco_ojets.push_back(signalNonBJets[i]);
	  }
	  
	  // Now loop over and find combinations
	  	  
	  double Chi2min = 999999999999999.;
	  double Chi2=0;
	  int W1j1_low = -1,W1j2_low = -1,W2j1_low = -1,W2j2_low = -1,b1_low = -1,b2_low = -1;
	  
	  HEPUtils::P4 topCand0, topCand1, WCand0, WCand1;
	  if(passPresel){
	    
	    for(int W1j1=0; W1j1<(int)topreco_ojets.size(); W1j1++) {
	      for(int W2j1=0;W2j1<(int)topreco_ojets.size(); W2j1++) {
		if (W2j1==W1j1) continue;
		for(int b1=0;b1<(int)topreco_bjets.size();b1++){
		  for(int b2=0;b2<(int)topreco_bjets.size();b2++){
		    if(b2==b1) continue;

		    double chi21, chi22, mW1, mW2, mt1, mt2;
		    //try bl,bl top candidates
		    if(W2j1>W1j1){
		      mW1 = topreco_ojets.at(W1j1)->mass();
		      mW2 = topreco_ojets.at(W2j1)->mass();
		      mt1 = (topreco_ojets.at(W1j1)->mom() + topreco_bjets.at(b1)->mom()).m();
		      mt2 = (topreco_ojets.at(W2j1)->mom() + topreco_bjets.at(b2)->mom()).m();
		      chi21 = (mW1-wMass)*(mW1-wMass)/wMass + (mt1-topMass)*(mt1-topMass)/topMass;
		      chi22 = (mW2-wMass)*(mW2-wMass)/wMass + (mt2-topMass)*(mt2-topMass)/topMass;

		      if(Chi2min > (chi21 + chi22)){
			Chi2min = chi21 + chi22;
			if(chi21 < chi22){
			  W1j1_low = W1j1;
			  W1j2_low = -1;
			  W2j1_low = W2j1;
			  W2j2_low = -1;
			  b1_low = b1;
			  b2_low = b2;
			}
			else{
			  W2j1_low = W1j1;
			  W2j2_low = -1;
			  W1j1_low = W2j1;
			  W1j2_low = -1;
			  b2_low = b1;
			  b1_low = b2;
			}
		      }
		    }

		    if(topreco_ojets.size() < 3) continue;
		    for(int W1j2=W1j1+1;W1j2<(int)topreco_ojets.size(); W1j2++) {
		      if(W1j2==W2j1) continue;

		      //try bll,bl top candidates
		      mW1 = (topreco_ojets.at(W1j1)->mom() + topreco_ojets.at(W1j2)->mom()).m();
		      mW2 = topreco_ojets.at(W2j1)->mass();
		      mt1 = (topreco_ojets.at(W1j1)->mom() + topreco_ojets.at(W1j2)->mom() + topreco_bjets.at(b1)->mom()).m();
		      mt2 = (topreco_ojets.at(W2j1)->mom() + topreco_bjets.at(b2)->mom()).m();
		      chi21 = (mW1-wMass)*(mW1-wMass)/wMass + (mt1-topMass)*(mt1-topMass)/topMass;
		      chi22 = (mW2-wMass)*(mW2-wMass)/wMass + (mt2-topMass)*(mt2-topMass)/topMass;

		      if(Chi2min > (chi21 + chi22)){
			Chi2min = chi21 + chi22;
			if(chi21 < chi22){
			  W1j1_low = W1j1;
			  W1j2_low = W1j2;
			  W2j1_low = W2j1;
			  W2j2_low = -1;
			  b1_low = b1;
			  b2_low = b2;
			}
			else{
			  W2j1_low = W1j1;
			  W2j2_low = W1j2;
			  W1j1_low = W2j1;
			  W1j2_low = -1;
			  b2_low = b1;
			  b1_low = b2;
			}
		      }

		      if(topreco_ojets.size() < 4)continue;
		      //try bll, bll top candidates
		      for(int W2j2=W2j1+1;W2j2<(int)topreco_ojets.size(); W2j2++){
			if((W2j2==W1j1) || (W2j2==W1j2)) continue;
			if(W2j1<W1j1) continue;  //runtime reasons, we don't want combinations checked twice <--------------------This line should be added
			mW1 = (topreco_ojets.at(W1j1)->mom() + topreco_ojets.at(W1j2)->mom()).m();
			mW2 = (topreco_ojets.at(W2j1)->mom() + topreco_ojets.at(W2j2)->mom()).m();
			mt1 = (topreco_ojets.at(W1j1)->mom() + topreco_ojets.at(W1j2)->mom() + topreco_bjets.at(b1)->mom()).m();
			mt2 = (topreco_ojets.at(W2j1)->mom() + topreco_ojets.at(W2j2)->mom() + topreco_bjets.at(b2)->mom()).m();
			chi21 = (mW1-wMass)*(mW1-wMass)/wMass + (mt1-topMass)*(mt1-topMass)/topMass;
			chi22 = (mW2-wMass)*(mW2-wMass)/wMass + (mt2-topMass)*(mt2-topMass)/topMass;

			if(Chi2min > (chi21 + chi22)){
			  Chi2min = chi21 + chi22;
			  if(chi21 < chi22){
			    W1j1_low = W1j1;
			    W1j2_low = W1j2;
			    W2j1_low = W2j1;
			    W2j2_low = W2j2;
			    b1_low = b1;
			    b2_low = b2;
			  }
			  else{
			    W2j1_low = W1j1;
			    W2j2_low = W1j2;
			    W1j1_low = W2j1;
			    W1j2_low = W2j2;
			    b2_low = b1;
			    b1_low = b2;
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }

	    if(W1j2_low == -1) WCand0 = topreco_ojets.at(W1j1_low)->mom();
	    else WCand0 = topreco_ojets.at(W1j1_low)->mom() + topreco_ojets.at(W1j2_low)->mom();
	    topCand0 = WCand0 + topreco_bjets.at(b1_low)->mom();
	    if(W2j2_low == -1) WCand1 = topreco_ojets.at(W2j1_low)->mom();
	    else WCand1 = topreco_ojets.at(W2j1_low)->mom() + topreco_ojets.at(W2j2_low)->mom();
	    topCand1 = WCand1 + topreco_bjets.at(b2_low)->mom();
	    //Chi2 = Chi2min;
	    
	    // We now have our top candidate four vectors, so calculate MT2
	    // Use only magnitude of pT and direction in phi, then set top mass by hand
	    // Need to replace unless we can confirm that we use this
	    
	    HEPUtils::P4 tempTop0=HEPUtils::P4::mkEtaPhiMPt(0.,topCand0.phi(),173.210,topCand0.pT());
	    HEPUtils::P4 tempTop1=HEPUtils::P4::mkEtaPhiMPt(0.,topCand1.phi(),173.210,topCand1.pT());

	    // Note that the first component here is the mass
	    // This must be the top mass (i.e. mass of our vectors) and not zero!
	    
	    double pa_a[3] = { tempTop0.m() , tempTop0.px(), tempTop0.py() };
	    double pb_a[3] = { tempTop1.m() , tempTop1.px(), tempTop1.py() };
	    double pmiss_a[3] = { 0, ptot.px(), ptot.py() };
	    double mn_a = 0.;
	    
	    mt2_bisect::mt2 mt2_event_a;
	    
	    mt2_event_a.set_momenta(pa_a,pb_a,pmiss_a);
	    mt2_event_a.set_mn(mn_a);
	    
	    mt2 = mt2_event_a.get_mt2();
	    
	  }
	}
	bool cutMT2Gt400 = false;
	if(mt2 > 400.)cutMT2Gt400=true;

	bool cutMT2Gt500 = false;
	if(mt2 > 500.)cutMT2Gt500=true;
	
	// Now implement deltaR b,b
	// Use the two randomly selected b jets that we chose earlier

	double dRBB=0;
	
	if(signalBJets.size()>1){
	  
	  HEPUtils::P4 bjet0mom=signalBJets.at(bjet0)->mom();
	  HEPUtils::P4 bjet1mom=signalBJets.at(bjet1)->mom();
	  dRBB=bjet0mom.deltaR_eta(bjet1mom);
        }

	bool cutdRBBGt1=false;
	if(dRBB>1)cutdRBBGt1=true;

	bool cutdRBBGt12=false;
	if(dRBB>1.2)cutdRBBGt12=true;

	bool cutdRBBGt08=false;
	if(dRBB>0.8)cutdRBBGt08=true;
	
	// Calculate HT
	// Scalar sum of pT of all R=0.4 jets in the event

	double HT=0;
	
	for(int jet=0;jet<signalJets.size();jet++)HT+=signalJets[jet]->pT();
	
	// Calculate met_sig
	// met/sqrt(HT)
	
	double metsig=met/sqrt(HT);

	bool devSkim = false;
	
	if( (HT > 150.) || (signalElectrons.size() > 0 && signalElectrons[0]->pT() > 100.) || (signalElectrons.size() > 1 && signalElectrons[0]->pT() > 20. && signalElectrons[1]->pT() > 20.) || (signalMuons.size() > 0 && signalMuons[0]->pT() > 100.) || (signalMuons.size() > 1 && signalMuons[0]->pT() > 20. && signalMuons[1]->pT() > 20.) || (signalPhotons.size() > 0 && signalPhotons[0]->pT() > 100.) || (signalPhotons.size() > 1 && signalPhotons[0]->pT() > 50. && signalPhotons[1]->pT() > 50.))devSkim=true;
	
	bool isSRA_TT=false;
	bool isSRA_TW=false;
	bool isSRA_T0=false;
	bool isSRB_TT=false;
	bool isSRB_TW=false;
	bool isSRB_T0=false;
	bool isSRC1=false;
	bool isSRC2=false;
	bool isSRC3=false;
	bool isSRC4=false;
	bool isSRC5=false;
	bool isSRD_low=false;
	bool isSRD_high=false;
	bool isSRE=false;
	
        /*if(cut_LeptonVeto && cut_Btag && cut_METGt150 && cut_dPhiJets && cut_mTbjetmetGt175 && cut_6jets && mbjj0.m() < 225. && mbjj1.m() < 250. && cut_tau && cut_METGt150)isSRA1=true;

        if(cut_LeptonVeto && cut_Btag && cut_METGt150 && cut_dPhiJets && cut_mTbjetmetGt175 && cut_6jets && mbjj0.m() < 225. && mbjj1.m() < 250. && cut_tau && cut_METGt250)isSRA2=true;

        if(cut_LeptonVeto && cut_Btag && cut_METGt150 && cut_dPhiJets && cut_mTbjetmetGt175 && cut_6jets && mbjj0.m() > 50. && mbjj0.m() < 250. && mbjj1.m() > 50. && mbjj1.m() < 400. && mtMin > 50. && cut_tau && cut_METGt300)isSRA3=true;

        if(cut_LeptonVeto && cut_Btag && cut_METGt150 && cut_dPhiJets && cut_mTbjetmetGt175 && cut_6jets && mbjj0.m() > 50. && mbjj0.m() < 250. && mbjj1.m() > 50. && mbjj1.m() < 400. && mtMin > 50. && cut_tau && cut_METGt350)isSRA4=true;*/

	cutFlowVector_str[0] = "No cuts ";
	cutFlowVector_str[1] = "Derivation skim";
        cutFlowVector_str[2] = "Lepton veto ";
        cutFlowVector_str[3] = "Njets >= 4 ";
        cutFlowVector_str[4] = "Nbjets >= 1 ";
        cutFlowVector_str[5] = "met > 250 GeV ";
        cutFlowVector_str[6] = "dPhi(jet,MET) > 0.4 ";
        cutFlowVector_str[7] = "pT jet 1 > 80 GeV ";
        cutFlowVector_str[8] = "pT jet 3 > 40 GeV ";
	cutFlowVector_str[9] = "m jet0, R=1.2 > 120 GeV ";
	cutFlowVector_str[10] = "SRA-TT: m jet1, R=1.2 > 120 GeV";
	cutFlowVector_str[11] = "SRA-TT: met > 400 GeV";	
	cutFlowVector_str[12] = "SRA-TT: m jet0, R=0.8 > 60 GeV ";
        cutFlowVector_str[13] = "SRA-TT: mT(b,MET) min > 200 ";
        cutFlowVector_str[14] = "SRA-TT: deltaR(b,b) > 1 ";
        cutFlowVector_str[15] = "SRA-TT: mT2 > 400 GeV";
	cutFlowVector_str[16] = "SRA-TT: Nbjets >=2 ";
        cutFlowVector_str[17] = "SRA-TW: m jet1, R=1.2 < 120 GeV";
	cutFlowVector_str[18] = "SRA-TW: m jet1, R=1.2 > 60 GeV";
	cutFlowVector_str[19] = "SRA-TW: met > 500 GeV ";
	cutFlowVector_str[20] = "SRA-TW: m jet0, R=0.8 > 60 GeV";
	cutFlowVector_str[21] = "SRA-TW: mT(b,MET) min > 200 GeV";
	cutFlowVector_str[22] = "SRA-TW: mT2 > 400 GeV ";
	cutFlowVector_str[23] = "SRA-TW: Nbjets >=2 ";
	cutFlowVector_str[24] = "SRA-T0: m jet1, R=1.2 < 60 GeV";
	cutFlowVector_str[25] = "SRA-T0: m jet0, R=0.8 > 60 GeV";
	cutFlowVector_str[26] = "SRA-T0: met > 550 GeV ";
	cutFlowVector_str[27] = "SRA-T0: mT(b,MET) min > 200 GeV";
	cutFlowVector_str[28] = "SRA-T0: mT2 > 500 GeV ";
	cutFlowVector_str[29] = "SRA-T0: Nbjets >=2 ";
	cutFlowVector_str[30] = "SRB-TT: m jet1, R=1.2 > 120 GeV";
	cutFlowVector_str[31] = "SRB-TT: deltaR(b,b) > 1.2";
	cutFlowVector_str[32] = "SRB-TT: mT(b,MET) max > 200 GeV";
	cutFlowVector_str[33] = "SRB-TT: mT(b,MET) min > 200 GeV";
	cutFlowVector_str[34] = "SRB-TT: Nbjets >=2 ";
	cutFlowVector_str[35] = "SRB-TW: m jet1, R=1.2 < 120 GeV";
	cutFlowVector_str[36] = "SRB-TW: m jet1, R=1.2 > 60 GeV";
	cutFlowVector_str[37] = "SRB-TW: deltaR(b,b) > 1.2";
	cutFlowVector_str[38] = "SRB-TW: mT(b,MET) max > 200 GeV";
	cutFlowVector_str[39] = "SRB-TW: mT(b,MET) min > 200 GeV";
	cutFlowVector_str[40] = "SRB-TW: Nbjets >=2 ";
	cutFlowVector_str[41] = "SRB-T0: m jet1, R=1.2 < 60 GeV";
	cutFlowVector_str[42] = "SRB-T0: mT(b,MET) min > 200 GeV";
	cutFlowVector_str[43] = "SRB-T0: deltaR(b,b) > 1.2";
	cutFlowVector_str[44] = "SRB-T0: mT(b,MET) max > 200 GeV";
	cutFlowVector_str[45] = "SRB-T0: met > 250 GeV ";
	cutFlowVector_str[46] = "SRB-T0: Nbjets >=2 ";

	// Cutflow for SRD
	cutFlowVector_str[47] = "SRD-high: No cuts ";
	cutFlowVector_str[48] = "SRD-high: Derivation skim";
        cutFlowVector_str[49] = "SRD-high: Lepton veto ";
        cutFlowVector_str[50] = "SRD-high: Njets >= 4 ";
        cutFlowVector_str[51] = "SRD-high: Nbjets >= 1 ";
        cutFlowVector_str[52] = "SRD-high: met > 250 GeV ";
        cutFlowVector_str[53] = "SRD-high: dPhi(jet,MET) > 0.4 ";
        cutFlowVector_str[54] = "SRD-high: pT jet 1 > 80 GeV ";
        cutFlowVector_str[55] = "SRD-high: pT jet 3 > 40 GeV ";
	cutFlowVector_str[56] = "SRD-high: Njets >= 5 ";
	cutFlowVector_str[57] = "SRD-high: pT jet 1 > 150 ";
	cutFlowVector_str[58] = "SRD-high: pT jet 3 > 80 ";
	cutFlowVector_str[59] = "SRD-high: pT jet 4 > 60 ";
	cutFlowVector_str[60] = "SRD-high: mT(b,MET) min > 350 GeV ";
	cutFlowVector_str[61] = "SRD-high: mT(b,MET) max > 450 GeV ";
	cutFlowVector_str[62] = "SRD-high: Nbjets >=2 ";
	cutFlowVector_str[63] = "SRD-high: met > 250 GeV ";
	cutFlowVector_str[64] = "SRD-high: deltaR(b,b) > 0.8";
	cutFlowVector_str[65] = "SRD-high: pT0b + pT1b > 400 GeV";
	cutFlowVector_str[66] = "SRD-low: Njets >=5";
	cutFlowVector_str[67] = "SRD-low: NBjets >=2";
	cutFlowVector_str[68] = "SRD-low: met > 250 GeV";
	cutFlowVector_str[69] = "SRD-low: mT(b,MET) min > 250 GeV ";
	cutFlowVector_str[70] = "SRD-low: mT(b,MET) max > 300 GeV ";
	cutFlowVector_str[71] = "SRD-low: deltaR(b,b) > 0.8";
	cutFlowVector_str[72] = "SRD-low: pT jet 1 > 150 GeV ";
	cutFlowVector_str[73] = "SRD-low: pT jet 3 > 100 GeV ";
	cutFlowVector_str[74] = "SRD-low: pT jet 4 > 60 GeV ";
	cutFlowVector_str[75] = "SRD-low: pT0b + pT1b > 300 GeV";
	
	// Cutflow for SRE
	cutFlowVector_str[76] = "SRE: met > 550 GeV";
	cutFlowVector_str[77] = "SRE: m jet0, R = 0.8 > 120 GeV";
	cutFlowVector_str[78] = "SRE: m jet1, R = 0.8 > 80 GeV";
	cutFlowVector_str[79] = "SRE: HT > 800 GeV";
	cutFlowVector_str[80] = "SRE: met/sqrt(HT) > 18 GeV^1/2";
	cutFlowVector_str[81] = "SRE: mT(b,MET) min > 200 GeV";
	cutFlowVector_str[82] = "SRE: NBjets >=2";

	
        for(int j=0;j<NCUTS;j++){
          if(
             (j==0) ||

	     (j==1 && devSkim) ||
	     
             (j==2 && devSkim && cut_LeptonVeto) ||
	     
             (j==3 && devSkim && cut_LeptonVeto && signalJets.size()>3) ||
	     
             (j==4 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0) ||
	     
             (j==5 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250) ||
	     
             (j==6 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB) ||

             (j==7 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80.) ||

	     (j==8 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. )  ||

	     (j==9 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120.)  ||

	     // SRA-TT
	     
	     (j==10 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()>120.)  ||

	     (j==11 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt400 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()>120.)  ||

	     (j==12 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt400 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()>120. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60.) ||

	     (j==13 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt400 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()>120. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_mT_bmin_Gt200) ||

	     (j==14 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt400 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()>120. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_mT_bmin_Gt200 && cutdRBBGt1) ||

	     (j==15 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt400 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()>120. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_mT_bmin_Gt200 && cutdRBBGt1 && cutMT2Gt400) ||

	     (j==16 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt400 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()>120. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_mT_bmin_Gt200 && cutdRBBGt1 && cutMT2Gt400) ||

	     // SRA-TW

	     /*cutFlowVector_str[16] = "SRA-TW: m jet1, R=1.2 < 120 GeV";
	     cutFlowVector_str[17] = "SRA-TW: m jet1, R=1.2 > 60 GeV";
	     cutFlowVector_str[18] = "SRA-TW: met > 500 GeV ";
	     cutFlowVector_str[19] = "SRA-TW: m jet0, R=0.8 > 60 GeV";
	     cutFlowVector_str[20] = "SRA-TW: mT(b,MET) min > 200 ";
	     cutFlowVector_str[21] = "SRA-TW: mT2 > 400 GeV ";
	     cutFlowVector_str[22] = "SRA-TW: Nbjets >=2 ";*/

	     (j==17 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<120.)  ||

	     (j==18 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60.)  ||

	     (j==19 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60. && cut_METGt500)  ||

	     (j==20 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60. && cut_METGt500 && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60.)  ||

	     (j==21 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60. && cut_METGt500 && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_mT_bmin_Gt200)  ||

	     (j==22 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60. && cut_METGt500 && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_mT_bmin_Gt200 && cutMT2Gt400)  ||

	     (j==23 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60. && cut_METGt500 && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_mT_bmin_Gt200 && cutMT2Gt400) ||

	     /*cutFlowVector_str[23] = "SRA-T0: m jet1, R=1.2 < 60 GeV";
	     cutFlowVector_str[24] = "SRA-T0: m jet0, R=0.8 > 60 GeV";
	     cutFlowVector_str[25] = "SRA-T0: met > 550 GeV ";
	     cutFlowVector_str[26] = "SRA-T0: mT(b,MET) min > 200 ";
	     cutFlowVector_str[27] = "SRA-T0: mT2 > 500 GeV ";
	     cutFlowVector_str[28] = "SRA-T0: Nbjets >=2 ";*/

	     (j==24 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<60.)  ||

	     (j==25 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<60. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60.)  ||

	     (j==26 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<60. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_METGt550)  ||

	     (j==27 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<60. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_METGt550 &&  cut_mT_bmin_Gt200)  ||

	     (j==28 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<60. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_METGt550 &&  cut_mT_bmin_Gt200 && cutMT2Gt500)  ||

	     (j==29 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<60. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_METGt550 &&  cut_mT_bmin_Gt200 && cutMT2Gt500) ||
	     
	     /* cutFlowVector_str[29] = "SRB-TT: m jet1, R=1.2 > 120 GeV";
	     cutFlowVector_str[30] = "SRB-TT: deltaR(b,b) > 1.2";
	     cutFlowVector_str[31] = "SRB-TT: mT(b,MET) max > 200 GeV";
	     cutFlowVector_str[32] = "SRB-TT: mT(b,MET) min > 200 GeV";
	     cutFlowVector_str[33] = "SRB-TT: Nbjets >=2 ";*/

	     (j==30 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()>120.)  ||

	     (j==31 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()>120. && cutdRBBGt12)  ||

	     (j==32 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()>120. && cutdRBBGt12 && cut_mT_bmax_Gt200)  ||

	     (j==33 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()>120. && cutdRBBGt12 && cut_mT_bmax_Gt200 && cut_mT_bmin_Gt200)  ||

	     (j==34 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()>120. && cutdRBBGt12 && cut_mT_bmax_Gt200 && cut_mT_bmin_Gt200)  ||
	     
	     /*cutFlowVector_str[34] = "SRB-TW: m jet1, R=1.2 < 120 GeV";
	     cutFlowVector_str[35] = "SRB-TW: m jet1, R=1.2 > 60 GeV";
	     cutFlowVector_str[36] = "SRB-TW: deltaR(b,b) > 1.2";
	     cutFlowVector_str[37] = "SRB-TW: mT(b,MET) max > 200 GeV";
	     cutFlowVector_str[38] = "SRB-TW: mT(b,MET) min > 200 GeV";
	     cutFlowVector_str[39] = "SRB-TW: Nbjets >=2 ";*/

	     (j==35 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<120.)  ||

	     (j==36 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60.)  ||

	     (j==37 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60. && cutdRBBGt12)  ||

	     (j==38 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60. && cutdRBBGt12 &&  cut_mT_bmax_Gt200)   ||

	     (j==39 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60. && cutdRBBGt12 &&  cut_mT_bmax_Gt200 && cut_mT_bmin_Gt200)   ||

	     (j==40 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60. && cutdRBBGt12 &&  cut_mT_bmax_Gt200 && cut_mT_bmin_Gt200)   ||


	     /* cutFlowVector_str[40] = "SRB-T0: m jet1, R=1.2 < 60 GeV";
		cutFlowVector_str[41] = "SRB-T0: mT(b,MET) min > 200 GeV";
		cutFlowVector_str[42] = "SRB-T0: deltaR(b,b) > 1.2";
		cutFlowVector_str[43] = "SRB-T0: mT(b,MET) max > 200 GeV";
		cutFlowVector_str[44] = "SRB-T0: met > 250 GeV ";
		cutFlowVector_str[45] = "SRB-T0: Nbjets >=2 "; */

	     (j==41 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<60.)  ||

	     (j==42 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<60. && cut_mT_bmin_Gt200)  ||

	     (j==43 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<60. && cut_mT_bmin_Gt200 && cutdRBBGt12)  ||

	     (j==44 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<60. && cut_mT_bmin_Gt200 && cutdRBBGt12 && cut_mT_bmax_Gt200)  ||

	     (j==45 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<60. && cut_mT_bmin_Gt200 && cutdRBBGt12 && cut_mT_bmax_Gt200)  ||

	     (j==46 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<60. && cut_mT_bmin_Gt200 && cutdRBBGt12 && cut_mT_bmax_Gt200) ||


	     /* cutFlowVector_str[47] = "SRD-high: No cuts ";
		cutFlowVector_str[48] = "SRD-high: Derivation skim";
		cutFlowVector_str[49] = "SRD-high: Lepton veto ";
		cutFlowVector_str[50] = "SRD-high: Njets >= 4 ";
		cutFlowVector_str[51] = "SRD-high: Nbjets >= 1 ";
		cutFlowVector_str[52] = "SRD-high: met > 250 GeV ";
		cutFlowVector_str[53] = "SRD-high: dPhi(jet,MET) > 0.4 ";
		cutFlowVector_str[54] = "SRD-high: pT jet 1 > 80 GeV ";
		cutFlowVector_str[55] = "SRD-high: pT jet 3 > 40 GeV ";
		cutFlowVector_str[56] = "SRD-high: Njets >= 5 ";
		cutFlowVector_str[57] = "SRD-high: pT jet 1 > 150 ";
		cutFlowVector_str[58] = "SRD-high: pT jet 3 > 80 ";
		cutFlowVector_str[59] = "SRD-high: pT jet 4 > 60 ";
		cutFlowVector_str[60] = "SRD-high: mT(b,MET) min > 350 GeV ";
		cutFlowVector_str[61] = "SRD-high: mT(b,MET) max > 450 GeV ";
		cutFlowVector_str[62] = "SRD-high: Nbjets >=2 ";
		cutFlowVector_str[63] = "SRD-high: met > 250 GeV ";
		cutFlowVector_str[64] = "SRD-high: deltaR(b,b) > 0.8";
		cutFlowVector_str[65] = "SRD-high: pT0b + pT1b > 400 GeV"; */

	     (j==47) ||
	     
	     (j==48 && devSkim) ||
	     
             (j==49 && devSkim && cut_LeptonVeto) ||
	     
             (j==50 && devSkim && cut_LeptonVeto && signalJets.size()>3) ||
	     
             (j==51 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0) ||
	     
             (j==52 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250) ||
	     
             (j==53 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB) ||
	     
             (j==54 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80.) ||
	     
	     (j==55 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. )  ||

	     (j==56 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. )  ||

	     (j==57 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>40. )  ||

	     (j==58 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>80. )  ||
	     
	     (j==59 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>80. && signalJets[4]->pT()>60.)  ||

	     (j==60 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>80. && signalJets[4]->pT()>60. && cut_mT_bmin_Gt350)  ||
	     
	     (j==61 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>80. && signalJets[4]->pT()>60. && cut_mT_bmin_Gt350 && cut_mT_bmax_Gt450)  ||

	     (j==62 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>80. && signalJets[4]->pT()>60. && cut_mT_bmin_Gt350 && cut_mT_bmax_Gt450)  ||

	     (j==63 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>80. && signalJets[4]->pT()>60. && cut_mT_bmin_Gt350 && cut_mT_bmax_Gt450)  ||

	     (j==64 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>80. && signalJets[4]->pT()>60. && cut_mT_bmin_Gt350 && cut_mT_bmax_Gt450 && cutdRBBGt08)  ||

	     (j==65 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>80. && signalJets[4]->pT()>60. && cut_mT_bmin_Gt350 && cut_mT_bmax_Gt450 && cutdRBBGt08 && ( (signalBJets[bjet0]->pT() + signalBJets[bjet1]->pT())>400.)) || 
	     
	     /*cutFlowVector_str[66] = "SRD-low: Njets >=5";
	     cutFlowVector_str[67] = "SRD-low: NBjets >=2";
	     cutFlowVector_str[68] = "SRD-low: met > 250 GeV";
	     cutFlowVector_str[69] = "SRD-low: mT(b,MET) min > 250 GeV ";
	     cutFlowVector_str[70] = "SRD-low: mT(b,MET) max > 300 GeV ";
	     cutFlowVector_str[71] = "SRD-low: deltaR(b,b) > 0.8";
	     cutFlowVector_str[72] = "SRD-low: pT jet 1 > 150 GeV ";
	     cutFlowVector_str[73] = "SRD-low: pT jet 3 > 100 GeV ";
	     cutFlowVector_str[74] = "SRD-low: pT jet 4 > 60 GeV ";
	     cutFlowVector_str[75] = "SRD-low: pT0b + pT1b > 300 GeV";*/
	     
	     
	     (j==66 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>0 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40.)  ||

	     (j==67 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40.)  ||

	     (j==68 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40.)  ||

	     (j==69 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && cut_mT_bmin_Gt250)  ||

	     (j==70 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && cut_mT_bmin_Gt250 && cut_mT_bmax_Gt300)  ||

	     (j==71 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && cut_mT_bmin_Gt250 && cut_mT_bmax_Gt300 && cutdRBBGt08)  ||

	     (j==72 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>40. && cut_mT_bmin_Gt250 && cut_mT_bmax_Gt300 && cutdRBBGt08)  ||

	     (j==73 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>100. && cut_mT_bmin_Gt250 && cut_mT_bmax_Gt300 && cutdRBBGt08)  ||

	     (j==74 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>100. && signalJets[4]->pT()>60. && cut_mT_bmin_Gt250 && cut_mT_bmax_Gt300 && cutdRBBGt08)  ||
	     
	     (j==75 && devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>100. && signalJets[4]->pT()>60. && cut_mT_bmin_Gt250 && cut_mT_bmax_Gt300 && cutdRBBGt08 &&  ( (signalBJets[bjet0]->pT() + signalBJets[bjet1]->pT())>300.))  ||

	     /* cutFlowVector_str[76] = "SRE: met > 550 GeV";
		cutFlowVector_str[77] = "SRE: m jet0, R = 0.8 > 120 GeV";
		cutFlowVector_str[78] = "SRE: m jet1, R = 0.8 > 80 GeV";
		cutFlowVector_str[79] = "SRE: HT > 800 GeV";
		cutFlowVector_str[80] = "SRE: met/sqrt(HT) > 18 GeV^1/2";
		cutFlowVector_str[81] = "SRE: mT(b,MET) min > 200 GeV";
		cutFlowVector_str[82] = "SRE: NBjets >=2";*/

	     (j==76 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt550 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40.)  ||

	     (j==77 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt550 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>120.)  ||

	     (j==78 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt550 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_0_8.size() > 1 && signalJets_0_8[0]->mass()>120. && signalJets_0_8[0]->mass()>80.)  ||

	     (j==79 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt550 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_0_8.size() > 1 && signalJets_0_8[0]->mass()>120. && signalJets_0_8[0]->mass()>80. && HT > 800.)  ||

	     (j==80 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt550 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_0_8.size() > 1 && signalJets_0_8[0]->mass()>120. && signalJets_0_8[0]->mass()>80. && HT > 800. && metsig > 18.)  ||

	     (j==81 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>0 && cut_METGt550 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_0_8.size() > 1 && signalJets_0_8[0]->mass()>120. && signalJets_0_8[0]->mass()>80. && HT > 800. && metsig > 18. && cut_mT_bmin_Gt200)  ||

	     (j==82 && devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt550 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_0_8.size() > 1 && signalJets_0_8[0]->mass()>120. && signalJets_0_8[0]->mass()>80. && HT > 800. && metsig > 18. && cut_mT_bmin_Gt200)
	     
	     // Still to do...
	     
	     
	     ){
	    
            cutFlowVector[j]++;
	    
          }
	}
	  

        //We're now ready to apply the cuts for each signal region
        //_numSR1, _numSR2, _numSR3;

        /*if(isSRA1)_numSRA1++;
        if(isSRA2)_numSRA2++;
        if(isSRA3)_numSRA3++;
        if(isSRA4)_numSRA4++;

        if(isSRC1)_numSRC1++;
        if(isSRC2)_numSRC2++;
        if(isSRC3)_numSRC3++;*/

	if(devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt400 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()>120. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_mT_bmin_Gt200 && cutdRBBGt1 && cutMT2Gt400)isSRA_TT=true;

	if(devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60. && cut_METGt500 && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_mT_bmin_Gt200 && cutMT2Gt400)isSRA_TW=true;

	if(devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && signalJets_1_2[1]->mass()<60. && signalJets_0_8.size() > 0 && signalJets_0_8[0]->mass()>60. && cut_METGt550 &&  cut_mT_bmin_Gt200 && cutMT2Gt500)isSRA_T0=true;

	if(devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[0]->mass()>120. && cutdRBBGt12 && cut_mT_bmax_Gt200 && cut_mT_bmin_Gt200)isSRB_TT=true;

	if(devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<120. && signalJets_1_2[1]->mass()>60. && cutdRBBGt12 &&  cut_mT_bmax_Gt200 && cut_mT_bmin_Gt200)isSRB_TW=true;

	if(devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_1_2.size()>0 && signalJets_1_2[0]->mass()>120. && signalJets_1_2.size()>1 && signalJets_1_2[1]->mass()<60. && cut_mT_bmin_Gt200 && cutdRBBGt12 && cut_mT_bmax_Gt200)isSRB_T0=true;

	if(devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>80. && signalJets[4]->pT()>60. && cut_mT_bmin_Gt350 && cut_mT_bmax_Gt450 && cutdRBBGt08 && ( (signalBJets[bjet0]->pT() + signalBJets[bjet1]->pT())>400.))isSRD_high=true;

	if(devSkim && cut_LeptonVeto && signalJets.size()>4 && signalBJets.size()>1 && cut_METGt250 && cut_dPhiJets_AB && signalJets[1]->pT()>150. && signalJets[3]->pT()>100. && signalJets[4]->pT()>60. && cut_mT_bmin_Gt250 && cut_mT_bmax_Gt300 && cutdRBBGt08 &&  ( (signalBJets[bjet0]->pT() + signalBJets[bjet1]->pT())>300.))isSRD_low=true;  

	if(devSkim && cut_LeptonVeto && signalJets.size()>3 && signalBJets.size()>1 && cut_METGt550 && cut_dPhiJets_AB && signalJets[1]->pT()>80. && signalJets[3]->pT()>40. && signalJets_0_8.size() > 1 && signalJets_0_8[0]->mass()>120. && signalJets_0_8[0]->mass()>80. && HT > 800. && metsig > 18. && cut_mT_bmin_Gt200)isSRE=true;
	

	if(isSRA_TT)_numSRA_TT++;
	if(isSRA_TW)_numSRA_TW++;
	if(isSRA_T0)_numSRA_T0++;
	if(isSRB_TT)_numSRB_TT++;
	if(isSRB_TW)_numSRB_TW++;
	if(isSRB_T0)_numSRB_T0++;
	if(isSRC1)_numSRC1++;
	if(isSRC2)_numSRC2++;
	if(isSRC3)_numSRC3++;
	if(isSRC4)_numSRC4++;
	if(isSRC5)_numSRC5++;
	if(isSRD_low)_numSRD_low++;
	if(isSRD_high)_numSRD_high++;
	if(isSRE)_numSRE++;
	
        return;

      }


      void add(BaseAnalysis* other) {
        // The base class add function handles the signal region vector and total # events.
        HEPUtilsAnalysis::add(other);

        Analysis_ATLAS_13TeV_0LEPStop_20invfb* specificOther
                = dynamic_cast<Analysis_ATLAS_13TeV_0LEPStop_20invfb*>(other);

        // Here we will add the subclass member variables:
        if (NCUTS != specificOther->NCUTS) NCUTS = specificOther->NCUTS;
        for (int j=0; j<NCUTS; j++) {
          cutFlowVector[j] += specificOther->cutFlowVector[j];
          cutFlowVector_str[j] = specificOther->cutFlowVector_str[j];
        }

	_numSRA_TT += specificOther->_numSRA_TT;
	_numSRA_TW += specificOther->_numSRA_TW;
	_numSRA_T0 += specificOther->_numSRA_T0;
	_numSRB_TT += specificOther->_numSRB_TT;
	_numSRB_TW += specificOther->_numSRB_TW;
	_numSRB_T0 += specificOther->_numSRB_T0;
	_numSRC1 += specificOther->_numSRC1;
	_numSRC2 += specificOther->_numSRC2;
	_numSRC3 += specificOther->_numSRC3;
	_numSRC4 += specificOther->_numSRC4;
	_numSRC5 += specificOther->_numSRC5;
	_numSRD_low += specificOther->_numSRD_low;
	_numSRD_high += specificOther->_numSRD_high;
	_numSRE += specificOther->_numSRE;
      }


      void collect_results() {

	double scale_by=1.;
	cout << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
	cout << "CUT FLOW: ATLAS 13 TeV 0 lep stop paper "<<endl;
	cout << "------------------------------------------------------------------------------------------------------------------------------"<<endl;
	cout<< right << setw(40) << "CUT" << setw(20) << "RAW" << setw(20) << "SCALED"
	    << setw(20) << "%" << setw(20) << "clean adj RAW"<< setw(20) << "clean adj %" << endl;
	for (unsigned int j=0; j<NCUTS; j++) {
	  cout << right << setw(40) << cutFlowVector_str[j].c_str() << setw(20)
	       << cutFlowVector[j] << setw(20) << cutFlowVector[j]*scale_by << setw(20)
	       << 100.*cutFlowVector[j]/cutFlowVector[0] << "%" << setw(20)
	       << cutFlowVector[j]*scale_by << setw(20) << 100.*cutFlowVector[j]/cutFlowVector[0]<< "%" << endl;
	}
	cout << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;

	/// Register results objects with the results for each SR; obs & bkg numbers from the paper

	static const string ANAME = "Analysis_ATLAS_13TeV_0LEPStop_13invfb";

	/*int _numSRA_TT, _numSRA_TW, _numSRA_T0;
	int _numSRB_TT, _numSRB_TW, _numSRB_T0;
	int _numSRC1, _numSRC2, _numSRC3, _numSRC4, _numSRC5;
	int _numSRD_low, _numSRD_high, _numSRE;*/
	
        add_result(SignalRegionData(ANAME, "SRA-TT", 11, {_numSRA_TT,  0.}, {8.6, 2.1}));
        add_result(SignalRegionData(ANAME, "SRA-TW", 9, {_numSRA_TW,  0.}, {9.3, 2.2}));
        add_result(SignalRegionData(ANAME, "SRA-T0",  18, {_numSRA_T0,  0.}, {18.7, 2.7}));
        add_result(SignalRegionData(ANAME, "SRB-TT",  38, {_numSRB_TT,  0.}, { 39.3,  7.6}));
        add_result(SignalRegionData(ANAME, "SRB-TW", 53, {_numSRB_TW,  0.}, {52.4, 7.4}));
        add_result(SignalRegionData(ANAME, "SRB-T0", 206, {_numSRB_T0,  0.}, { 179.,  26.}));
        //add_result(SignalRegionData(ANAME, "SRC1", 20, {_numSRC1,  0.}, { 20.6,  6.5}));
        //add_result(SignalRegionData(ANAME, "SRC2", 22, {_numSRC2,  0.}, { 27.6,  4.9}));
        //add_result(SignalRegionData(ANAME, "SRC3", 22, {_numSRC3,  0.}, {  18.9, 3.4}));
        //add_result(SignalRegionData(ANAME, "SRC4", 1, {_numSRC4,  0.}, {  7.7, 1.2}));
        //add_result(SignalRegionData(ANAME, "SRC5", 0, {_numSRC5, 0.}, { 0.91,  0.73}));
        add_result(SignalRegionData(ANAME, "SRD-low", 27, {_numSRD_low, 0.}, {  25.1, 6.2}));
        add_result(SignalRegionData(ANAME, "SRD-high", 11, {_numSRD_high, 0.}, {  8.5,1.5}));
	add_result(SignalRegionData(ANAME, "SRE", 3, {_numSRE, 0.}, {  3.64,0.79}));

	
        /*SignalRegionData results_SRA_TT;
        results_SRA_TT.analysis_name = "Analysis_ATLAS_13TeV_0LEPStop_20invfb";
        results_SRA_TT.sr_label = "SRA_TT";
        results_SRA_TT.n_observed = 11.;
        results_SRA_TT.n_background = 15.8;
        results_SRA_TT.background_sys = 1.9;
        results_SRA_TT.signal_sys = 0.;
        results_SRA_TT.n_signal = _numSRA1;

        add_result(results_SRA_TT);*/

        return;
      }

    };


    DEFINE_ANALYSIS_FACTORY(ATLAS_13TeV_0LEPStop_20invfb)


  }
}
