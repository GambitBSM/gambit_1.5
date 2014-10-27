/* ReadDelphesTree.cc was original the code used to create templates from Delphes FastSim.
 * The original code has been modified to just read the truth and reco particles and
 * fill histograms for further comparison and testing
 * A. Saavedra 10-10-2014
 *
 * Below is the original information.
 *
 * Extracts a template in MET-Njets parameter space from a given ROOT file,
 * with a specified cut placed on the jet transverse momentum. Along the Njets
 * axis there are two bins, one containing events with no jets above the cut,
 * and the other containing events with at least one jet above the cut. Along
 * the MET axis there are 20 bins containing events with energies from 0 to 200
 * GeV in increments of 10 GeV, with overflow (MET > 200 GeV) assigned to the 
 * last bin.
 *
 * The parameter space is only populated with the event if the fastsim (reconstructed)
 * leptons are matched with truth_prompt leptons
 *
 * A number of histograms show fastsim and truth distributions
 * for the different physics objects. 
 * From DecayChannel the acceptance and efficiency can be calculated
 *
 * 
 * Takes as input five arguments, in the following order:
 *  - Input file containing truth record of the process of interest (string)
 *      [N.B. This must be in the Delphes ROOT ntuple output format]
 *  - Name of the file to which the template is written (string)
 *  - Desired normalisation of the histogram (float)
 *  - Name of the template that will be created (string)
 *
 * Created:   Jan   2014  C. Loi
 * Modified:  Feb   2014  C. Suster
 * Modified:  June  2014  A. Saavedra
 *
 */

#include "ReadDelphesTree.h"

#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include <set>

#include <sstream>

#define ELECTRON_PT_MIN 25.0
#define ELECTRON_ABSETA_MAX 2.47
#define ELECTRON_ABSETA_CRACKMIN 1.37
#define ELECTRON_ABSETA_CRACKMAX 1.52
#define MUON_PT_MIN 20.0
#define MUON_ABSETA_MAX 2.5
#define JET_PT_MIN 25.0
#define JET_MULTIPLICITY_PT_MIN 30.0
#define JET_ABSETA_MAX 2.5

#define DR_ELECTRON_JET_REMOVE_JET 0.2
#define DR_ELECTRON_JET_REMOVE_ELECTRON 0.4
#define DR_MUON_JET_REMOVE_JET 0.2
#define DR_MUON_JET_REMOVE_MUON 0.4

#define DELPHES_FLAG_ELECTRON_EFF_MEDIUM 0
#define DELPHES_FLAG_ELECTRON_EFF_TIGHT 1
#define DELPHES_FLAG_ELECTRON_ISO_YES 0
#define DELPHES_FLAG_MUON_ISO_YES 0

/* event categories */
#define NONE 0
#define TWO_SS_ELEC 4
#define TWO_OS_ELEC 1

#define TWO_SS_MUON 5
#define TWO_OS_MUON 2
#define OS_ELEC_MUON 3
#define ONE_ELEC_P 6
#define ONE_ELEC_M 7
#define ONE_MUON_P 8
#define ONE_MUON_M 9
#define THREE_LEP 10
#define OTHER 11

/* stop decays */
#define CHRGINO1_NUTRINO1_P 1
#define CHRGINO1_NUTRINO1_N 2
#define CHRGINO1_NUTRINO2_P 3
#define CHRGINO1_NUTRINO2_N 4
#define CHRGINO2_NUTRINO1_P 5
#define CHRGINO2_NUTRINO1_N 6
#define CHRGINO2_NUTRINO2_P 7
#define CHRGINO2_NUTRINO2_N 8
#define NUTRINO1_TOP_P 9
#define NUTRINO1_TOP_N 10
#define NUTRINO2_TOP_P 11
#define NUTRINO2_TOP_N 12
#define OTHER_DECAY 13

#define DEBUG 1
#define DEBUG3 0
#define TRUTH_DEBUG 1
#define TRUTH_DEBUG1 0
#define OLD_DELPHES 0



void CreateHistograms(std::string histoName) {

  // Book histograms
  h_METNjets = new TH2D(histoName.c_str(), "MET-Njets histo", 20, 0, 200, 2, -0.5, 1.5);
  the_histos->Add(h_METNjets);

  h_nJets = new TH1F("nJet","Number of accepted Jets (all)",10,-0.5,9.5);
  the_histos->Add(h_nJets);

  h_nJets_emu = new TH1F("nJet_emu","Number of accepted Jets (emu events only)",10,-0.5,9.5);
  the_histos->Add(h_nJets_emu);

  h_nJets_overlap = new TH1F("nJet_overlap","Number of accepted Jets (after overlap removal)",10,-0.5,9.5);
  the_histos->Add(h_nJets_overlap);

  h_MET = new TH1F("MET","Transverse MET",100,0,500);
  the_histos->Add(h_MET);
  h_MET_emu = new TH1F("MET_emu","Transverse MET (emu pairs)",100,0,500);
  the_histos->Add(h_MET_emu);

  h_truth_MET_SumAll = new TH1F("truth_MET_SumAll","Truth Transverse MET from all truth visible stable particles",100,0,500);
  the_histos->Add(h_truth_MET_SumAll);
  h_truth_MET_SumNeutrals = new TH1F("truth_MET_SumNeutrals","Truth Transverse MET from all neutrals",100,0,500);
  the_histos->Add(h_truth_MET_SumNeutrals);

  h_truth_MET_SumAll_emu = new TH1F("truth_MET_SumAll_emu","Truth Transverse MET from all truth visible stable particles (emu pairs)",100,0,500);
  the_histos->Add(h_truth_MET_SumAll_emu);
  h_truth_MET_SumNeutrals_emu = new TH1F("truth_MET_SumNeutrals_emu","Truth Transverse MET from all neutrals (emu pairs)",100,0,500);
  the_histos->Add(h_truth_MET_SumNeutrals_emu);


  h_truth_METPhi_SumAll = new TH1F("truth_METPhi_SumAll","Phi of Truth Transverse MET from All truth visible stable particles",100,-6,6);
  h_truth_METPhi_SumNeutrals = new TH1F("truth_METPhi_SumNeutrals","Phi of Truth Transverse MET from Neutrals",100,-6,6);
  the_histos->Add(h_truth_METPhi_SumAll);
  the_histos->Add(h_truth_METPhi_SumNeutrals);

  h_JetPT = new TH1F("JetPT", "PT of accepted jets", 100, 0, 600);
  the_histos->Add(h_JetPT);
  h_ElecPT = new TH1F("ElecPT", "PT of accepted electrons", 100, 0, 400);
  the_histos->Add(h_ElecPT);
  h_MuonPT = new TH1F("MuonPT", "PT of accepted muons", 100, 0, 400);
  the_histos->Add(h_MuonPT);

  h_JetEta = new TH1F("JetEta", "Eta of accepted jets", 100, -5,5);
  the_histos->Add(h_JetEta);
  h_ElecEta = new TH1F("ElecEta", "Eta of accepted electrons", 100, -5,5);
  the_histos->Add(h_ElecEta);
  h_MuonEta = new TH1F("MuonEta", "Eta of accepted muons", 100, -5, 5);
  the_histos->Add(h_MuonEta);

  h_JetPhi = new TH1F("JetPhi", "Phi of accepted jets", 100, -6,6);
  the_histos->Add(h_JetPhi);
  h_ElecPhi = new TH1F("ElecPhi", "Phi of accepted electrons", 100, -6,6);
  the_histos->Add(h_ElecPhi);
  h_MuonPhi = new TH1F("MuonPhi", "Phi of accepted muons", 100, -6, 6);
  the_histos->Add(h_MuonPhi);



  h_JetPT_emu = new TH1F("JetPT_emu", "PT of accepted jets (emu pairs)", 100, 0, 600);
  the_histos->Add(h_JetPT_emu);
  h_ElecPT_emu = new TH1F("ElecPT_emu", "PT of accepted electrons (emu pairs)", 100, 0, 400);
  the_histos->Add(h_ElecPT_emu);
  h_MuonPT_emu = new TH1F("MuonPT_emu", "PT of accepted muons (emu pairs)", 100, 0, 400);
  the_histos->Add(h_MuonPT_emu);

  h_JetEta_emu = new TH1F("JetEta_emu", "Eta of accepted jets (emu pairs)", 100, -5,5);
  the_histos->Add(h_JetEta_emu);
  h_ElecEta_emu = new TH1F("ElecEta_emu", "Eta of accepted electrons (emu pairs)", 100, -5,5);
  the_histos->Add(h_ElecEta_emu);
  h_MuonEta_emu = new TH1F("MuonEta_emu", "Eta of accepted muons (emu pairs)", 100, -5, 5);
  the_histos->Add(h_MuonEta_emu);

  h_JetPhi_emu = new TH1F("JetPhi_emu", "Phi of accepted jets (emu pairs)", 100, -6,6);
  the_histos->Add(h_JetPhi_emu);
  h_ElecPhi_emu = new TH1F("ElecPhi_emu", "Phi of accepted electrons (emu pairs)", 100, -6,6);
  the_histos->Add(h_ElecPhi_emu);
  h_MuonPhi_emu = new TH1F("MuonPhi_emu", "Phi of accepted muons (emu pairs)", 100, -6, 6);
  the_histos->Add(h_MuonPhi_emu);

  h_decay = new TH1F("DecayChannels", "Different Decay channels of the stop", 14, -0.5, 13.5);
  h_decayfid = new TH1F("DecayChannelsFID", "Different Decay channels of the stop Fid", 14, -0.5, 13.5);
  h_decaymatched = new TH1F("DecayChannelsMatched", "Different Decay channels of the stop Matched", 14, -0.5, 13.5);
  the_histos->Add(h_decay);
  the_histos->Add(h_decayfid);
  the_histos->Add(h_decaymatched);

  h_truth_ElecPT = new TH1F("truth_ElecPT", "PT of truth prompt electrons", 100, 0, 400);
  the_histos->Add(h_truth_ElecPT);
  h_truth_MuonPT = new TH1F("truth_MuonPT", "PT of truth prompt muons", 100, 0, 400);
  the_histos->Add(h_truth_MuonPT);

  h_truth_MuonEta = new TH1F("truth_MuonETA", "ETA of truth prompt Muons", 100, -5, 5);
  h_truth_ElecEta = new TH1F("truth_ElecETA", "ETA of truth prompt Electrons", 100, -5, 5);
  the_histos->Add(h_truth_ElecEta);
  the_histos->Add(h_truth_MuonEta);
 
  h_truth_MuonPhi = new TH1F("truth_MuonPhi", "Phi of truth prompt Muons", 100, -6, 6);
  h_truth_ElecPhi = new TH1F("truth_ElecPhi", "Phi of truth prompt Electrons", 100, -6, 6);
  the_histos->Add(h_truth_ElecPhi);
  the_histos->Add(h_truth_MuonPhi);

  h_truth_ElecPTFid = new TH1F("truth_ElecPTFid", "PT of truth prompt electrons+Fid Cuts", 100, 0, 400);
  h_truth_MuonPTFid = new TH1F("truth_MuonPTFid", "PT of truth prompt muons+Fid Cuts", 100, 0, 400);
  the_histos->Add(h_truth_ElecPTFid);
  the_histos->Add(h_truth_MuonPTFid);

  h_truth_MuonEtaFid = new TH1F("truth_MuonETAFid", "ETA of truth prompt Muons+Fid Cuts", 100, -5, 5);
  h_truth_ElecEtaFid = new TH1F("truth_ElecETAFid", "ETA of truth prompt Electrons+Fid Cuts", 100, -5, 5);
  the_histos->Add(h_truth_ElecEtaFid);
  the_histos->Add(h_truth_MuonEtaFid);
 
  h_truth_MuonPhiFid = new TH1F("truth_MuonPhiFid", "Phi of truth prompt Muons+Fid Cuts", 100, -6, 6);
  h_truth_ElecPhiFid = new TH1F("truth_ElecPhiFid", "Phi of truth prompt Electrons+Fid Cuts", 100, -6, 6);
  the_histos->Add(h_truth_ElecPhiFid);
  the_histos->Add(h_truth_MuonPhiFid);

  h_truth_nElec = new TH1F("truth_nelecs", "Number of prompt truth electrons", 6, -0.5, 5.5); 
  the_histos->Add(h_truth_nElec);
  h_truth_nMuon = new TH1F("truth_nmuons", "Number of prompt truth muons", 6, -0.5, 5.5); 
  the_histos->Add(h_truth_nMuon);
  h_truth_nElecFid = new TH1F("truth_nelecsFid", "Number of prompt truth electrons+Fid Cuts", 6, -0.5, 5.5); 
  h_truth_nMuonFid = new TH1F("truth_nmuonsFid", "Number of prompt truth muons+ Fid Cuts", 6, -0.5, 5.5); 
  the_histos->Add(h_truth_nElecFid);
  the_histos->Add(h_truth_nMuonFid);

  h_lepclass_truth = new TH1F("truth_lepclass","Classification of leptonic events",12,-0.5,11.5);
  h_lepclass_fid = new TH1F("fid_lepclass","Classification of leptonic events",12,-0.5,11.5);
  h_lepclass_matched = new TH1F("matched_lepclass","Classification of leptonic events",12,-0.5,11.5);
  the_histos->Add(h_lepclass_truth);
  the_histos->Add(h_lepclass_fid);
  the_histos->Add(h_lepclass_matched);

  h_MElecPT = new TH1F("MElecPT", "PT of Truth Matched electrons", 100, 0, 400);
  h_MMuonPT = new TH1F("MMuonPT", "PT of Truth Matched muons", 100, 0, 400);
  the_histos->Add(h_MElecPT);
  the_histos->Add(h_MMuonPT);

  h_MMuonEta = new TH1F("MMuonETA", "ETA of Truth Matched ID Muons", 100, -5, 5);
  h_MElecEta = new TH1F("MElecETA", "ETA of Truth Matched ID Electrons", 100, -5, 5);
  the_histos->Add(h_MElecEta);
  the_histos->Add(h_MMuonEta);
 
  h_MMuonPhi = new TH1F("MMuonPhi", "Phi of Truth Matched ID Muons", 100, -6, 6);
  h_MElecPhi = new TH1F("MElecPhi", "Phi of Truth Matched ID Electrons", 100, -6, 6);
  the_histos->Add(h_MElecPhi);
  the_histos->Add(h_MMuonPhi);

  h_nMMuon = new TH1F("nMmuons", "Number of reco+id muons matched to truth", 6, -0.5, 5.5); 
  h_nMElec = new TH1F("nMelecs", "Number of reco+id electrons matched to truth", 6, -0.5, 5.5); 

  the_histos->Add(h_nMElec);
  the_histos->Add(h_nMMuon);


  // resolutions
  h_MElecPtRes = new TH1F("MElecPtRes", "(truthE - E)/truth E ", 100, -0.1, 0.1);
  the_histos->Add(h_MElecPtRes);
  
}

void CleanHistograms() {

  TH1F *histo;
  for (int k = 0; k < the_histos->GetEntries(); k++) {
    histo = (TH1F*)the_histos->At(k);
    delete histo;
  }
}


void WriteHistograms() {


  TH1F *histo;
  for (int k = 0; k < the_histos->GetEntries(); k++) {
    histo = (TH1F*)the_histos->At(k);
    histo->Write();
  }
}


void Fill_TruthLeptons_Histograms(TClonesArray *brTruth, std::vector<int> truth_electrons,std::vector<int> truth_muons ) {


    for (unsigned int j=0; j< truth_electrons.size(); j++) {
  
      GenParticle* p = (GenParticle*)brTruth->At(truth_electrons[j]);
      h_truth_ElecPT->Fill(p->PT);
      h_truth_ElecEta->Fill(p->Eta);
      h_truth_ElecPhi->Fill(p->Phi);
    }

    for (unsigned int j=0; j< truth_muons.size(); j++) {
  
      GenParticle* p = (GenParticle*)brTruth->At(truth_muons[j]);
      h_truth_MuonPT->Fill(p->PT);
      h_truth_MuonEta->Fill(p->Eta);
      h_truth_MuonPhi->Fill(p->Phi);
    }

    h_truth_nMuon->Fill(truth_muons.size());
    h_truth_nElec->Fill(truth_electrons.size());
}

void Fill_TruthLeptonsFid_Histograms(TClonesArray *brTruth, std::vector<int> truth_electrons,std::vector<int> truth_muons ) {

    for (unsigned int j=0; j< truth_electrons.size(); j++) {
  
      GenParticle* p = (GenParticle*)brTruth->At(truth_electrons[j]);

      h_truth_ElecPTFid->Fill(p->PT);
      h_truth_ElecEtaFid->Fill(p->Eta);
      h_truth_ElecPhiFid->Fill(p->Phi);
    }

    for (unsigned int j=0; j< truth_muons.size(); j++) {
  
      GenParticle* p = (GenParticle*)brTruth->At(truth_muons[j]);

      h_truth_MuonPTFid->Fill(p->PT);
      h_truth_MuonEtaFid->Fill(p->Eta);
      h_truth_MuonPhiFid->Fill(p->Phi);
    }
    h_truth_nMuonFid->Fill(truth_muons.size());
    h_truth_nElecFid->Fill(truth_electrons.size());
}


void Fill_MatchedLeptons_Histograms(TClonesArray *brElec,TClonesArray *brMuon, std::vector<int> matched_electrons, std::vector<int> matched_muons){

    for (unsigned int j=0; j< matched_electrons.size(); j++) {
  
      Electron* p = (Electron*)brElec->At(matched_electrons[j]);
      h_MElecPT->Fill(p->PT);
      h_MElecEta->Fill(p->Eta);
      h_MElecPhi->Fill(p->Phi);
    }

    for (unsigned int j=0; j< matched_muons.size(); j++) {
  
      Muon* p = (Muon*)brMuon->At(matched_muons[j]);
      h_MMuonPT->Fill(p->PT);
      h_MMuonEta->Fill(p->Eta);
      h_MMuonPhi->Fill(p->Phi);
    }
    h_nMMuon->Fill(matched_muons.size());
    h_nMElec->Fill(matched_electrons.size());
}

int Event_CategoryReco(TClonesArray *brElectron,TClonesArray *brMuon, std::vector<int> electrons, std::vector<int> muons) {
  // determine which category is the event - ee, mumu, emu, other

    int category = -1;
    int e_plus = 0;
    int e_minus = 0;
    int mu_plus = 0;
    int mu_minus = 0;
    for (unsigned int j=0; j< electrons.size(); j++) {
  
      Electron* p = (Electron*)brElectron->At(electrons[j]);
      if (p->Charge > 0)
        e_plus++;
      if (p->Charge < 0)
        e_minus++;
    }

    for (unsigned int j=0; j< muons.size(); j++) {
  
      Muon* p = (Muon*)brMuon->At(muons[j]);
      if (p->Charge > 0)
        mu_plus++;
      if (p->Charge < 0)
        mu_minus++;
    }

    category = Which_Event_Category(e_plus,e_minus,mu_plus,mu_minus);

    return category;
}




int Event_CategoryTruth(TClonesArray *brTruth, std::vector<int> electrons, std::vector<int> muons) {
  // determine which category is the event - ee, mumu, emu, other

    int category;
    int e_plus = 0;
    int e_minus = 0;
    int mu_plus = 0;
    int mu_minus = 0;
    for (unsigned int j=0; j< electrons.size(); j++) {
  
      GenParticle* p = (GenParticle*)brTruth->At(electrons[j]);
      if (p->Charge > 0)
        e_plus++;
      if (p->Charge < 0)
        e_minus++;
    }

    for (unsigned int j=0; j< muons.size(); j++) {
  
      GenParticle* p = (GenParticle*)brTruth->At(muons[j]);
      if (p->Charge > 0)
        mu_plus++;
      if (p->Charge < 0)
        mu_minus++;
    }


    category = Which_Event_Category(e_plus,e_minus,mu_plus,mu_minus);

    return category;
}

int Which_Event_Category(int e_plus, int e_minus, int mu_plus, int mu_minus) {

  int category = -1;

    if ((e_plus == 0) && (e_minus == 0) && (mu_plus == 0) && (mu_minus == 0)) 
      category = NONE;
    else if ((e_plus == 1) && (e_minus == 0) && (mu_plus == 0) && (mu_minus == 0))
      category = ONE_ELEC_P;
    else if ((e_plus == 0) && (e_minus == 1) && (mu_plus == 0) && (mu_minus == 0))
      category = ONE_ELEC_M;
    else if ((e_plus == 0) && (e_minus == 0) && (mu_plus == 1) && (mu_minus == 0))
      category = ONE_MUON_P;
    else if ((e_plus == 0) && (e_minus == 0) && (mu_plus == 0) && (mu_minus == 1))
      category = ONE_MUON_M;
    else if ((e_plus == 2) && (e_minus == 0) && (mu_plus == 0) && (mu_minus == 0))
      category = TWO_SS_ELEC;
    else if  ((e_plus == 0) && (e_minus == 2) && (mu_plus == 0) && (mu_minus == 0))
      category = TWO_SS_ELEC;
    else if ((e_plus == 1) && (e_minus == 1) && (mu_plus == 0) && (mu_minus == 0))
      category = TWO_OS_ELEC;
    else if  ((e_plus == 0) && (e_minus == 0) && (mu_plus == 2) && (mu_minus == 0))
      category = TWO_SS_MUON;
    else if  ((e_plus == 0) && (e_minus == 0) && (mu_plus == 0) && (mu_minus == 2))
      category = TWO_SS_MUON;
    else if  ((e_plus == 0) && (e_minus == 0) && (mu_plus == 1) && (mu_minus == 1))
      category = TWO_OS_MUON;
    else if ((e_plus == 1) && (e_minus == 0) && (mu_plus == 0) && (mu_minus == 1))
      category = OS_ELEC_MUON;
     else if ((e_plus == 0) && (e_minus == 1) && (mu_plus == 1) && (mu_minus == 0))
      category = OS_ELEC_MUON;
     else if ((e_plus > 2) || (e_minus > 2) || (mu_plus > 2) || (mu_minus > 2))
      category = THREE_LEP;
    else
      category = OTHER;

    return category;
 
}





void Truth_Prompt_LeptonFid(TClonesArray *brTruth, std::vector<int> prompt_e_truth, std::vector<int> prompt_mu_truth,
                             std::vector<int> &fid_e_truth,std::vector<int> &fid_mu_truth,bool explain) {

    fid_e_truth.clear();
    fid_mu_truth.clear();
    for (unsigned int j=0; j< prompt_e_truth.size(); j++) {
  
      GenParticle* p = (GenParticle*)brTruth->At(prompt_e_truth[j]);
      if (isGoodTruthElectron(p->PT,p->Eta,explain)) {
        fid_e_truth.push_back(prompt_e_truth[j]);
      }
    }

    for (unsigned int j=0; j< prompt_mu_truth.size(); j++) {
  
      GenParticle* p = (GenParticle*)brTruth->At(prompt_mu_truth[j]);
      if (isGoodTruthMuon(p->PT,p->Eta,explain)) {
        fid_mu_truth.push_back(prompt_mu_truth[j]);
      }
    }
}


bool isGoodTruthElectron(float pt, float eta, bool explain) {

   if (pt < ELECTRON_PT_MIN) {
    if (DEBUG && explain) std::cerr << "[isGoodElectron] PT=" << pt << " < " << ELECTRON_PT_MIN << std::endl;
    return false;
  }

  if (TMath::Abs(eta) > ELECTRON_ABSETA_MAX) {
    if (DEBUG && explain) std::cerr << "[isGoodElectron] |ETA|=" << TMath::Abs(eta) << " > " << ELECTRON_ABSETA_MAX << std::endl;
    return false;
  }

  /*
  if (TMath::Abs(eta) > ELECTRON_ABSETA_CRACKMIN &&
      TMath::Abs(eta) < ELECTRON_ABSETA_CRACKMAX) {
    if (DEBUG && explain) std::cerr << "[isGoodElectron] |ETA|=" << TMath::Abs(eta) << " in crack" << std::endl;
    return false;
  }
  */

  return true;
}


bool isGoodTruthMuon(float pt, float eta, bool explain) {

  if (pt < MUON_PT_MIN) {
    if (DEBUG && explain) std::cerr << "[isGoodMuon] PT=" << pt << " < " << MUON_PT_MIN << std::endl;
    return false;
  }

  if (TMath::Abs(eta) > MUON_ABSETA_MAX) {
    if (DEBUG && explain) std::cerr << "[isGoodMuon] |ETA|=" << TMath::Abs(eta) << " > " << MUON_ABSETA_MAX << std::endl;
    return false;
  }
  return true;
}

bool isGoodElectron(Electron* electron, bool explain)
{
  if (electron->PT < ELECTRON_PT_MIN) {
    if (DEBUG && explain) std::cerr << "[isGoodElectron] PT=" << electron->PT << " < " << ELECTRON_PT_MIN << std::endl;
    return false;
  }

  /*
  if (! OLD_DELPHES) {
    if (! electron->EfficiencyFlags & (1 << DELPHES_FLAG_ELECTRON_EFF_TIGHT)) {
      if (DEBUG && explain) std::cerr << "[isGoodElectron] EfficiencyFlags=" << electron->EfficiencyFlags << std::endl;
      return false;
    }
  }
  */

  if (TMath::Abs(electron->Eta) > ELECTRON_ABSETA_MAX) {
    if (DEBUG && explain) std::cerr << "[isGoodElectron] |ETA|=" << TMath::Abs(electron->Eta) << " > " << ELECTRON_ABSETA_MAX << std::endl;
    return false;
  }

  /*
  if (TMath::Abs(electron->Eta) > ELECTRON_ABSETA_CRACKMIN &&
      TMath::Abs(electron->Eta) < ELECTRON_ABSETA_CRACKMAX) {
    if (DEBUG && explain) std::cerr << "[isGoodElectron] |ETA|=" << TMath::Abs(electron->Eta) << " in crack" << std::endl;
    return false;
  }

  
  if (! OLD_DELPHES) {
    if (! electron->IsolationFlags & (1 << DELPHES_FLAG_ELECTRON_ISO_YES)) {
      if (DEBUG && explain) std::cerr << "[isGoodElectron] IsolationFlags=" << electron->IsolationFlags << std::endl;
      return false;
    }
  }
  */

  return true;
}

bool isGoodMuon(Muon *muon, bool explain)
{
  if (muon->PT < MUON_PT_MIN) {
    if (DEBUG && explain) std::cerr << "[isGoodMuon] PT=" << muon->PT << " < " << MUON_PT_MIN << std::endl;
    return false;
  }

  if (TMath::Abs(muon->Eta) > MUON_ABSETA_MAX) {
    if (DEBUG && explain) std::cerr << "[isGoodMuon] |ETA|=" << TMath::Abs(muon->Eta) << " > " << MUON_ABSETA_MAX << std::endl;
    return false;
  }

  /*
  if (! OLD_DELPHES) {
    if (! muon->IsolationFlags & (1 << DELPHES_FLAG_MUON_ISO_YES)) {
      if (DEBUG && explain) std::cerr << "[isGoodMuon] IsolationFlags=" << muon->IsolationFlags << std::endl;
      return false;
    }
  }
  */

  return true;
}

bool isGoodJet(Jet *jet, bool explain)
{
  if (jet->PT < JET_PT_MIN) {
    if (DEBUG && explain) std::cerr << "[isGoodJet] PT=" << jet->PT << " < " << JET_PT_MIN << std::endl;
    return false;
  }

  if (TMath::Abs(jet->Eta) > JET_ABSETA_MAX) {
    if (DEBUG && explain) std::cerr << "[isGoodJet] |ETA|=" << TMath::Abs(jet->Eta) << " > " << JET_ABSETA_MAX << std::endl;
    return false;
  }

  return true;
}

void MET_Truth_SumAll(TClonesArray *brTruth, float &truth_met, float &truth_met_phi) {
  // this function adds the momentum of all the visible truth stable particles = status == 1
  // and from x and y components the 

  float px=0;
  float py=0;
  truth_met_phi = 0;
  truth_met = 0;
  for (Int_t iTruth = 0; iTruth < brTruth->GetEntries(); ++iTruth) {

    GenParticle* p = (GenParticle*)brTruth->At(iTruth);

    if (p->Status == 1) {
      switch (p->PID) {// do not count the neutrals

        case 12: // electron neutrino
        case -12:
        case 14: // muon neutrino
        case -14:
        case 16: // tau neutrino
        case -16:
        case 	1000022: // neutralino 1
        case 	1000023: // neutralino 2
          break;
        default: 
          px += p->Px;
          py += p->Py;

    //      std::cout << " All P " << p->PID << " " << p->Status << " " << p->Px << "  py " << p->Py << std::endl;
      }
    }
  }

  //std::cout << "SUm px " << px << "  py " << py << std::endl;
  truth_met = TMath::Sqrt(px*px + py*py);
  truth_met_phi=atan2(-py,-px);
}

void MET_Truth_SumNeutrals(TClonesArray *brTruth, float &truth_met, float &truth_met_phi) {
  // this function adds the momentum of all the stable particles = status == 1
  // and from x and y components the 

  float px=0;
  float py=0;
  truth_met_phi = 0;
  truth_met = 0;
  for (Int_t iTruth = 0; iTruth < brTruth->GetEntries(); ++iTruth) {

    GenParticle* p = (GenParticle*)brTruth->At(iTruth);

    switch ((int)p->PID) {

      case 12: // electron neutrino
      case -12:
      case 14: // muon neutrino
      case -14:
      case 16: // tau neutrino
      case -16:
      case 	1000022: // neutralino 1
      case 	1000023: // neutralino 2

    //    std::cout << " Neutral P " << p->Px << "  py " << p->Py << std::endl;
        px += p->Px;
        py += p->Py;
      break;
      default: ;
    }
  }

  //std::cout << " Neutrals px " << px << "  py " << py << std::endl;
  truth_met = TMath::Sqrt(px*px + py*py);
  truth_met_phi=atan2(py,px);
}




void PrintTruthParticle(std::string prefix,GenParticle *part, int index) {

  std::cout << prefix << " " << index << " pid " << part->PID << " pt " << part->PT << " " << part->Eta << " " << part->Phi << " mass " << part->Mass << " D1 D2 " << part->D1 << " " << part->D2 << " M1 M2 " << part->M1 << " " << part->M2 << std::endl;
}

void PrintList_TruthParticles(std::string prefix,TClonesArray *brTruth,std::vector<int> leptons, std::vector<int> bosons) {

    for (unsigned int j=0; j< leptons.size(); j++) {
      GenParticle* pl = (GenParticle*)brTruth->At(leptons[j]);
      GenParticle* pb = (GenParticle*)brTruth->At(bosons[j]);
      PrintTruthParticle(prefix,pl,leptons[j]);
      PrintTruthParticle(prefix,pb,bosons[j]);
    }
}




int Truth_Prompt_Lepton(TClonesArray *brTruth, int PID,std::vector<int> &truth_leptons,std::vector<int> &truth_bosons ) {

  //std::vector<int> truth_leptons,truth_bosons;

  truth_leptons.empty();
  truth_bosons.empty();

  for (Int_t iTruth = 0; iTruth < brTruth->GetEntries(); ++iTruth)
  {
    GenParticle *part = (GenParticle*) brTruth->At(iTruth);
    //std::cout << "particle " << part->PID << std::endl;
    //std::cout << iTruth << std::endl;
    
    if ((fabs(part->PID) == PID)) {


      if (TRUTH_DEBUG) PrintTruthParticle("prompt_test ",part,iTruth);
      if (isLeptonPrompt(brTruth,iTruth,truth_bosons)) {// it is prompt
        truth_leptons.push_back(iTruth);
      }
    }
  }

  /*
  // check that the boson is from either a top, stop or chargino
  for (unsigned int j=0; j< truth_leptons.size(); j++) {


    GenParticle* pb = (GenParticle*)brTruth->At(truth_bosons[j]);
    GenParticle* gb = (GenParticle*)brTruth->At(pb->M1);

    while ((fabs(gb->PID) == 24)) {
      gb = (GenParticle*)brTruth->At(gb->M1);
    }

    if ((fabs(gb->PID) != 6) && (fabs(gb->PID) < 1000000)) {
      std::cout << " removing prompt lepton Woson grand parent is " << gb->PID << std::endl;
      truth_leptons.erase(truth_leptons.begin()+j);
    }
    */

    /*
    GenParticle* pl = (GenParticle*)brTruth->At(truth_leptons[j]);


    std::cout << " prompt lepton found: " <<std::endl;
    PrintTruthParticle(pl);
    PrintTruthParticle(pb);

    GenParticle* gb = (GenParticle*)brTruth->At(pb->M1);
    PrintTruthParticle(gb);
    
  }
  */

  return 0;
}

 



bool isLeptonPrompt(TClonesArray *brTruth, int ipos,std::vector<int> &boson_i) {

  int W_pid = 24;
  int Z_pid = 23;
  int chargino_pid = 1000024;
  int k= ipos;
  GenParticle *part = (GenParticle*)brTruth->At(k);

  return true;

  /*
  if (TRUTH_DEBUG) printf("Testing whether %d %d is prompt\n",k,part->PID);

  GenParticle *mother = (GenParticle*)brTruth->At(part->M1);
  if ((fabs(mother->PID) == W_pid) || ((fabs(mother->PID)) == Z_pid)|| ((fabs(mother->PID)) == chargino_pid)) {
      if (TRUTH_DEBUG) printf("Prompt! %d %d\n",ipos,part->PID);
      boson_i.push_back(part->M1);
      return true;
  }
  else {
    return false;
  }
  */



  /*
  while ((fabs(part->PID) != W_pid) && (fabs(part->PID) != Z_pid) && (k != -1)) {

    k = part->M1;
    if (k != -1) {
      part = (GenParticle*) brTruth->At(k);
      printf("Testing up %d %d is prompt\n",k,part->PID);
    }
  }
  // we found the position where the particle is not PID

  if (k != -1) {
    std::cout << " it is prompt" << std::endl;
    // now check whether it is directly after the boson.. not the last one
    GenParticle *part2 = (GenParticle*)brTruth->At(ipos);

    PrintTruthParticle("prompt_test ",part2,ipos);
    if ((int)part2->D1 != -1) {

      printf("Prompt! %d %d\n",ipos,part2->PID);

      boson_i.push_back(k);
      return true; // parent of the lepton is a boson - making it prompt
    }
    else
      return false;
  }
  else 
  {
    //std::cout << " not prompt" << std::endl;
    return false; // the lepton is not a prompt lepton
  }
  */

}

int FindLastParticle(TClonesArray *brTruth, int ipos, int PID) {

  int k= ipos;
  GenParticle *part = (GenParticle*)brTruth->At(k);
  while (fabs(part->PID) == PID) {

    k = part->D1;
    part = (GenParticle*) brTruth->At(k);
  }
  // we found the position where the particle is not PID
  return (part->M1);
}


void DecayChannel(TClonesArray *brTruth, std::vector<int> & decays) {

  // lets find the stop that decays
  // then find chargino that decays
  int stop_pid = 1000006;
  int chargino1_pid = 1000024;
  int chargino2_pid = 1000025;
  int neutralino1_pid = 1000022;
  int neutralino2_pid = 1000023;

  int b_pid = 5;
  int t_pid = 6;
  int stop_pos;
  int nstop_pos;

  bool found_stop = false;
  bool found_nstop = false;
  bool is_chargino1, is_chargino2;

  decays.clear();
  int nstops = 0; // need to find two
  for (Int_t iTruth = 0; iTruth < brTruth->GetEntries(); ++iTruth)
  {
    GenParticle *part = (GenParticle*) brTruth->At(iTruth);
    //std::cout << "particle " << part->PID << std::endl;
    //std::cout << iTruth << std::endl;
    
    if ((fabs(part->PID) == stop_pid) && (found_stop == false)) {
      nstops++;
      stop_pos = FindLastParticle(brTruth,iTruth,stop_pid);
      if (TRUTH_DEBUG1) PrintTruthParticle(" ",(GenParticle*)brTruth->At(stop_pos),stop_pos);
//      std::cout << " Hi stop" << std::endl;
      //found_stop = true;

      GenParticle *stop_part = (GenParticle*) brTruth->At(stop_pos);
      // daughters 1 and 2
      GenParticle *part1 = (GenParticle*) brTruth->At(stop_part->D1);
      GenParticle *part2 = (GenParticle*) brTruth->At(stop_part->D2);
      int chargino_pos = -1;
      if ((fabs(part1->PID) == chargino1_pid) || (fabs(part1->PID) == chargino2_pid)) {
        chargino_pos = stop_part->D1;
        if (fabs(part1->PID) == chargino1_pid) 
          is_chargino1 = true;
        else
          is_chargino2 = true;
      }
      else if ((fabs(part2->PID) == chargino1_pid) || (fabs(part2->PID) == chargino2_pid)) {
        chargino_pos = stop_part->D2;
        if (fabs(part1->PID) == chargino1_pid) 
          is_chargino1 = true;
        else
          is_chargino2 = true;
      }

      if (chargino_pos != -1) {
        if (is_chargino1)
          chargino_pos = FindLastParticle(brTruth,chargino_pos,chargino1_pid);
        else
          chargino_pos = FindLastParticle(brTruth,chargino_pos,chargino2_pid);

        if (TRUTH_DEBUG1) PrintTruthParticle(" ",(GenParticle*) brTruth->At(chargino_pos),chargino_pos);

        GenParticle *chargino_part = (GenParticle*) brTruth->At(chargino_pos);
        GenParticle *d1 = (GenParticle*) brTruth->At(chargino_part->D1);
        GenParticle *d2 = (GenParticle*) brTruth->At(chargino_part->D2);

        if (chargino_part->Charge > 0) {
          if (is_chargino1) {
            if (d1->PID == neutralino1_pid)
              decays.push_back(CHRGINO1_NUTRINO1_P);
            else if (d1->PID == neutralino2_pid)
              decays.push_back(CHRGINO1_NUTRINO2_P);
            else if (d2->PID == neutralino1_pid)
              decays.push_back(CHRGINO1_NUTRINO1_P);
            else if (d2->PID == neutralino2_pid)
              decays.push_back(CHRGINO1_NUTRINO2_P);
            else
            {
              PrintTruthParticle(" ",d1,chargino_part->D1);
              PrintTruthParticle(" ",d2,chargino_part->D2);
              decays.push_back(OTHER_DECAY);
            }
          }
          else if  (is_chargino2) {
            if (d1->PID == neutralino1_pid)
              decays.push_back(CHRGINO2_NUTRINO1_P);
            else if (d1->PID == neutralino2_pid)
              decays.push_back(CHRGINO2_NUTRINO2_P);
            else if (d2->PID == neutralino1_pid)
              decays.push_back(CHRGINO2_NUTRINO1_P);
            else if (d2->PID == neutralino2_pid)
              decays.push_back(CHRGINO2_NUTRINO2_P);
            else {
              PrintTruthParticle("other decay  ",d1,chargino_part->D1);
              PrintTruthParticle("other decay  ",d2,chargino_part->D2);
              decays.push_back(OTHER_DECAY);
            }
          }
        }
        else{
//          std::cout << " here first " << d1->PID <<  " " << d2->PID << std::endl;
          if (is_chargino1) {
            if (d1->PID == neutralino1_pid) {
//              std::cout << " here aad " << std::endl;
              decays.push_back(CHRGINO1_NUTRINO1_N);
            }
            else if (d1->PID == neutralino2_pid)
              decays.push_back(CHRGINO1_NUTRINO2_N);
            else if (d2->PID == neutralino1_pid)
              decays.push_back(CHRGINO1_NUTRINO1_N);
            else if (d2->PID == neutralino2_pid)
              decays.push_back(CHRGINO1_NUTRINO2_N);
            else {
              PrintTruthParticle("other decay ",d1,chargino_part->D1);
              PrintTruthParticle("other decay ",d2,chargino_part->D2);
              decays.push_back(OTHER_DECAY);
            }
          }
          else if  (is_chargino2) {
            if (d1->PID == neutralino1_pid)
              decays.push_back(CHRGINO2_NUTRINO1_N);
            else if (d1->PID == neutralino2_pid)
              decays.push_back(CHRGINO2_NUTRINO2_N);
            else if (d2->PID == neutralino1_pid)
              decays.push_back(CHRGINO2_NUTRINO1_N);
            else if (d2->PID == neutralino2_pid)
              decays.push_back(CHRGINO2_NUTRINO2_N);
            else {
              PrintTruthParticle("other decay  ",d1,chargino_part->D1);
              PrintTruthParticle("other decay  ",d2,chargino_part->D2);
              decays.push_back(OTHER_DECAY);
            }
          }
        }
      }
      else {
        if ((fabs(part1->PID) == neutralino1_pid) || (fabs(part2->PID) == neutralino1_pid)) {
          if (stop_part->Charge > 0) 
            decays.push_back(NUTRINO1_TOP_P);
          else 
            decays.push_back(NUTRINO1_TOP_N);
        }
        else if  ((fabs(part1->PID) == neutralino2_pid) || (fabs(part2->PID) == neutralino2_pid)) {
          if (stop_part->Charge > 0) 
            decays.push_back(NUTRINO2_TOP_P);
          else 
            decays.push_back(NUTRINO2_TOP_N);
        }
        else {

          PrintTruthParticle("other decay  ",part1,stop_part->D1);
          PrintTruthParticle("other decay  ",part2,stop_part->D2);
          decays.push_back(OTHER_DECAY);
        }


        /*
        if (TRUTH_DEBUG1) PrintTruthParticle("truth decays ",(GenParticle*) brTruth->At(chargino_part->D1),chargino_part->D1);
        if (TRUTH_DEBUG1) PrintTruthParticle("truth decays ",(GenParticle*) brTruth->At(chargino_part->D2),chargino_part->D2);

          GenParticle *next_part = (GenParticle*) brTruth->At(chargino_part->D1);

        // check for daighters
        if (TRUTH_DEBUG1)  {
          GenParticle *next_part = (GenParticle*) brTruth->At(chargino_part->D1);
          if (next_part->D1 != -1)
            PrintTruthParticle("truth decays daughters 1 D1 ",(GenParticle*) brTruth->At(next_part->D1),next_part->D1);

          if (next_part->D2 != -1)
            PrintTruthParticle("truth decays daughters 1 D2 ",(GenParticle*) brTruth->At(next_part->D2),next_part->D2);

          next_part = (GenParticle*) brTruth->At(chargino_part->D2);
          if (next_part->D1 != -1)
            PrintTruthParticle("truth decays daughters 2 D1 ",(GenParticle*) brTruth->At(next_part->D1),next_part->D1);

          if (next_part->D2 != -1)
            PrintTruthParticle("truth decays daughters 2 D2 ",(GenParticle*) brTruth->At(next_part->D2),next_part->D2);
        }
        */
      }
    }

    if (nstops ==2) {
      found_stop = true;
      break;
    }
    
//    if ((part->PID == (-1.0*stop_pid)) && (found_nstop == false)) {
//      nstop_pos = FindLastParticle(brTruth,iTruth,-stop_pid);


//      std::cout << " Hi nstop" << std::endl;
//      found_nstop = true;
  }


//    if ((found_stop) && (found_nstop)) { 
//      break;
//    }
    
//  }
  
}

//int main(
//  std::string inputFile,
//  std::string outputFile,
//  Float_t norm,
//  std::string histoName
//  )

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


int main(int argc, char *argv[])
{

  std::cout << " argc " << argc << std::endl;
  if (argc != 4) {
    std::cout << " error " << argc << " " << argv[1] << " "<< argv[2]<< " "<< argv[3]  << std::endl;
    return 1;
  }

  std::vector<std::string> input_delphes_files;
  std::string inputFiles(argv[1]);
 // now the files delimited by '\0'
 
  input_delphes_files = split(inputFiles,',');

  std::string outputFile(argv[2]);
  /*
  Float_t norm;
  norm = atof(argv[3]);
  */

  std::string histoName(argv[3]);

//  printf("%s %s %.2f %s\n",argv[1],argv[2],atof(argv[3]),argv[4]);
  printf("reading %s histo_labels %s\n",inputFiles.c_str(),histoName.c_str());

  //std::cout << inputFile << std::endl;
  // Load input data into a TChain object
  TChain chain("Delphes");
  for (unsigned int l=0;l<input_delphes_files.size();l++) {
    chain.Add(input_delphes_files[l].c_str());
  }

  // Create pointers for tree and branch reading
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  TClonesArray *brJet;
  TClonesArray *brElec; 
  TClonesArray *brMuon; 
  TClonesArray *brMET;
  TClonesArray *brTruth;

  Long64_t nEvents = treeReader->GetEntries();
  // if (nEvents > 1000) nEvents = 1000; // XXX
  //if (OLD_DELPHES) {
    brJet = treeReader->UseBranch("Jet");
    brElec = treeReader->UseBranch("Electron");
    brMuon = treeReader->UseBranch("Muon");
    brMET = treeReader->UseBranch("MissingET");
    brTruth = treeReader->UseBranch("Particle");
  /*
  } else {
    brJet = treeReader->UseBranch("Jet");
    brElec = treeReader->UseBranch("ElectronATLAS");
    brMuon = treeReader->UseBranch("MuonATLAS");
    brMET = treeReader->UseBranch("MissingETATLAS");
    brTruth = treeReader->UseBranch("Particle");
  }
  */


  // Book histograms
  the_histos = new TObjArray();
  CreateHistograms(histoName);

  // Open output file for writing
  TFile f(outputFile.c_str(), "recreate");

  // Initialise variables
  int tallyFilledEvents = 0;
  std::set<int> goodJets, goodMuons, goodElectrons;
  Double_t DeltaR;
  float truth_met_all,truth_met_phi_all;
  float truth_met_neutrals,truth_met_phi_neutrals;
  Float_t MET;

  // Loop over events
  //nEvents = 5;
  for (Int_t evt = 0; evt < nEvents; ++evt)
  {
    // Load selected branches for current event
    if (DEBUG) std::cerr << "Loading event #" << evt << "..." << std::endl;
    treeReader->ReadEntry(evt);
    if (DEBUG) std::cerr << "Object candidate counts: Jets=" << brJet->GetEntries() << " Electrons=" << brElec->GetEntries() << " Muons=" << brMuon->GetEntries() << std::endl;

    goodJets.clear();
    goodElectrons.clear();
    goodMuons.clear();
    int tallyPosElectrons = 0;
    int tallyNegElectrons = 0;
    int tallyPosMuons = 0;
    int tallyNegMuons = 0;

    // Check MET exists
    // if (brMET->GetEntriesFast() < 1) {
    //   if (DEBUG) std::err << "Rejected event #" << evt << ": no MissingEt object" << std::endl;
    //   continue;
    // }

    /*
    std::vector<int> stop_decays;
    if (DEBUG3) std::cout << evt << std::endl;
    DecayChannel(brTruth,stop_decays);
    if (stop_decays.size() == 2) {
      h_decay->Fill(stop_decays[0]);
      h_decay->Fill(stop_decays[1]);
    }
    else {
      std::cout << " Problem with the stop decays " << std::endl;
    }
    */

    std::vector<int> prompt_mu_truth,prompt_e_truth;
    std::vector<int> fid_mu_truth,fid_e_truth;
    std::vector<int> mubosons_truth,ebosons_truth;

    std::vector<int> matched_elecs,matched_muons;

    // Find the truth prompt electron and muon - either decay from W from top or stop or chargino
    Truth_Prompt_Lepton(brTruth,11,prompt_e_truth,ebosons_truth);
    Truth_Prompt_Lepton(brTruth,13,prompt_mu_truth,mubosons_truth);

    // Find the truth prompt leptons that pass the fiducial cuts
    Truth_Prompt_LeptonFid(brTruth,prompt_e_truth,prompt_mu_truth,fid_e_truth,fid_mu_truth,true);

//    PrintList_TruthParticles("truth_e    ",brTruth,prompt_e_truth_p,ebosons_truth_p);
//    PrintList_TruthParticles("e_minus",brTruth,prompt_e_truth_m,ebosons_truth_m);

//    PrintList_TruthParticles("truth_muon ",brTruth,prompt_mu_truth_p,mubosons_truth_p);
//    PrintList_TruthParticles("mu_minus",brTruth,prompt_mu_truth_m,mubosons_truth_m);

    Fill_TruthLeptons_Histograms(brTruth,prompt_e_truth,prompt_mu_truth);
    Fill_TruthLeptonsFid_Histograms(brTruth,fid_e_truth,fid_mu_truth);

    /*
    // the truth MET
    MET_Truth_SumAll(brTruth,truth_met_all,truth_met_phi_all);
    h_truth_MET_SumAll->Fill(truth_met_all); 
    h_truth_METPhi_SumAll->Fill(truth_met_phi_all); 
    MET_Truth_SumNeutrals(brTruth,truth_met_neutrals,truth_met_phi_neutrals);
    h_truth_MET_SumNeutrals->Fill(truth_met_neutrals); 
    h_truth_METPhi_SumNeutrals->Fill(truth_met_phi_neutrals); 
    */

    // Look for acceptable jets (meet PT, eta requirements and do not overlap 
    // with any good electrons or good muons)
    h_nJets->Fill(brJet->GetEntries());
    for (Int_t iJet = 0; iJet < brJet->GetEntries(); ++iJet)
    {
      bool jetRejected = false;

      Jet *jet = (Jet*) brJet->At(iJet);
      if (! isGoodJet(jet, true)) {
        if (DEBUG) std::cerr << "Rejected jet #" << iJet << ": isGoodJet" << std::endl;
        continue;
      }

      // only check against prompt leptons

      for (Int_t iElec = 0; iElec < brElec->GetEntriesFast(); ++iElec) {

        Electron *elec = (Electron*) brElec->At(iElec);
        if (calcDeltaR(elec->Phi, jet->Phi, elec->Eta, jet->Eta) < DR_ELECTRON_JET_REMOVE_JET) {
          if (DEBUG) std::cerr << "Rejected jet #" << iJet << ": DR_ELECTRON_JET_REMOVE_JET" << std::endl;
          jetRejected = true;
          break;
        }
      }
      if (jetRejected) continue;

      for (Int_t iMuon = 0; iMuon < brMuon->GetEntriesFast(); ++iMuon) {

        Muon *muon = (Muon*) brMuon->At(iMuon);
        if (calcDeltaR(muon->Phi, jet->Phi, muon->Eta, jet->Eta) < DR_MUON_JET_REMOVE_JET) {
          if (DEBUG) std::cerr << "Rejected jet #" << iJet << ": DR_MUON_JET_REMOVE_JET" << std::endl;
          jetRejected = true;
          break;
        }
      }
      if (jetRejected) continue;

      goodJets.insert(iJet);
    } // end loop over jets

    // Look for good muons (meet PT, eta requirements and not within vicinity
    // of unvetoed jet)
    for (Int_t iMuon = 0; iMuon < brMuon->GetEntriesFast(); ++iMuon) {
      Muon *muon = (Muon*) brMuon->At(iMuon);
      if (! isGoodMuon(muon, true)) {
        if (DEBUG) std::cerr << "Rejected muon #" << iMuon << ": isGoodMuon" << std::endl;
        continue;
      }

      // Reset closest distance between jet and muon
      Double_t closestDist = TMath::Infinity();

      for (std::set<int>::iterator itJet=goodJets.begin(); itJet != goodJets.end(); ++itJet) {
        Jet *jet = (Jet*) brJet->At(*itJet);
        DeltaR = calcDeltaR(muon->Phi, jet->Phi, muon->Eta, jet->Eta);
        if (DeltaR < closestDist) closestDist = DeltaR;
      }

      if (closestDist < DR_MUON_JET_REMOVE_MUON) {
        if (DEBUG) std::cerr << "Rejected muon #" << iMuon << ": DR_MUON_JET_REMOVE_MUON" << std::endl;
        continue;
      }

      goodMuons.insert(iMuon);
    } // end loop over muons

    // Look for good electrons (meet PT, eta requirements and are not within
    // vicinity of unvetoed jet)
    for (Int_t iElec = 0; iElec < brElec->GetEntriesFast(); ++iElec) {
      Electron *elec = (Electron*) brElec->At(iElec);
      if (! isGoodElectron(elec, true)) {
        if (DEBUG) std::cerr << "Rejected electron #" << iElec << ": isGoodElectron" << std::endl;
        continue;
      }

      // Reset closest distance between jet and electron
      Double_t closestDist = TMath::Infinity();

      for (std::set<int>::iterator itJet=goodJets.begin(); itJet != goodJets.end(); ++itJet) {
        Jet *jet = (Jet*) brJet->At(*itJet);
        DeltaR = calcDeltaR(elec->Phi, jet->Phi, elec->Eta, jet->Eta);
        if (DeltaR < closestDist) closestDist = DeltaR;
      }

      if (closestDist < DR_ELECTRON_JET_REMOVE_ELECTRON) {
        if (DEBUG) std::cerr << "Rejected electron #" << iElec << ": DR_ELECTRON_JET_REMOVE_ELECTRON " << closestDist << " < " << DR_ELECTRON_JET_REMOVE_ELECTRON << std::endl;
        continue;
      }

      goodElectrons.insert(iElec);
    } // end loop over electrons
 
    // MET for all the events
    MET = ((MissingET*) brMET->At(0))->MET;
    h_MET->Fill(MET);
    // lets histograms the jets,electrons and muons of all events
    h_nJets_overlap->Fill(goodJets.size());
    for (std::set<int>::iterator itJet=goodJets.begin(); itJet != goodJets.end(); ++itJet) {
      Jet *jet = (Jet*) brJet->At(*itJet);
      h_JetPT->Fill(jet->PT);
      h_JetEta->Fill(jet->Eta);
      h_JetPhi->Fill(jet->Phi);
    }

    for (std::set<int>::iterator itElec=goodElectrons.begin(); itElec != goodElectrons.end(); ++itElec) {
      Electron *elec = (Electron*) brElec->At(*itElec);
      h_ElecPT->Fill(elec->PT);
      h_ElecEta->Fill(elec->Eta);
      h_ElecPhi->Fill(elec->Phi);
    }

    for (std::set<int>::iterator itMuon=goodMuons.begin(); itMuon != goodMuons.end(); ++itMuon) {
      Muon *muon = (Muon*) brMuon->At(*itMuon);
      h_MuonPT->Fill(muon->PT);
      h_MuonEta->Fill(muon->Eta);
      h_MuonPhi->Fill(muon->Phi);
    }

// lets match prompt truth leptons to reco+id leptons
    matched_elecs.clear();
    matched_muons.clear();

    for (unsigned int l=0;l< fid_e_truth.size();l++) {

      GenParticle *truth_el = (GenParticle*) brTruth->At(fid_e_truth[l]);
      for (std::set<int>::iterator it=goodElectrons.begin(); it != goodElectrons.end(); ++it) {
        Electron *elec = (Electron*) brElec->At(*it);
        if (calcDeltaR(elec->Phi, truth_el->Phi, elec->Eta, truth_el->Eta) < DR_ELECTRON_JET_REMOVE_JET) {
          matched_elecs.push_back(*it);

          if (truth_el->E > 0) {
            double px = elec->PT*TMath::Cos(elec->Phi);
            double py = elec->PT*TMath::Sin(elec->Phi);
            double pz = elec->PT*TMath::SinH(elec->Eta);
            
            double E = sqrt(px*px + py*py + pz*pz);
            h_MElecPtRes->Fill((truth_el->E - E)/truth_el->E);
              }
          if (elec->Charge < 0) {
            tallyNegElectrons++;
          } else {
            tallyPosElectrons++;
          }
          break;
        }
      }
    }

    for (unsigned int l=0;l< fid_mu_truth.size();l++) {

      GenParticle *truth_mu = (GenParticle*) brTruth->At(fid_mu_truth[l]);

      for (std::set<int>::iterator itMuon=goodMuons.begin(); itMuon != goodMuons.end(); ++itMuon) {
        Muon *muon = (Muon*) brMuon->At(*itMuon);
        if (calcDeltaR(muon->Phi, truth_mu->Phi, muon->Eta, truth_mu->Eta) < DR_MUON_JET_REMOVE_JET) {
//          std::cout << " mu index " << *itMuon << std::endl;
          matched_muons.push_back(*itMuon);
          if (muon->Charge < 0) {
            tallyNegMuons++;
          } else {
            tallyPosMuons++;
          }
          break;
        }
      }
    }
    //std::cout << "n matched " << matched_elecs.size() << " muons " << matched_muons.size() << std::endl;
    Fill_MatchedLeptons_Histograms(brElec,brMuon,matched_elecs,matched_muons);

    /*
    int truth_category = Event_CategoryTruth(brTruth,prompt_e_truth,prompt_mu_truth);
    h_lepclass_truth->Fill(truth_category);
    int fid_category = Event_CategoryTruth(brTruth,fid_e_truth,fid_mu_truth);
    h_lepclass_fid->Fill(fid_category);
    if (fid_category == OS_ELEC_MUON) {
      if (stop_decays.size() == 2) {
        h_decayfid->Fill(stop_decays[0]);
        h_decayfid->Fill(stop_decays[1]);
      }
    }
    int reco_category = Event_CategoryReco(brElec,brMuon,matched_elecs,matched_muons);
    h_lepclass_matched->Fill(reco_category);
    if (reco_category == OS_ELEC_MUON) {
      if (stop_decays.size() == 2) {
        h_decaymatched->Fill(stop_decays[0]);
        h_decaymatched->Fill(stop_decays[1]);
      }
    }

    // Throw out events without opposite-sign electron-muon pairs
    if (! ((tallyNegElectrons == 1 && tallyPosMuons == 1)
        || (tallyPosElectrons == 1 && tallyNegMuons == 1))) {
      if (DEBUG) {
        if (tallyNegElectrons + tallyNegMuons
          + tallyPosElectrons + tallyPosMuons != 2) {
          std::cerr << "Rejected event #" << evt << ": LEPTON_COUNT ";
        } else if ((tallyNegElectrons + tallyPosElectrons > 1)
                || (tallyNegMuons     + tallyPosMuons     > 1)) {
          std::cerr << "Rejected event #" << evt << ": LEPTON_FLAVOUR ";
        } else {
          std::cerr << "Rejected event #" << evt << ": LEPTON_CHARGE ";
        }
        std::cerr << tallyPosMuons     << "m+ " << tallyNegMuons     << "m- "
                  << tallyPosElectrons << "e+ " << tallyNegElectrons << "e-"
                  << std::endl;
      }
      continue;
    }


    // Fill histograms
    if (DEBUG) std::cerr << "Filling histograms for event #" << evt << "..." << std::endl;

    Int_t jetMultiplicity = 0;
    for (std::set<int>::iterator itJet=goodJets.begin(); itJet != goodJets.end(); ++itJet) {
      Jet *jet = (Jet*) brJet->At(*itJet);
      if (jet->PT > JET_MULTIPLICITY_PT_MIN) {
        jetMultiplicity++;
        h_JetPT_emu->Fill(jet->PT);
        h_JetEta_emu->Fill(jet->Eta);
        h_JetPhi_emu->Fill(jet->Phi);
      }
    }
    */

    for (unsigned int j=0; j< matched_elecs.size(); j++) {
      Electron* elec = (Electron*)brElec->At(matched_elecs[j]);
      h_ElecPT_emu->Fill(elec->PT);
      h_ElecEta_emu->Fill(elec->Eta);
      h_ElecPhi_emu->Fill(elec->Phi);
    }

    for (unsigned int j=0; j< matched_muons.size(); j++) {
      Muon *muon = (Muon*) brMuon->At(matched_muons[j]);
      h_MuonPT_emu->Fill(muon->PT);
      h_MuonEta_emu->Fill(muon->Eta);
      h_MuonPhi_emu->Fill(muon->Phi);
    }

    h_truth_MET_SumAll_emu->Fill(truth_met_all); 
    h_truth_MET_SumNeutrals_emu->Fill(truth_met_neutrals); 


    MET = ((MissingET*) brMET->At(0))->MET;
    h_MET_emu->Fill(MET);

    /*
    h_nJets_emu->Fill(jetMultiplicity);
    if (MET > 200) MET = 200-1e-5;
    if (jetMultiplicity > 1) jetMultiplicity = 1;
    h_METNjets->Fill(MET, jetMultiplicity);
    */
    tallyFilledEvents++;
  }

  Long64_t effecEvents = h_METNjets->GetEffectiveEntries();
  //if (DEBUG) std::cerr << "Rescaling by " << norm << "/" << effecEvents << std::endl;

  //h_METNjets->Scale(norm/effecEvents);
  WriteHistograms();

  std::cerr << "Number of valid dilepton events: " << tallyFilledEvents;
  std::cerr << " from " << nEvents << std::endl;

  CleanHistograms();
  delete the_histos;

  return 0;
}

Double_t calcDeltaR(Double_t phi1, Double_t phi2, Double_t eta1, Double_t eta2) 
{
  Double_t deltaPhi = abs(phi1 - phi2);
  Double_t DeltaR;
  if (deltaPhi > TMath::Pi()) deltaPhi = 2.0*TMath::Pi() - deltaPhi;

  return DeltaR = TMath::Sqrt(TMath::Power(eta1 - eta2, 2) + TMath::Power(deltaPhi, 2));
}
