#ifndef EXTRACT_AIDA_TEMPLATE_H
#define EXTRACT_AIDA_TEMPLATE_H

#include <iostream>
#include <math.h>

#include "TMath.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#include "classes/DelphesClasses.h"

// Forward declaration of Delphes classes
//class Electron;
//class Muon;
//class Jet;
//class GenParticle;

void extractAidaTemplate(
  std::string inputFile,
  std::string outputFile,
  Float_t normalisation,
  std::string histoName
  );

Double_t calcDeltaR(
  Double_t phi1,
  Double_t phi2,
  Double_t eta1,
  Double_t eta2
  );

TH2D *h_METNjets;
TH1F *h_nJets;
TH1F *h_nJets_emu;
TH1F *h_nJets_overlap;
TH1F *h_MET;
TH1F *h_MET_emu;
TH1F *h_truth_METPhi_SumAll;
TH1F *h_truth_METPhi_SumNeutrals; 
TH1F *h_truth_MET_SumAll;
TH1F *h_truth_MET_SumNeutrals; 
TH1F *h_truth_MET_SumAll_emu;
TH1F *h_truth_MET_SumNeutrals_emu; 

TH1F *h_JetPT;
TH1F *h_ElecPT;
TH1F *h_MuonPT;
TH1F *h_JetEta;
TH1F *h_ElecEta;
TH1F *h_MuonEta;
TH1F *h_JetPhi;
TH1F *h_ElecPhi;
TH1F *h_MuonPhi;
TH1F *h_JetPT_emu;
TH1F *h_ElecPT_emu;
TH1F *h_MuonPT_emu;
TH1F *h_JetEta_emu;
TH1F *h_ElecEta_emu;
TH1F *h_MuonEta_emu;
TH1F *h_JetPhi_emu;
TH1F *h_ElecPhi_emu;
TH1F *h_MuonPhi_emu;

TH1F *h_decay;
TH1F *h_decayfid;
TH1F *h_decaymatched;

TH1F *h_truth_nElec; 
TH1F *h_truth_nMuon; 
TH1F *h_truth_nElecFid; 
TH1F *h_truth_nMuonFid; 

TH1F *h_truth_ElecPT; 
TH1F *h_truth_MuonPT;
TH1F *h_truth_MuonEta;
TH1F *h_truth_ElecEta;
TH1F *h_truth_MuonPhi;
TH1F *h_truth_ElecPhi;

TH1F *h_truth_ElecPTFid; 
TH1F *h_truth_MuonPTFid;
TH1F *h_truth_MuonEtaFid;
TH1F *h_truth_ElecEtaFid;
TH1F *h_truth_MuonPhiFid;
TH1F *h_truth_ElecPhiFid;

TH1F *h_lepclass_truth;
TH1F *h_lepclass_fid;
TH1F *h_lepclass_matched;
TH1F *h_MElecPT; 
TH1F *h_MMuonPT;
TH1F *h_MMuonEta;
TH1F *h_MElecEta;
TH1F *h_MMuonPhi;
TH1F *h_MElecPhi;
TH1F *h_nMMuon; 
TH1F *h_nMElec;


TH1F *h_MElecPtRes;

TObjArray *the_histos;


void CreateHistograms(std::string histo_name);
void CleanHistograms();
void WriteHistograms();

bool isGoodElectron(Electron *electron, bool explain=false);
bool isGoodMuon(Muon *muon, bool explain=false);
bool isGoodJet(Jet *jet, bool explain=false);


bool isGoodTruthElectron(float pt, float eta, bool explain);
bool isGoodTruthMuon(float pt, float eta, bool explain);


int Which_Event_Category(int e_plus, int e_minus, int mu_plus, int mu_minus);
int Event_CategoryTruth(TClonesArray *brTruth, std::vector<int> electrons, std::vector<int> muons);
int Event_CategoryReco(TClonesArray *brElectron,TClonesArray *brMuon, std::vector<int> electrons, std::vector<int> muons);


void MET_Truth_SumAll(TClonesArray *brTruth, float &truth_met, float &truth_met_phi);
void MET_Truth_SumNeutrals(TClonesArray *brTruth, float &truth_met, float &truth_met_phi);

int FindLastParticle(TClonesArray *brTruth, int ipos, int PID);
void DecayChannel(TClonesArray *brTruth);
void PrintTruthParticle(std::string prefix,GenParticle *part,int index);
void PrintList_TruthParticles(std::string prefix,TClonesArray *brTruth,std::vector<int> leptons, std::vector<int> bosons);
void Fill_TruthLeptons_Histograms(TClonesArray *brTruth, std::vector<int> truth_electrons,std::vector<int> truth_muons );
void Fill_TruthLeptonsFid_Histograms(TClonesArray *brTruth, std::vector<int> truth_electrons,std::vector<int> truth_muons );
void Fill_MatchedLeptons_Histograms(TClonesArray *brElec,TClonesArray *brMuon, std::vector<int> matched_elecs, std::vector<int> matched_muons);

bool isLeptonPrompt(TClonesArray *brTruth, int ipos,std::vector<int> &boson_i);
int Truth_Prompt_Lepton(TClonesArray *brTruth, int PID,std::vector<int> &leptons,std::vector<int> &boson_i);
void Truth_Prompt_LeptonFid(TClonesArray *brTruth, std::vector<int> prompt_e_truth, std::vector<int> prompt_mu_truth,
                             std::vector<int> &fid_e_truth,std::vector<int> &fid_mu_truth,bool explain);



#endif
