
#include <vector>
#include <algorithm>
#include <iostream>
// Pythia8
#include "Pythia8/Pythia.h"

// ROOT
#include "TROOT.h"
#include "TTree.h"
#include "TH1.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"

template class std::vector<float>;

//#include "Particle.hpp"
#include "FastSim.hpp"


//////////////////////////////////////////////////


class MyAnalysis {

public:

  // Constructor can be empty.
  MyAnalysis() {}

  ~MyAnalysis();

  // Initialization actions.
  void init(std::string fastsim_init, std::string output_rootfilename, int debug_level);

  // clears the vectors, should be used after each event
  void clear();
  void clear_vectors(); // clear the vectors associated with the TTree
  // Analysis of each new event.
  void analyze(Pythia8::Event& event);
  void FillTree(hep_simple_lib::Event &fastevent,vector<hep_simple_lib::Particle*> truth_el, 
                                                vector<hep_simple_lib::Particle*> truth_mu,
                                                vector<hep_simple_lib::Particle*> truth_ph,
                                                vector<hep_simple_lib::Particle*> truth_bjets,
                                                vector<hep_simple_lib::Particle*> truth_htau);
  void FillHistograms_Z(hep_simple_lib::Event &fastevent);
  void FillHistograms_W(hep_simple_lib::Event &fastevent);
  // Show final results.
  void finish();

  TFile *m_ROOToutFile;

  void SelectParticles(Pythia8::Event& event);

  void PrintParticles( vector<hep_simple_lib::Particle*> list);

private:

  // Declare variables and objects that span init - analyze - finish.
  int  nEvt;

  std::string _outroot_filename;
  std::string _fastsim_initfilename;
  fast_sim::FastSim _sim;

  // the list of particles that are input to the detector response
  vector<hep_simple_lib::Particle*> _electrons;
  vector<hep_simple_lib::Particle*> _muons;
  vector<hep_simple_lib::Particle*> _photons;
  vector<hep_simple_lib::Particle*> _bjets;
  vector<hep_simple_lib::Particle*> _tauhads;
  vector<hep_simple_lib::Particle*> _chargedhads;
  vector<hep_simple_lib::Particle*> _nonprompt_leptons;
  vector<hep_simple_lib::Particle*> _weakly_interacting; // stdm neutrinos, susy neutralinos

  
  TH1F *_hBosonPt, *_hBosoneta, *_hBosonphi, *_hBosonMass_truth,*_hBosonMass_reco;
  TH1F *_hlep1Pt, *_hlep1eta, *_hlep1phi, *_hlep1iso; 
  TH1F *_hlep2Pt, *_hlep2eta, *_hlep2phi, *_hlep2iso; 
  TH1F *_hlep1Pt_truth, *_hlep1eta_truth, *_hlep1phi_truth; 
  TH1F *_hlep2Pt_truth, *_hlep2eta_truth, *_hlep2phi_truth; 
  TH1F *_hNelec,*_hNelec_truth;
  TH1F *_hNjet,*_hjetpT, *_hjeteta, *_hjetphi;   
  TH1F *_hinv, *_hmet;
  TH1F *_hinv_truth, *_hmet_truth;

  TTree *_tree;

  // ntuple variables
  // electrons
  int _nelecs;
  std::vector<int> _el_charge;
  std::vector<float> _el_eta;
  std::vector<float> _el_phi;
  std::vector<float> _el_pt;
  std::vector<float> _el_px;
  std::vector<float> _el_py;
  std::vector<float> _el_pz;
  std::vector<float> _el_iso;
  std::vector<int>   _el_type;
  // truth electrons
  int _truth_nelecs;
  std::vector<int> _truth_el_charge;
  std::vector<float> _truth_el_eta;
  std::vector<float> _truth_el_phi;
  std::vector<float> _truth_el_pt;
  std::vector<float> _truth_el_px;
  std::vector<float> _truth_el_py;
  std::vector<float> _truth_el_pz;
  std::vector<float> _truth_el_iso;
  int _nphotons;
  std::vector<int> _ph_charge;
  std::vector<float> _ph_eta;
  std::vector<float> _ph_phi;
  std::vector<float> _ph_pt;
  std::vector<float> _ph_px;
  std::vector<float> _ph_py;
  std::vector<float> _ph_pz;
  std::vector<float> _ph_iso;
  std::vector<int>   _ph_type;
  // truth _phectrons
  int _truth_nphotons;
  std::vector<int> _truth_ph_charge;
  std::vector<float> _truth_ph_eta;
  std::vector<float> _truth_ph_phi;
  std::vector<float> _truth_ph_pt;
  std::vector<float> _truth_ph_px;
  std::vector<float> _truth_ph_py;
  std::vector<float> _truth_ph_pz;
  std::vector<float> _truth_ph_iso;
  // muons
  int _nmuons;
  std::vector<int> _mu_charge;
  std::vector<float> _mu_eta;
  std::vector<float> _mu_phi;
  std::vector<float> _mu_pt;
  std::vector<float> _mu_px;
  std::vector<float> _mu_py;
  std::vector<float> _mu_pz;
  std::vector<float> _mu_iso;
  std::vector<int>   _mu_type;
  int _truth_nmuons;
  std::vector<int> _truth_mu_charge;
  std::vector<float> _truth_mu_eta;
  std::vector<float> _truth_mu_phi;
  std::vector<float> _truth_mu_pt;
  std::vector<float> _truth_mu_px;
  std::vector<float> _truth_mu_py;
  std::vector<float> _truth_mu_pz;
  std::vector<float> _truth_mu_iso;
   // had taus
  int _nhtaus;
  std::vector<int> _htau_charge;
  std::vector<float> _htau_eta;
  std::vector<float> _htau_phi;
  std::vector<float> _htau_pt;
  std::vector<float> _htau_px;
  std::vector<float> _htau_py;
  std::vector<float> _htau_pz;
  std::vector<float> _htau_iso;
  std::vector<int>   _htau_type;
  int _truth_nhtaus;
  std::vector<int> _truth_htau_charge;
  std::vector<float> _truth_htau_eta;
  std::vector<float> _truth_htau_phi;
  std::vector<float> _truth_htau_pt;
  std::vector<float> _truth_htau_px;
  std::vector<float> _truth_htau_py;
  std::vector<float> _truth_htau_pz;
  std::vector<float> _truth_htau_iso;
 
  // jets
  int _njets;
  std::vector<float> _jets_pt;
  std::vector<float> _jets_eta;
  std::vector<float> _jets_phi;
  std::vector<int>   _jets_btag;
  std::vector<float>   _jets_mass;
  // met
  float _met_val;
  float _met_phi;
  float _truth_met_val;
  float _truth_met_phi;
  float _truth_sumet;
  // sumet
  float _sumet;

  // _truth important particles
  int _truth_nparticles;
  std::vector<float> _truth_particles_pt;
  std::vector<float> _truth_particles_eta;
  std::vector<float> _truth_particles_phi;
  std::vector<int>    _truth_particles_charge;
  std::vector<float> _truth_particles_mass;
  std::vector<int> _truth_particles_pid;
 

};

//--------------------------------------------------------------------------

// The initialization code.

MyAnalysis::~MyAnalysis() {

  clear();

  delete m_ROOToutFile;
}

void MyAnalysis::clear() {

  
  // it will clear the vectors at the end of each event
  #define DELETE_PTRVEC(vec) for (size_t i = 0; i < vec.size();i++) delete vec[i]; vec.clear()
 
  DELETE_PTRVEC(_electrons);
  DELETE_PTRVEC(_muons);
  DELETE_PTRVEC(_photons);
  DELETE_PTRVEC(_bjets);
  DELETE_PTRVEC(_tauhads);
  DELETE_PTRVEC(_chargedhads);
  DELETE_PTRVEC(_weakly_interacting);
  DELETE_PTRVEC(_nonprompt_leptons);

  #undef DELETE_PTRVEC
  
}

void MyAnalysis::clear_vectors() {

  _nelecs = 0;
  _el_charge.clear();
  _el_eta.clear();

  _el_phi.clear();
  _el_pt.clear();
  _el_px.clear();
  _el_py.clear();
  _el_pz.clear();
  _el_iso.clear();
  _el_type.clear();
  // truth electrons
  _truth_nelecs = 0;
  _truth_el_charge.clear();
  _truth_el_eta.clear();
  _truth_el_phi.clear();
  _truth_el_pt.clear();
  _truth_el_px.clear();
  _truth_el_py.clear();
  _truth_el_pz.clear();
  _truth_el_iso.clear();
  _nphotons = 0;
  _ph_charge.clear();
  _ph_eta.clear();
  _ph_phi.clear();
  _ph_pt.clear();
  _ph_px.clear();
  _ph_py.clear();
  _ph_pz.clear();
  _ph_iso.clear();
  _ph_type.clear();
  // truth _phectrons
  _truth_nphotons = 0;
  _truth_ph_charge.clear();
  _truth_ph_eta.clear();
  _truth_ph_phi.clear();
  _truth_ph_pt.clear();
  _truth_ph_px.clear();
  _truth_ph_py.clear();
  _truth_ph_pz.clear();
  _truth_ph_iso.clear();
  // muons
  _nmuons = 0;
  _mu_charge.clear();
  _mu_eta.clear();
  _mu_phi.clear();
  _mu_pt.clear();
  _mu_px.clear();
  _mu_py.clear();
  _mu_pz.clear();
  _mu_iso.clear();
  _mu_type.clear();
  _truth_nmuons = 0;
  _truth_mu_charge.clear();
  _truth_mu_eta.clear();
  _truth_mu_phi.clear();
  _truth_mu_pt.clear();
  _truth_mu_px.clear();
  _truth_mu_py.clear();
  _truth_mu_pz.clear();
  _truth_mu_iso.clear();
   // had taus
  _nhtaus = 0;
  _htau_charge.clear();
  _htau_eta.clear();
  _htau_phi.clear();
  _htau_pt.clear();
  _htau_px.clear();
  _htau_py.clear();
  _htau_pz.clear();
  _htau_iso.clear();
  _htau_type.clear();
  _truth_nhtaus = 0;
  _truth_htau_charge.clear();
  _truth_htau_eta.clear();
  _truth_htau_phi.clear();
  _truth_htau_pt.clear();
  _truth_htau_px.clear();
  _truth_htau_py.clear();
  _truth_htau_pz.clear();
  _truth_htau_iso.clear();
 
  // jets
  _njets = 0;
  _jets_pt.clear();
  _jets_eta.clear();
  _jets_phi.clear();
  _jets_btag.clear();
  _jets_mass.clear();
  // met
  _met_val = (float)0.0;
  _met_phi = (float)0.0;
  _truth_met_val = (float)0.0;
  _truth_met_phi = (float)0.0;
  // sumet
  _sumet = (float)0.0;
  _truth_sumet = (float)0.0;

  // _truth important particles
  _truth_nparticles = 0;
  _truth_particles_pt.clear();
  _truth_particles_eta.clear();
  _truth_particles_phi.clear();
  _truth_particles_charge.clear();
  _truth_particles_mass.clear();
 
}


void MyAnalysis::init(std::string fastsim_initfilename, std::string output_rootfilename, int verbosity_level) {

  // Initialize counter for number of events.
  nEvt = 0;

  clear_vectors();
  
  m_ROOToutFile= new TFile(output_rootfilename.c_str(),"RECREATE");

  _tree = new TTree("FastSim","FastSim and truth particles");

  // the integer variables
  _tree->Branch("nelecs",&_nelecs,"nelecs/I");
  _tree->Branch("truth_nelecs",&_truth_nelecs,"truth_nelecs/I");
  _tree->Branch("nmuons",&_nmuons,"nmuons/I");
  _tree->Branch("truth_nmuons",&_truth_nmuons,"truth_nmuons/I");
  _tree->Branch("nphotons",&_nphotons,"nphotons/I");
  _tree->Branch("truth_nphotons",&_truth_nphotons,"truth_nphotons/I");
  _tree->Branch("nhtaus",&_nhtaus,"nhtaus/I");
  _tree->Branch("truth_nhtaus",&_truth_nhtaus,"truth_nhtaus/I");
  _tree->Branch("njets",&_njets,"njets/I");
  _tree->Branch("truth_nparticles",&_truth_nparticles,"truth_nparticles/I");

  // met                             
  _tree->Branch("met_val",&_met_val,"met_val/F");
  _tree->Branch("met_phi",&_met_phi,"met_phi/F");
  _tree->Branch("truth_met_val",&_truth_met_val,"truth_met_val/F");
  _tree->Branch("truth_met_phi",&_truth_met_phi,"truth_met_phi/F");
  // sumet
  _tree->Branch("sumet",&_sumet,"sumet/F");
  _tree->Branch("truth_sumet",&_truth_sumet,"truth_sumet/F");




  // fastsim electrons
  //_tree->Branch("el_charge","std::vector<float>",&_el_charge); 
  _tree->Branch("el_charge",&_el_charge); 
  _tree->Branch("el_eta",&_el_eta);
  _tree->Branch("el_phi",&_el_phi);
  _tree->Branch("el_pt",&_el_pt);
  _tree->Branch("el_px",&_el_px);
  _tree->Branch("el_py",&_el_py);
  _tree->Branch("el_pz",&_el_pz);
  _tree->Branch("el_iso",&_el_iso);
  _tree->Branch("el_type",&_el_type);

  // truth electrons
  _tree->Branch("truth_el_charge",&_truth_el_charge); 
  _tree->Branch("truth_el_eta",   &_truth_el_eta);
  _tree->Branch("truth_el_phi",   &_truth_el_phi);
  _tree->Branch("truth_el_pt",    &_truth_el_pt);
  _tree->Branch("truth_el_px",    &_truth_el_px);
  _tree->Branch("truth_el_py",    &_truth_el_py);
  _tree->Branch("truth_el_pz",    &_truth_el_pz);
  _tree->Branch("truth_el_iso",   &_truth_el_iso);
                                  
 //fastsim photons
  _tree->Branch("ph_charge",&_ph_charge); 
  _tree->Branch("ph_eta",   &_ph_eta);
  _tree->Branch("ph_phi",   &_ph_phi);
  _tree->Branch("ph_pt",    &_ph_pt);
  _tree->Branch("ph_px",    &_ph_px);
  _tree->Branch("ph_py",    &_ph_py);
  _tree->Branch("ph_pz",    &_ph_pz);
  _tree->Branch("ph_iso",   &_ph_iso);
  _tree->Branch("ph_type",  &_ph_type);

  // truth photons
  _tree->Branch("truth_ph_charge",&_truth_ph_charge);
  _tree->Branch("truth_ph_eta",   &_truth_ph_eta);
  _tree->Branch("truth_ph_phi",   &_truth_ph_phi);
  _tree->Branch("truth_ph_pt",    &_truth_ph_pt);
  _tree->Branch("truth_ph_px",    &_truth_ph_px);
  _tree->Branch("truth_ph_py",    &_truth_ph_py);
  _tree->Branch("truth_ph_pz",    &_truth_ph_pz);
  _tree->Branch("truth_ph_iso",   &_truth_ph_iso);

  // fastsim muons
  _tree->Branch("mu_charge",&_mu_charge);
  _tree->Branch("mu_eta",   &_mu_eta);
  _tree->Branch("mu_phi",   &_mu_phi);
  _tree->Branch("mu_pt",    &_mu_pt);
  _tree->Branch("mu_px",    &_mu_px);
  _tree->Branch("mu_py",    &_mu_py);
  _tree->Branch("mu_pz",    &_mu_pz);
  _tree->Branch("mu_iso",   &_mu_iso);
  _tree->Branch("mu_type",  &_mu_type);

  // truth prompt muons
  _tree->Branch("truth_mu_charge",&_truth_mu_charge);
  _tree->Branch("truth_mu_eta",   &_truth_mu_eta);
  _tree->Branch("truth_mu_phi",   &_truth_mu_phi);
  _tree->Branch("truth_mu_pt",    &_truth_mu_pt);
  _tree->Branch("truth_mu_px",    &_truth_mu_px);
  _tree->Branch("truth_mu_py",    &_truth_mu_py);
  _tree->Branch("truth_mu_pz",    &_truth_mu_pz);
  _tree->Branch("truth_mu_iso",   &_truth_mu_iso);

  // fastsim hadronic taus
  _tree->Branch("htau_charge",&_htau_charge);
  _tree->Branch("htau_eta",   &_htau_eta);
  _tree->Branch("htau_phi",   &_htau_phi);
  _tree->Branch("htau_pt",    &_htau_pt);
  _tree->Branch("htau_px",    &_htau_px);
  _tree->Branch("htau_py",    &_htau_py);
  _tree->Branch("htau_pz",    &_htau_pz);
  _tree->Branch("htau_iso",   &_htau_iso);
  _tree->Branch("htau_type",  &_htau_type);

  // truth hadronic taus
  _tree->Branch("truth_htau_charge", &_truth_htau_charge);
  _tree->Branch("truth_htau_eta",    &_truth_htau_eta);
  _tree->Branch("truth_htau_phi",    &_truth_htau_phi);
  _tree->Branch("truth_htau_pt",     &_truth_htau_pt);
  _tree->Branch("truth_htau_px",     &_truth_htau_px);
  _tree->Branch("truth_htau_py",     &_truth_htau_py);
  _tree->Branch("truth_htau_pz",     &_truth_htau_pz);
  _tree->Branch("truth_htau_iso",    &_truth_htau_iso);
 
  // fastsim jets
  _tree->Branch("jets_pt",   &_jets_pt);
  _tree->Branch("jets_eta",  &_jets_eta);
  _tree->Branch("jets_phi",  &_jets_phi);
  _tree->Branch("jets_btag", &_jets_btag);

  // truth important particles
  _tree->Branch("truth_particles_pt",     &_truth_particles_pt);
  _tree->Branch("truth_particles_eta",    &_truth_particles_eta);
  _tree->Branch("truth_particles_phi",    &_truth_particles_phi);
  _tree->Branch("truth_particles_charge", &_truth_particles_charge);
  _tree->Branch("truth_particles_mass",   &_truth_particles_mass);
                                         

  /*
  // Book histograms.
  _hBosonPt = new TH1F("BosonPt"," Boson Generated Pt;GeV;",100, 0., 200.);
  _hBosoneta = new TH1F("Bosoneta"," Boson Generated eta;",100, -5., 5.);
  _hBosonphi = new TH1F( "BosonPhi","Boson Generated Phi;",100, -6.0, 6.0);

  _hlep1Pt_truth = new TH1F("lep1PtTruth","Leading lep Pt (Truth);GeV;",100, 0., 200.);
  _hlep1eta_truth = new TH1F("lep1etaTruth","Leading lep eta (Truth);",100, -5., 5.);
  _hlep1phi_truth = new TH1F( "lep1PhiTruth","Leading lep Phi (Truth);",100, -6.0, 6.0);

  _hlep2Pt_truth = new TH1F("lep2PtTruth","SubLeading lep Pt (Truth);GeV;",100, 0., 200.);
  _hlep2eta_truth = new TH1F("lep2etaTruth","SubLeading  lep pseudorapidity (Truth);eta",100, -5., 5.);
  _hlep2phi_truth = new TH1F( "lep2PhiTruth","SubLeading lep Phi (Truth);",100, -6.0, 6.0);

  _hlep1Pt = new TH1F("lep1Pt","Leading lep Pt;GeV;",100, 0., 200.);
  _hlep1eta = new TH1F("lep1eta","Leading lep eta;",100, -5., 5.);
  _hlep1phi = new TH1F( "lep1Phi","Leading lep Phi;",100, -6.0, 6.0);
  _hlep1iso = new TH1F( "lep1iso","Leading lep Iso;",50, 0,10);

  _hlep2Pt = new TH1F("lep2Pt","SubLeading lep Pt;GeV;",100, 0., 200.);
  _hlep2eta = new TH1F("lep2eta","SubLeading  lep eta;",100, -5., 5.);
  _hlep2phi = new TH1F( "lep2Phi","SubLeading lep Phi;",100, -6.0, 6.0);
  _hlep2iso = new TH1F( "lep2iso","SubLeading lep Iso;",50, 0, 10.0);

  _hNelec = new TH1F("Nelec","Number of Isolated Electrons;Number/Event",5,-0.5,4.5);
  _hNelec_truth = new TH1F("NelecTruth","Number of Electrons (Truth);Number/Event",5,-0.5,4.5);

  _hNjet = new TH1F("Njet","Number of Jets;Number/Event",10,-0.5,9.5);
  _hjetpT = new TH1F("jetPt","Jet Pt;GeV;",100, 0., 200.);
  _hjeteta = new TH1F("jeteta","Jet eta;",100, -5., 5.);
  _hjetphi = new TH1F( "jetPhi","Jet Phi;",100, -6.0, 6.0);

  _hmet = new TH1F( "MET","MET (Calorimeter Sum);GeV",100, 0, 200);

  _hBosonMass_truth = new TH1F( "BosonMassTruth","Z Invariant Mass (Truth);GeV",100, 0, 200);
  _hBosonMass_reco = new TH1F( "BosonMassReco","Z Invariant Mass (FastSim);GeV",100, 0, 200);
  _hmet_truth = new TH1F( "METTruth","MET (Truth);GeV",100, 0, 200);
  */


  //_sim.init(fast_sim::NOMINAL);
  _sim.init(fastsim_initfilename.c_str(),verbosity_level);
  //_sim.init(fast_sim::ACERDET);

}

// -- Histogramming

//--------------------------------------------------------------------------

void MyAnalysis::SelectParticles(Pythia8::Event& event) {
  // this method selects and categorizes the particles into the respective vectors
  //
  hep_simple_lib::Particle* chosen;

  std::cout << " event " << std::endl;
  // iterate through each of the particles, select and sort them into the different vectors
  for (int i = 0; i < event.size(); ++i) {

    //check whether the leptons are prompt first
    if (event[i].isLepton()) {

      if (event[i].isCharged()) {

          int mother1 = event[i].mother1();
          //std::cout << " mother1 " << mother1 << " mother2 " << mother2 << std::endl;
          // is the mother of the lepton a W or a Z or a tau lepton
          if (((int(fabs(event[mother1].id()))) == 23) || ((int(fabs(event[mother1].id()))) == 24) ||  ((int(fabs(event[mother1].id()))) == 15)) {

            std::cout << " prompt " << event[i].id() << " mother " << event[mother1].id() <<  ' ' << event[mother1].status() << std::endl;
            switch (int(fabs(event[i].id()))) {
              /// @todo This needs to change, it should be only prompt leptons.. not just any lepton
              case 11: // electron
                chosen = new hep_simple_lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
                _electrons.push_back(chosen);
                break;
              case 13: // muon
                chosen = new hep_simple_lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
                _muons.push_back(chosen);
                break;
              case 15: 
                {// tau - if the tau decays hadronically then we want to save the visible tau vector
                  hep_simple_lib::Particle* neu = NULL;
                  hep_simple_lib::Particle* tau = NULL;
                  int daughter1 = event[i].daughter1();
                  int daughter2 = event[i].daughter2();
                  bool leptonic = false;
                  for (int d= daughter1;d<=daughter2;d++) {
                    if (event[d].isLepton())
                      leptonic = true;
                    if ((int(fabs(event[d].id()))) == 16) 
                      neu = new hep_simple_lib::Particle(event[d].px(), event[d].py(), event[d].pz(), event[d].e(), event[d].id());
                  }
                  if ((!leptonic) && (neu != NULL)) {// we should save this as a hadronic tau decay tau_vis = tau - tau_neu

                    tau = new hep_simple_lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
                    hep_simple_lib::P4 visible_tau = tau->mom() - neu->mom();
                    chosen = new hep_simple_lib::Particle(visible_tau,event[i].id());
                    _tauhads.push_back(chosen);
                  }
                }
                break;
              default: std::cout << " what kind of prompt charged lepton are you " << std::endl;
            }
          }
          else if (event[i].isFinal()) {

            chosen = new hep_simple_lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].id());
            _nonprompt_leptons.push_back(chosen);
              std::cout << " is Final but not prompt lepton " << event[i].id() << std::endl;
          }
      }
      else { // electron neutrinos // muon neutrinos // tau neutrinos
        chosen = new hep_simple_lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].id());
        _weakly_interacting.push_back(chosen);
      }
    }
    else if (event[i].isFinal()) {

      chosen = new hep_simple_lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].id());
       switch (int(fabs(event[i].id()))) {
         case 22: // photon
           _photons.push_back(chosen);
           break;
         case 1000022: // neutralino 1
           _weakly_interacting.push_back(chosen);
           break;
          default:
           //default: // every other hadronic charged particle - for the jets or neutral pion or kaon
           _chargedhads.push_back(chosen);
       }
    }
    else {

      //if (event[i].isFinal()) {
      /*
      //       printf("neu particle %d\n",event[i].id());
      switch (int(fabs(event[i].id()))) {
      case 22: // photon

      //            P2 = event[i].px()*event[i].px() + event[i].py()*event[i].py() + event[i].pz()*event[i].pz();

      //            printf("** P2 %f E %f\n",(float)P2,(float)(event[i].e()*event[i].e()));


      chosen = new hep_simple_lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].id());
      _photons.push_back(chosen);
      break;
      case 12: // electron neutrinos
      case 14: // muon neutrinos
      case 16: // tau neutrinos
      chosen = new hep_simple_lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].id());
      _weakly_interacting.push_back(chosen);
      //            printf("neutrinos %d\n",event[i].id());

      break;
      case -2112:
      case 2112: // neutrons
      case 130: // neutral kaon - not sure what to do
      chosen = new hep_simple_lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
      _chargedhads.push_back(chosen);
      break;
      default: //printf("neutral final particle missed %d\n",event[i].id());
      }
      */
      // }
      //}

      if ((event[i].isQuark()) && (fabs(event[i].id()) == 5)) {

        int daughter1 = event[i].daughter1();
        int daughter2 = event[i].daughter2();
        bool final_b = true;
        for (int d= daughter1;d<=daughter2;d++) {
          if (fabs(event[d].id()) == 5) {
            final_b = false;
            break;
          }
        }
        if (final_b) {
          chosen = new hep_simple_lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
          _bjets.push_back(chosen);
        }
      }
    }
  }

/*
    // lets get our taus, they are unstable so we need a special list and case for them
    if ((event[i].isLepton()) && (fabs(event[i].id()) == 15)) {

      int daughter1 = event[i].daughter1();
      int daughter2 = event[i].daughter2();

      for (int j=daughter1; j<=daughter2;j++) {
        // printf("tau daughters (%d) are %d pdgId %d pt %.2f\n",j,event[j].id(),event[j].pT());
        cout << "tau daughters (" << j << ") are pdgId = " << event[j].id() << ", pt = " << event[j].pT() << " GeV" << endl;
      }

      // we need to remove the leptonically decaying taus - perhaps do an overlap removal with electrons
      chosen = new hep_simple_lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
      _tauhads.push_back(chosen);
    }
*/
}


// The event analysis code.
void MyAnalysis::analyze(Pythia8::Event& event) {

  // select the different particles for the fast simulator
  SelectParticles(event);
 
  //PrintParticles(_electrons);

//  hep_simple_lib::P4 z_p = _electrons[0]->mom()+_electrons[1]->mom();
//  _hBosonMass_truth->Fill(z_p.m());
//  cout << " The mass is " << z_p.m() << std::endl;

  //fast_sim::FastSim sim;
  //sim.init(fast_sim::ACERDET);
  //

  hep_simple_lib::Event reco_event;
  _sim.setParticles(_electrons, _muons, _photons, _nonprompt_leptons,_chargedhads, _bjets, _tauhads, _weakly_interacting);

  _sim.doDetectorResponse();
  _sim.getRecoEvent(reco_event,"medium",0.0,"medium",0.0,"medium");

  FillTree(reco_event,_electrons, _muons, _photons,_bjets, _tauhads);
  //FillHistograms_Z(reco_event);
  //FillHistograms_W(reco_event);

  //sim.printElectrons();



  // clear the vectors in the example
  clear();
  clear_vectors();
}

/*
  // Increase counter.
  ++nEvt;


  for (int i = 0; i < (int)event.size(); ++i) {
    if((event[i].isFinal())&&((abs(event[i].id())==12)||(abs(event[i].id())==14)||(abs(event[i].id())==16))){

      pxxnues += event[i].px();
      pyynues += event[i].py();

    }

  // Plot pseudorapidity distribution. Sum up charged multiplicity.
  for (int i = 0; i < event.size(); ++i) {

    switch (event[i].id()) {

      case 23: // z boson

        m_hBosonPt->Fill(event[i].pT());
        m_hBosoneta->Fill(event[i].eta());
        m_hBosonphi->Fill(event[i].phi());

        printf(" Boson %.2f (GeV) %.2f %.2f\n",event[i].pT(),event[i].eta(),event[i].phi());
        break;

      case 11: // electron
      case -11: // electron

        m_hElectronPt->Fill(event[i].pT());
        m_hElectroneta->Fill(event[i].eta());
        m_hElectronphi->Fill(event[i].phi());

        printf(" Electron %.2f (GeV) %.2f %.2f\n",event[i].pT(),event[i].eta(),event[i].phi());
        break;
      default:
        printf(" what other particles are there id %d %.2f (GeV) %.2f %.2f\n",event[i].id(),event[i].pT(),event[i].eta(),event[i].phi());

    }

  }

}
*/


void MyAnalysis::PrintParticles( vector<hep_simple_lib::Particle*> list)
{
  for (int j=0;j<(int)list.size();j++) {
    cout << "Particle " << j << " "<< list[j]->pid() << " P: "<< list[j]->mom().px() << " " << list[j]->mom().py() << " " << list[j]->mom().pz() << " " << list[j]->mom().E()
      << " Eta: " << list[j]->mom().eta() << " Phi: " << list[j]->mom().phi()  << " Pt: " << list[j]->pT() << endl;
  }

}

void MyAnalysis::FillHistograms_W(hep_simple_lib::Event &fastevent) {
  // this function fills the histograms

  _hNelec->Fill(fastevent.electrons().size());
  if ((int)fastevent.electrons().size() > 0) {

    _hmet->Fill(fastevent.met());
    _hmet_truth->Fill(fastevent.met_truth());
    _hlep1Pt->Fill(fastevent.electrons()[0]->pT());
    _hlep1eta->Fill(fastevent.electrons()[0]->mom().eta());
    _hlep1phi->Fill(fastevent.electrons()[0]->mom().phi());
    _hlep1iso->Fill(fastevent.electrons()[0]->isol());


    _hNjet->Fill((int)fastevent.jets().size());
    for (size_t i = 0; i < fastevent.jets().size();i++) {

      _hjetpT->Fill(fastevent.jets()[i]->pT());
      _hjeteta->Fill(fastevent.jets()[i]->eta());
      _hjetphi->Fill(fastevent.jets()[i]->phi());

    }
  }

}

void MyAnalysis::FillTree(hep_simple_lib::Event &fastevent,vector<hep_simple_lib::Particle*> truth_el, 
                                                           vector<hep_simple_lib::Particle*> truth_mu,
                                                           vector<hep_simple_lib::Particle*> truth_ph,
                                                           vector<hep_simple_lib::Particle*> truth_bjets,
                                                           vector<hep_simple_lib::Particle*> truth_htau
    ) {

  // lets fill the tree with the fastsim particles
  _nelecs = fastevent.electrons().size();
  _nmuons = fastevent.muons().size();
  _njets = fastevent.jets().size();

  for (size_t kk=0;kk<fastevent.electrons().size();kk++) {
    _el_pt.push_back(fastevent.electrons()[kk]->pT());
    _el_px.push_back(fastevent.electrons()[kk]->mom().px());
    _el_py.push_back(fastevent.electrons()[kk]->mom().py());
    _el_pz.push_back(fastevent.electrons()[kk]->mom().pz());
    _el_eta.push_back(fastevent.electrons()[kk]->eta());
    _el_phi.push_back(fastevent.electrons()[kk]->phi());
    _el_iso.push_back(fastevent.electrons()[kk]->isol());
    if (fastevent.electrons()[kk]->pid() > 0)
      _el_charge.push_back(1);
    else
      _el_charge.push_back(-1);
  }

  for (size_t kk=0;kk<fastevent.muons().size();kk++) {
    _mu_pt.push_back(fastevent.muons()[kk]->pT());
    _mu_px.push_back(fastevent.muons()[kk]->mom().px());
    _mu_py.push_back(fastevent.muons()[kk]->mom().py());
    _mu_pz.push_back(fastevent.muons()[kk]->mom().pz());
    _mu_eta.push_back(fastevent.muons()[kk]->eta());
    _mu_phi.push_back(fastevent.muons()[kk]->phi());
    _mu_iso.push_back(fastevent.muons()[kk]->isol());
    if (fastevent.muons()[kk]->pid() > 0)
      _mu_charge.push_back(1);
    else
      _mu_charge.push_back(-1);
 
  }

  for (size_t kk=0;kk<fastevent.jets().size();kk++) {
    _jets_pt.push_back(fastevent.jets()[kk]->pT());
    _jets_eta.push_back(fastevent.jets()[kk]->eta());
    _jets_phi.push_back(fastevent.jets()[kk]->phi());
    _jets_btag.push_back(fastevent.jets()[kk]->isBJet());
  }

  // fill the truth particles
  _truth_nelecs = truth_el.size();
  _truth_nmuons = truth_mu.size();
  _truth_nphotons = truth_ph.size();
  _truth_nhtaus = truth_htau.size();

  for (size_t kk=0;kk<truth_el.size();kk++) {
    _truth_el_pt.push_back(truth_el[kk]->pT());
    _truth_el_px.push_back(truth_el[kk]->mom().px());
    _truth_el_py.push_back(truth_el[kk]->mom().py());
    _truth_el_pz.push_back(truth_el[kk]->mom().pz());
    _truth_el_eta.push_back(truth_el[kk]->eta());
    _truth_el_phi.push_back(truth_el[kk]->phi());
    _truth_el_iso.push_back(truth_el[kk]->isol());
    if (truth_el[kk]->pid() > 0)
      _truth_el_charge.push_back(1);
    else
      _truth_el_charge.push_back(-1);
  }

  // the EtMiss
  _met_val = fastevent.met();
  _truth_met_val = fastevent.met_truth();
  _met_phi = fastevent.missingmom().phi();
  _truth_met_phi = fastevent.missingmom_truth().phi();


  _tree->Fill();

  // clear the vectors for the next event
  clear_vectors();


}

void MyAnalysis::FillHistograms_Z(hep_simple_lib::Event &fastevent) {
  // this function fills the histograms

  _hNelec->Fill(fastevent.electrons().size());



  //cout << " the size of electrons " << fastevent.electrons().size() << endl;
  if ((int)fastevent.electrons().size() == 2) {
    hep_simple_lib::P4 temp;
    temp = fastevent.electrons()[0]->mom() + fastevent.electrons()[1]->mom();
    _hBosonMass_reco->Fill(temp.m());

    _hmet->Fill(fastevent.met());
    _hmet_truth->Fill(fastevent.met_truth());

    if (fastevent.electrons()[0]->pT() >  fastevent.electrons()[1]->pT()) {
      _hlep1Pt->Fill(fastevent.electrons()[0]->pT());
      _hlep1eta->Fill(fastevent.electrons()[0]->mom().eta());
      _hlep1phi->Fill(fastevent.electrons()[0]->mom().phi());
      _hlep1iso->Fill(fastevent.electrons()[0]->isol());
      _hlep2Pt->Fill(fastevent.electrons()[1]->pT());
      _hlep2eta->Fill(fastevent.electrons()[1]->mom().eta());
      _hlep2phi->Fill(fastevent.electrons()[1]->mom().phi());
      _hlep2iso->Fill(fastevent.electrons()[1]->isol());
    }
    else {
      _hlep1Pt->Fill(fastevent.electrons()[1]->pT());
      _hlep1eta->Fill(fastevent.electrons()[1]->mom().eta());
      _hlep1phi->Fill(fastevent.electrons()[1]->mom().phi());
      _hlep2iso->Fill(fastevent.electrons()[1]->isol());
      _hlep2Pt->Fill(fastevent.electrons()[0]->pT());
      _hlep2eta->Fill(fastevent.electrons()[0]->mom().eta());
      _hlep2phi->Fill(fastevent.electrons()[0]->mom().phi());
      _hlep1iso->Fill(fastevent.electrons()[0]->isol());
    }
  }

  _hNjet->Fill((int)fastevent.jets().size());
  for (size_t i = 0; i < fastevent.jets().size();i++) {

    _hjetpT->Fill(fastevent.jets()[i]->pT());
    _hjeteta->Fill(fastevent.jets()[i]->eta());
    _hjetphi->Fill(fastevent.jets()[i]->phi());

  }


  if ((int)fastevent.muons().size() == 2) {
    hep_simple_lib::P4 temp;
    temp = fastevent.muons()[0]->mom() + fastevent.muons()[1]->mom();
    _hBosonMass_reco->Fill(temp.m());

    _hmet->Fill(fastevent.met());
    _hmet_truth->Fill(fastevent.met_truth());

    if (fastevent.muons()[0]->pT() >  fastevent.muons()[1]->pT()) {
      _hlep1Pt->Fill(fastevent.muons()[0]->pT());
      _hlep1eta->Fill(fastevent.muons()[0]->mom().eta());
      _hlep1phi->Fill(fastevent.muons()[0]->mom().phi());
      _hlep1iso->Fill(fastevent.muons()[0]->isol());
      _hlep2Pt->Fill(fastevent.muons()[1]->pT());
      _hlep2eta->Fill(fastevent.muons()[1]->mom().eta());
      _hlep2phi->Fill(fastevent.muons()[1]->mom().phi());
      _hlep2iso->Fill(fastevent.muons()[1]->isol());
    }
    else {
      _hlep1Pt->Fill(fastevent.muons()[1]->pT());
      _hlep1eta->Fill(fastevent.muons()[1]->mom().eta());
      _hlep1phi->Fill(fastevent.muons()[1]->mom().phi());
      _hlep2iso->Fill(fastevent.muons()[1]->isol());
      _hlep2Pt->Fill(fastevent.muons()[0]->pT());
      _hlep2eta->Fill(fastevent.muons()[0]->mom().eta());
      _hlep2phi->Fill(fastevent.muons()[0]->mom().phi());
      _hlep1iso->Fill(fastevent.muons()[0]->isol());
    }

  }

}




//--------------------------------------------------------------------------

// The finishing code.

void MyAnalysis::finish() {

  // Normalize histograms.
//  double binFactor = 5. / nEvt;
//  yH     *= binFactor;
//  etaChg *= binFactor;

  // Print histograms.
  //  cout << m_hBosonPt << mBosoneta << Bosonphi;
 
  cout << " finishing my analysis " << endl;
  m_ROOToutFile->Write();


}

std::string getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {

      std::string val((*itr),strlen((*itr)));
 
        return val;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

int cmdReadOptions(int argc, char* argv[], std::string &cmd_filename,  std::string &json_filename,  std::string &out_filename, int &debug_level) 
{

  if ((cmdOptionExists(argv, argv+argc, "-h")) || (argc == 0))
  {
    std::cout << " usage ./example -cmd pythia_file -j fastsim_description.json -o output_rootfilename < -d debug level > \n"
              << " -h to display this message \n"
              << " the debug level is 0 so if the option is omitted, the program will not print any debug statements\n"
              << " possible values can be 1,2, or 3 " << std::endl;

    return 1;
    
  }

  if (cmdOptionExists(argv, argv+argc, "-cmd")) 
     cmd_filename = getCmdOption(argv, argv + argc, "-cmd");
  else {
     std::cout << " usage ./example -cmd pythia_file -j fastsim_description.json -o output_rootfilename < -d debug level > \n"
              << " the debug level is 0 so if the option is omitted, the program will not print any debug statements\n"
              << " possible values can be 1,2, or 3 " << std::endl;
    return 1;
    
  }
  std::cout << " command file " << cmd_filename << std::endl;

  if (cmdOptionExists(argv, argv+argc, "-j")) 
     json_filename = getCmdOption(argv, argv + argc, "-j");
  else {
     std::cout << " usage ./example -cmd pythia_file -j fastsim_description.json -o output_rootfilename < -d debug level > \n"
              << " the debug level is 0 so if the option is omitted, the program will not print any debug statements\n"
              << " possible values can be 1,2, or 3 " << std::endl;
    return 1;
    
  }
  std::cout << " json file " << json_filename << std::endl;

  if (cmdOptionExists(argv, argv+argc, "-o")) 
     out_filename = getCmdOption(argv, argv + argc, "-o");
  else {
     std::cout << " usage ./example -cmd pythia_file -j fastsim_description.json -o output_rootfilename < -d debug level > \n"
              << " the debug level is 0 so if the option is omitted, the program will not print any debug statements\n"
              << " possible values can be 1,2, or 3 " << std::endl;
    return 1;
    
  }


  std::cout << " output file " << out_filename << std::endl;

  if (cmdOptionExists(argv, argv+argc, "-d")) {
    std::string val = getCmdOption(argv, argv + argc, "-d");
     debug_level = stol(val.c_str());
  }


  return 0;

}


//==========================================================================

// You should not need to touch the main program: its actions are
// determined by the .cmnd file and the rest belongs in MyAnalysis.

int main(int argc, char* argv[]) {

  std::string cmd_filename;
  std::string json_filename;
  std::string out_filename;
  int debug_level=0;

  if (cmdReadOptions(argc,argv,cmd_filename,json_filename,out_filename,debug_level) != 0)
    return 1;

  // Check that the provided file name corresponds to an existing file.
  ifstream is1(cmd_filename);
  if (!is1) {
    std::cout << " Command-line file " << cmd_filename << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  ifstream is2(json_filename);
  if (!is2) {
    std::cout << " fastsim init file " << json_filename << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // root dictionary loading
  gROOT->ProcessLine("#include <vector>");
  // Enable autoloading
//  gInterpreter->EnableAutoLoading();
  //ifstream is3(argv[3]);
  //if (!is3) {
  //  cerr << " root file " << argv[3] << " was not found. \n"
  //       << " Program stopped! " << endl;
  //  return 1;
  //}

  // Confirm to the user what the inputs are
  cout << " PYTHIA settings will be read from file " << cmd_filename << endl;
  cout << " FastSim settings will be read from file " << json_filename << endl;
  cout << " ROOT histograms will be saved to file " << out_filename << endl;


  // Declare generator. Read in commands from external file.
  Pythia8::Pythia pythia;
  pythia.readFile(cmd_filename);


  // Initialization.
  pythia.init();

  // Declare user analysis class. Do initialization part of it.
  MyAnalysis myAnalysis;
  myAnalysis.init(json_filename,out_filename,debug_level);

  // Read in number of event and maximal number of aborts.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    //printf("Event %d\n",iEvent);

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // User Analysis of current event.
    myAnalysis.analyze( pythia.event);

  // End of event loop.
  }

  // Final statistics.
  pythia.stat();

  // User finishing.
  myAnalysis.finish();

  // Done.
  return 0;
}
