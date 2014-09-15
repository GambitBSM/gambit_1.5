#include "omp.h"

// Pythia8
#include "Pythia8/Pythia.h"

// ROOT
#include "TH1.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"

#include "Particle.hpp"
#include "FastSim.hpp"

#include <iostream>

//////////////////////////////////////////////////


class MyAnalysis {

public:

  // Constructor can be empty.
  MyAnalysis() {}

  ~MyAnalysis();

  // Initialization actions.
  void init();

  // clears the vectors, should be used after each event
  void clear();
  // Analysis of each new event.
  void analyze(Pythia8::Event& event);

  void FillHistograms(HEP_Simple_Lib::Event &fastevent);
  // Show final results.
  void finish();

  TFile *m_ROOToutFile;

  void SelectParticles(Pythia8::Event& event);

  void PrintParticles( vector<HEP_Simple_Lib::Particle*> list);

private:

  // Declare variables and objects that span init - analyze - finish.
  int  nEvt;

  //fast_sim::FastSim _sim;

  // the list of particles that are input to the detector response
  vector<HEP_Simple_Lib::Particle*> _electrons;
  vector<HEP_Simple_Lib::Particle*> _muons;
  vector<HEP_Simple_Lib::Particle*> _photons;
  vector<HEP_Simple_Lib::Particle*> _bjets;
  vector<HEP_Simple_Lib::Particle*> _tauhads;
  vector<HEP_Simple_Lib::Particle*> _chargedhads;
  vector<HEP_Simple_Lib::Particle*> _weakly_interacting; // stdm neutrinos, susy neutralinos

  
  TH1F *_hBosonPt, *_hBosoneta, *_hBosonphi, *_hBosonMass_truth,*_hBosonMass_reco;
  TH1F *_hlep1Pt, *_hlep1eta, *_hlep1phi, *_hlep1iso; 
  TH1F *_hlep2Pt, *_hlep2eta, *_hlep2phi, *_hlep2iso; 
  TH1F *_hlep1Pt_truth, *_hlep1eta_truth, *_hlep1phi_truth; 
  TH1F *_hlep2Pt_truth, *_hlep2eta_truth, *_hlep2phi_truth; 
  TH1F *_hNelec,*_hNelec_truth;
  TH1F *_hNjet,*_hjetpT, *_hjeteta, *_hjetphi;   
  TH1F *_hinv, *_hmet;
  TH1F *_hinv_truth, *_hmet_truth;
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

  #undef DELETE_PTRVEC
  
}


void MyAnalysis::init() {

  // Initialize counter for number of events.
  nEvt = 0;


  m_ROOToutFile= new TFile("FastSim_Pythia8_output.root","RECREATE");

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

  _hmet = new TH1F( "MET","MET;GeV",100, 0, 200);

  _hBosonMass_truth = new TH1F( "BosonMassTruth","Z Invariant Mass (Truth);GeV",100, 0, 200);
  _hBosonMass_reco = new TH1F( "BosonMassReco","Z Invariant Mass (FastSim);GeV",100, 0, 200);
  _hmet_truth = new TH1F( "METTruth","MET (Truth);GeV",100, 0, 200);


  //_sim.init(fast_sim::NOMINAL);
  //_sim.init(fast_sim::ACERDET);

}

// -- Histogramming

//--------------------------------------------------------------------------

void MyAnalysis::SelectParticles(Pythia8::Event& event) {
  // this method selects and categorizes the particles into the respective vectors
  //
  HEP_Simple_Lib::Particle* chosen;

  // iterate through each of the particles, select and sort them into the different vectors
  for (int i = 0; i < event.size(); ++i) {

    if (event[i].isFinal()) {

      if (event[i].isCharged()) {
        switch (int(fabs(event[i].id()))) {
          /// @todo This needs to change, it should be only prompt leptons.. not just any lepton
          case 11: // electron
            chosen = new HEP_Simple_Lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
            _electrons.push_back(chosen);
            break;
          case 13: // muon
            chosen = new HEP_Simple_Lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
            _muons.push_back(chosen);
            break;
          default: // every other hadronic charged particle - for the jets
            chosen = new HEP_Simple_Lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
            _chargedhads.push_back(chosen);
//           printf("charged final particle missed %d\n",event[i].id());
        }
      }
      else {

 //       printf("neu particle %d\n",event[i].id());
        switch (int(fabs(event[i].id()))) {
          case 22: // photon
            
//            P2 = event[i].px()*event[i].px() + event[i].py()*event[i].py() + event[i].pz()*event[i].pz();

//            printf("** P2 %f E %f\n",(float)P2,(float)(event[i].e()*event[i].e()));


            chosen = new HEP_Simple_Lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].id());
            _photons.push_back(chosen);
            break;
          case 12: // electron neutrinos
          case 14: // muon neutrinos
          case 16: // tau neutrinos
            chosen = new HEP_Simple_Lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].id());
            _weakly_interacting.push_back(chosen);
//            printf("neutrinos %d\n",event[i].id());

            break;
          case -2112:
          case 2112: // neutrons
          case 130: // neutral kaon - not sure what to do
            chosen = new HEP_Simple_Lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
            _chargedhads.push_back(chosen);
            break;
          default: printf("neutral final particle missed %d\n",event[i].id());
        }
      }
    }

    if ((event[i].isQuark()) && (fabs(event[i].id()) == 6)) {

      chosen = new HEP_Simple_Lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
      _bjets.push_back(chosen);
    }


    // lets get our taus, they are unstable so we need a special list and case for them
    if ((event[i].isLepton()) && (fabs(event[i].id()) == 15)) {

      int daughter1 = event[i].daughter1();
      int daughter2 = event[i].daughter2();

      for (int j=daughter1; j<=daughter2;j++) {
        // printf("tau daughters (%d) are %d pdgId %d pt %.2f\n",j,event[j].id(),event[j].pT());
        std::cout << "tau daughters (" << j << ") are pdgId = " << event[j].id() << ", pt = " << event[j].pT() << " GeV" << endl;
      }

      // we need to remove the leptonically decaying taus - perhaps do an overlap removal with electrons
      chosen = new HEP_Simple_Lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
      _tauhads.push_back(chosen);
    }


  }
}


// The event analysis code.
void MyAnalysis::analyze(Pythia8::Event& event) {

  // select the different particles for the fast simulator
  SelectParticles(event);
 
  //PrintParticles(_electrons);

  HEP_Simple_Lib::P4 z_p = _electrons[0]->mom()+_electrons[1]->mom();
  _hBosonMass_truth->Fill(z_p.m());
  //cout << " The mass is " << z_p.m() << std::endl;

  /*
  //fast_sim::FastSim sim;
  //sim.init(fast_sim::ACERDET);
  //

  HEP_Simple_Lib::Event reco_event;
  _sim.setParticles(_electrons, _muons, _photons, _chargedhads, _bjets, _tauhads, _weakly_interacting);

  //cout << " do detector response " << std::endl;
  _sim.doDetectorResponse();
  _sim.getRecoEvent(reco_event);
  FillHistograms(reco_event);

  //sim.printElectrons();

  */

  // clear the vectors in the example
  clear();
}

/*
  // Increase counter.
  ++nEvt;


  for (int i = 0; i < (int)event.size(); ++i) {
    if((event[i].isFinal())&&((abs(event[i].id())==12)||(abs(event[i].id())==14)||(abs(event[i].id())==16))){

      pxxnues += event[i].px();
      pyynues += event[i].py();

    }

  // Plot pseudorapidity distribution. Sum upp charged multiplicity.
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


void MyAnalysis::PrintParticles( vector<HEP_Simple_Lib::Particle*> list)
{
  for (int j=0;j<(int)list.size();j++) {
    std::cout << "Particle " << j << " "<< list[j]->pid() << " P: "<< list[j]->mom().px() << " " << list[j]->mom().py() << " " << list[j]->mom().pz() << " " << list[j]->mom().E()
      << " Eta: " << list[j]->mom().eta() << " Phi: " << list[j]->mom().phi()  << " Pt: " << list[j]->pT() << endl;
  }

}

void MyAnalysis::FillHistograms(HEP_Simple_Lib::Event &fastevent) {
  // this function fills the histograms

  _hNelec->Fill(fastevent.electrons().size());



  //cout << " the size of electrons " << fastevent.electrons().size() << endl;
  if ((int)fastevent.electrons().size() == 2) {
    HEP_Simple_Lib::P4 temp;
    temp = fastevent.electrons()[0]->mom() + fastevent.electrons()[1]->mom();
    _hBosonMass_reco->Fill(temp.m());

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

    _hNjet->Fill((int)fastevent.jets().size());
    for (size_t i = 0; i < fastevent.jets().size();i++) {

      _hjetpT->Fill(fastevent.jets()[i]->pT());
      _hjeteta->Fill(fastevent.jets()[i]->eta());
      _hjetphi->Fill(fastevent.jets()[i]->phi());

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
 
  std::cout << " finishing my analysis " << endl;
  m_ROOToutFile->Write();


}

//==========================================================================

// You should not need to touch the main program: its actions are
// determined by the .cmnd file and the rest belongs in MyAnalysis.

int main(int argc, char* argv[]) {

  // Check that correct number of command-line arguments
  if (argc != 2) {
    cerr << " Unexpected number of command-line arguments. \n"
         << " You are expected to provide a file name and nothing else. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Check that the provided file name corresponds to an existing file.
  ifstream is(argv[1]);
  if (!is) {
    cerr << " Command-line file " << argv[1] << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Confirm that external file will be used for settings..
  std::cout << " PYTHIA settings will be read from file " << argv[1] << endl;

//  TApplication theApp("hist", &argc, argv);
  std::cout << " PYTHIA settings will be read from file " << argv[1] << endl;

  // Declare generator. Read in commands from external file.
  //Pythia8::Pythia pythia;
  //pythia.readFile(argv[1]);


  //fast_sim::FastSim _sim2;
  //_sim2.init(fast_sim::NOMINAL);

  // Initialization.
  //pythia.init();

  // Declare user analysis class. Do initialization part of it.
  //MyAnalysis myAnalysis;
  //myAnalysis.init();

  // Read in number of event and maximal number of aborts.
//  int nEvent = pythia.mode("Main:numberOfEvents");
//  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Begin event loop.
  int iAbort = 0;
  int iEvent = 0;

  int nelectrons=0;
  int nmuons = 0;
  int nchargedhads = 0;
  int nneutrinos = 0;
  int nphotons = 0;
  int nneutrals = 0;

  fast_sim::FastSim *_sim2; 
  Pythia8::Pythia *pythia;

//  #pragma omp parallel for reduction(+:nelectrons,nmuons,nchargedhads,nneutrinos,nphotons,nneutrals) default(none) shared(argv) private(iEvent,pythia,_sim2)

  std::vector<Pythia8::Pythia*> gen_list;
  std::vector<fast_sim::FastSim*>sim_list;

  std::vector<int> counts;
  std::vector<int> nevents;


  for (int k=0;k<6;k++) {
    pythia = new Pythia8::Pythia();
    pythia->readFile(argv[1]);
    pythia->init();
    pythia->readString("Random:setSeed=on");

    switch(k) {
      case 0: pythia->readString("Random:set=12345"); break;
      case 1: pythia->readString("Random:set=12346"); break;
      case 2: pythia->readString("Random:set=12347"); break;
      case 3: pythia->readString("Random:set=12348"); break;
      case 4: pythia->readString("Random:set=12349"); break;
      case 5: pythia->readString("Random:set=12350"); break;
      default: pythia->readString("Random:set=12310"); break;
    }

    gen_list.push_back(pythia);

    nevents.push_back(pythia->mode("Main:numberOfEvents"));
    counts.push_back(0);

    _sim2 = new fast_sim::FastSim();
    _sim2->init(fast_sim::NOMINAL);
    sim_list.push_back(_sim2);
  }


  TFile *m_ROOToutFile= new TFile("FastSim_Pythia8_output.root","RECREATE");
  TH1F *_hBosonMass_truth = new TH1F( "BosonMassTruth","Z Invariant Mass (Truth);GeV",100, 0, 200);
  TH1F *_hBosonMass_reco = new TH1F( "BosonMassReco","Z Invariant Mass (FastSim);GeV",100, 0, 200);
  TH1F *_hNelec_truth = new TH1F("NelecTruth","Number of Electrons (Truth);Number/Event",5,-0.5,4.5);
  TH1F *_hNelec = new TH1F("Nelec","Number of Isolated Electrons;Number/Event",5,-0.5,4.5);
  TH1F *_hlep1Pt = new TH1F("lep1Pt","Leading lep Pt;GeV;",100, 0., 200.);
  TH1F *_hlep1eta = new TH1F("lep1eta","Leading lep eta;",100, -5., 5.);
  TH1F *_hlep1phi = new TH1F( "lep1Phi","Leading lep Phi;",100, -6.0, 6.0);
  TH1F *_hlep1iso = new TH1F( "lep1iso","Leading lep Iso;",50, 0,10);

  TH1F *_hlep2Pt = new TH1F("lep2Pt","SubLeading lep Pt;GeV;",100, 0., 200.);
  TH1F *_hlep2eta = new TH1F("lep2eta","SubLeading  lep eta;",100, -5., 5.);
  TH1F *_hlep2phi = new TH1F( "lep2Phi","SubLeading lep Phi;",100, -6.0, 6.0);
  TH1F *_hlep2iso = new TH1F( "lep2iso","SubLeading lep Iso;",50, 0, 10.0);
  TH1F *_hNjet = new TH1F("Njet","Number of Jets;Number/Event",10,-0.5,9.5);
  TH1F *_hjetpT = new TH1F("jetPt","Jet Pt;GeV;",100, 0., 200.);
  TH1F *_hjeteta = new TH1F("jeteta","Jet eta;",100, -5., 5.);
  TH1F *_hjetphi = new TH1F( "jetPhi","Jet Phi;",100, -6.0, 6.0);

  vector<HEP_Simple_Lib::Particle*> _electrons;
  vector<HEP_Simple_Lib::Particle*> _muons;
  vector<HEP_Simple_Lib::Particle*> _photons;
  vector<HEP_Simple_Lib::Particle*> _bjets;
  vector<HEP_Simple_Lib::Particle*> _tauhads;
  vector<HEP_Simple_Lib::Particle*> _chargedhads;
  vector<HEP_Simple_Lib::Particle*> _weakly_interacting; // stdm neutrinos, susy neutralinos

  HEP_Simple_Lib::Particle* chosen;

  #pragma omp parallel num_threads(6) shared(_hBosonMass_truth,_hNelec_truth,_hBosonMass_reco,_hlep1Pt,_hlep1eta,_hlep1phi,_hlep1iso,_hlep2Pt,_hlep2eta,_hlep2phi,_hlep2iso,_hNjet,_hjetpT,_hjeteta,_hjetphi) private(_electrons,_muons,_photons,_bjets,_tauhads,_chargedhads,_weakly_interacting,chosen)
  {

    //vector<HEP_Simple_Lib::Particle*> _electrons;
 

  //int nEvent = pythia->mode("Main:numberOfEvents");
  //int nAbort = pythia->mode("Main:timesAllowErrors");


  //_sim2 = new fast_sim::FastSim();
  //_sim2->init(fast_sim::NOMINAL);

  std::cout << " number of threads is " << omp_get_thread_num() << std::endl;

// #pragma omp parallel for reduction(+:nelectrons,nmuons,nchargedhads,nneutrinos,nphotons,nneutrals) default(none) shared(argv,nEvent,pythia) private(iEvent)
  for (iEvent = 0; iEvent < nevents[0]; ++iEvent) {

    //printf("Event %d\n",iEvent);
    int thread = omp_get_thread_num(); 

    counts[thread] += 1; 
    // Generate events. Quit if too many failures.
    gen_list[thread]->next();

    //if (!pythia.next()) {
    //  if (++iAbort < nAbort) continue;
    //  cout << " Event generation aborted prematurely, owing to error!\n";
    //  break;
    //}

    // User Analysis of current event.

    //std::cout << " number of threads is " << omp_get_thread_num() << std::endl;
    for (int i = 0; i < gen_list[thread]->event.size(); ++i) {

      if (gen_list[thread]->event[i].isFinal()) {

        if (gen_list[thread]->event[i].isCharged()) {
          switch (int(fabs(gen_list[thread]->event[i].id()))) {
            /// @todo This needs to change, it should be only prompt leptons.. not just any lepton
            case 11: // electron
              
              nelectrons += 1;
              chosen = new HEP_Simple_Lib::Particle(gen_list[thread]->event[i].px(), gen_list[thread]->event[i].py(), gen_list[thread]->event[i].pz(), gen_list[thread]->event[i].e(), gen_list[thread]->event[i].id());
              _electrons.push_back(chosen);
              break;
            case 13: // muon

              nmuons += 1;
              chosen = new HEP_Simple_Lib::Particle(gen_list[thread]->event[i].px(), gen_list[thread]->event[i].py(), gen_list[thread]->event[i].pz(), gen_list[thread]->event[i].e(), gen_list[thread]->event[i].id());
              _muons.push_back(chosen);
              break;
            default: ;// every other hadronic charged particle - for the jets
              nchargedhads += 1; 
              chosen = new HEP_Simple_Lib::Particle(gen_list[thread]->event[i].px(), gen_list[thread]->event[i].py(), gen_list[thread]->event[i].pz(), gen_list[thread]->event[i].e(), gen_list[thread]->event[i].id());
              _chargedhads.push_back(chosen);
              //           printf("charged final particle missed %d\n",event[i].id());
          }
        }
        else {

          //       printf("neu particle %d\n",event[i].id());
          switch (int(fabs(gen_list[thread]->event[i].id()))) {
            case 22: // photon

              nphotons += 1;

              //            P2 = event[i].px()*event[i].px() + event[i].py()*event[i].py() + event[i].pz()*event[i].pz();

              //            printf("** P2 %f E %f\n",(float)P2,(float)(event[i].e()*event[i].e()));


              chosen = new HEP_Simple_Lib::Particle(gen_list[thread]->event[i].px(), gen_list[thread]->event[i].py(), gen_list[thread]->event[i].pz(), gen_list[thread]->event[i].id());
              _photons.push_back(chosen);
              break;
            case 12: // electron neutrinos
            case 14: // muon neutrinos
            case 16: // tau neutrinos
              nneutrinos += 1;
              chosen = new HEP_Simple_Lib::Particle(gen_list[thread]->event[i].px(), gen_list[thread]->event[i].py(), gen_list[thread]->event[i].pz(), gen_list[thread]->event[i].id());
              _weakly_interacting.push_back(chosen);
              //            printf("neutrinos %d\n",event[i].id());

              break;
            case -2112:
            case 2112: // neutrons
            case 130: // neutral kaon - not sure what to do
              nneutrals += 1;
              chosen = new HEP_Simple_Lib::Particle(gen_list[thread]->event[i].px(), gen_list[thread]->event[i].py(), gen_list[thread]->event[i].pz(), gen_list[thread]->event[i].e(), gen_list[thread]->event[i].id());
              _chargedhads.push_back(chosen);
              break;
            default: printf("neutral final particle missed %d\n",gen_list[thread]->event[i].id());
          }
        }
      }

      if ((gen_list[thread]->event[i].isQuark()) && (fabs(gen_list[thread]->event[i].id()) == 6)) {

        chosen = new HEP_Simple_Lib::Particle(gen_list[thread]->event[i].px(), gen_list[thread]->event[i].py(), gen_list[thread]->event[i].pz(), gen_list[thread]->event[i].e(), gen_list[thread]->event[i].id());
        _bjets.push_back(chosen);
      }


      /*
      // lets get our taus, they are unstable so we need a special list and case for them
      if ((gen_list[thread]->event[i].isLepton()) && (fabs(gen_list[thread]->event[i].id()) == 15)) {

        int daughter1 = pythia->event[i].daughter1();
        int daughter2 = pythia->event[i].daughter2();

        for (int j=daughter1; j<=daughter2;j++) {
          // printf("tau daughters (%d) are %d pdgId %d pt %.2f\n",j,event[j].id(),event[j].pT());
          std::cout << "tau daughters (" << j << ") are pdgId = " << pythia->event[j].id() << ", pt = " << pythia->event[j].pT() << " GeV" << endl;
        }

        // we need to remove the leptonically decaying taus - perhaps do an overlap removal with electrons
        // chosen = new HEP_Simple_Lib::Particle(event[i].px(), event[i].py(), event[i].pz(), event[i].e(), event[i].id());
        // _tauhads.push_back(chosen);
      }
      */
    }
    
    HEP_Simple_Lib::Event reco_event;
    sim_list[thread]->setParticles(_electrons, _muons, _photons, _chargedhads, _bjets, _tauhads, _weakly_interacting);
    sim_list[thread]->doDetectorResponse();
    sim_list[thread]->getRecoEvent(reco_event);

    #pragma omp critical
    {
       HEP_Simple_Lib::P4 z_p = _electrons[0]->mom()+_electrons[1]->mom();
       _hBosonMass_truth->Fill(z_p.m());
       _hNelec_truth->Fill((int)_electrons.size());

       std::cout << " number of threads is " << omp_get_thread_num() << " size of electrons " << _electrons.size() << std::endl;
       for (size_t kk = 0; kk < _electrons.size();kk++) 
         delete _electrons[kk]; 
       _electrons.clear();

       for (size_t kk = 0; kk < _muons.size();kk++) 
         delete _muons[kk]; 
       _muons.clear();

       for (size_t kk = 0; kk < _photons.size();kk++) 
         delete _photons[kk]; 
       _photons.clear();

       for (size_t kk = 0; kk < _chargedhads.size();kk++) 
         delete _chargedhads[kk]; 
       _chargedhads.clear();

       for (size_t kk = 0; kk < _bjets.size();kk++) 
         delete _bjets[kk]; 
       _bjets.clear();

        for (size_t kk = 0; kk < _tauhads.size();kk++) 
         delete _tauhads[kk]; 
       _tauhads.clear();

        for (size_t kk = 0; kk < _weakly_interacting.size();kk++) 
         delete _weakly_interacting[kk]; 
       _weakly_interacting.clear();


       _hNelec->Fill(reco_event.electrons().size());

       //cout << " the size of electrons " << reco_event.electrons().size() << endl;
       if ((int)reco_event.electrons().size() == 2) {
         HEP_Simple_Lib::P4 temp;
         temp = reco_event.electrons()[0]->mom() + reco_event.electrons()[1]->mom();
         _hBosonMass_reco->Fill(temp.m());

         if (reco_event.electrons()[0]->pT() >  reco_event.electrons()[1]->pT()) {
           _hlep1Pt->Fill(reco_event.electrons()[0]->pT());
           _hlep1eta->Fill(reco_event.electrons()[0]->mom().eta());
           _hlep1phi->Fill(reco_event.electrons()[0]->mom().phi());
           _hlep1iso->Fill(reco_event.electrons()[0]->isol());
           _hlep2Pt->Fill(reco_event.electrons()[1]->pT());
           _hlep2eta->Fill(reco_event.electrons()[1]->mom().eta());
           _hlep2phi->Fill(reco_event.electrons()[1]->mom().phi());
           _hlep2iso->Fill(reco_event.electrons()[1]->isol());
         }
         else {
           _hlep1Pt->Fill(reco_event.electrons()[1]->pT());
           _hlep1eta->Fill(reco_event.electrons()[1]->mom().eta());
           _hlep1phi->Fill(reco_event.electrons()[1]->mom().phi());
           _hlep2iso->Fill(reco_event.electrons()[1]->isol());
           _hlep2Pt->Fill(reco_event.electrons()[0]->pT());
           _hlep2eta->Fill(reco_event.electrons()[0]->mom().eta());
           _hlep2phi->Fill(reco_event.electrons()[0]->mom().phi());
           _hlep1iso->Fill(reco_event.electrons()[0]->isol());
         }

         _hNjet->Fill((int)reco_event.jets().size());
         for (size_t kk = 0; kk < reco_event.jets().size();kk++) {

           _hjetpT->Fill(reco_event.jets()[kk]->pT());
           _hjeteta->Fill(reco_event.jets()[kk]->eta());
           _hjetphi->Fill(reco_event.jets()[kk]->phi());

         }
       }
    }
    


//    myAnalysis.analyze( pythia->event);

  // End of event loop.
  }

  std::cout << "electrons " << nelectrons << " muons " << nmuons << std::endl;

  }


  m_ROOToutFile->Write();
  // Final statistics.

  for (int k=0;k<6;k++) {
    std::cout << " the counts for thread " << k << " is " << counts[k] << std::endl;
    gen_list[k]->stat();
    delete gen_list[k];
  }
  gen_list.clear();

  //delete pythia;
  //delete _sim2;

  // User finishing.
  //myAnalysis.finish();

  // Done.
  return 0;
}
