//  GAMBIT: Global and Modular BSM Inference Tool
//  //  ********************************************
//  //
//  //  Functions for FastSim
//  //
//  //  ********************************************
//  //
//  //  Authors
//  //  =======
//  //
//  //  (add name and date if you modify)
//  //  Kevin Le
//  //  2013 Feb
//  //  Aldo Saavedra
//  //  2013 March
//  //  2013 June 13,20
//  //  2014 Feb
//  //
//  //  ********************************************
//
//

#include "FastSim.hpp"
#include "fastjet/ClusterSequence.hh"
#include "simple_hep_lib/MathUtils.hpp"
#include "simple_hep_lib/MCUtils.hpp"
#include <algorithm>
#include <iostream>
#include <gsl/gsl_randist.h>

#include "FastSim_Logger.hpp"

static logging::logger log_inst(0);


using namespace std;
//using namespace fastjet;

using namespace hep_simple_lib;

namespace fast_sim {

  FastSim::FastSim() {
    _metx = -1.0;
    _mety = -1.0;
    _metx_truth = -1.0;
    _mety_truth = -1.0;

    _electron = NULL;
    _muon = NULL;
    _tau = NULL;
    _photon = NULL;
    _bjet = NULL;
    _jet = NULL;
    _rest = NULL;

    _random_num = gsl_rng_alloc(gsl_rng_mt19937);

    LOG_INFO("FastSim started:");
    //srand(1);
    // srand(time(NULL));

    //init(fast_sim::NOMINAL);
  }


  FastSim::~FastSim() {

    clear();

   for (size_t i = 0; i < _detector_perf.size(); i++) {

     if (_electron == _detector_perf[i])
       _electron = NULL;

     if (_muon == _detector_perf[i])
       _muon = NULL;

     if (_photon == _detector_perf[i])
       _photon = NULL;

     if (_tau == _detector_perf[i])
       _tau = NULL;

     if (_bjet == _detector_perf[i])
       _bjet = NULL;

     if (_jet == _detector_perf[i])
       _jet = NULL;

     if (_rest == _detector_perf[i])
       _rest = NULL;

     delete _detector_perf[i];
   }

   _detector_perf.clear();

   if (_electron != NULL)
     delete _electron;

   if (_muon != NULL)
     delete _muon;

   if (_tau != NULL)
     delete _tau;

   if (_jet != NULL)
     delete _jet;

   if (_bjet != NULL)
     delete _bjet;

   if (_rest != NULL)
     delete _rest;

   if (_photon != NULL)
     delete _photon;


  gsl_rng_free(_random_num);


    /*
    // for (size_t i = 0; i < _prompt_electrons.size();i++) delete _prompt_electrons[i];
    // _prompt_electrons.clear();
#define DELETE_PTRVEC(vec) for (size_t i = 0; i < vec.size();i++) delete vec[i]; vec.clear()

    DELETE_PTRVEC(_prompt_electrons);
    DELETE_PTRVEC(_prompt_muons);
    DELETE_PTRVEC(_prompt_photons);
    DELETE_PTRVEC(_bquarks);
    DELETE_PTRVEC(_bjets);
    DELETE_PTRVEC(_jets);
    DELETE_PTRVEC(_tauhads);
    DELETE_PTRVEC(_chargedhads);
    DELETE_PTRVEC(_weakly_interacting);

    DELETE_PTRVEC(_iso_electrons);
    DELETE_PTRVEC(_iso_muons);
    DELETE_PTRVEC(_iso_photons);

    DELETE_PTRVEC(_noniso_electrons);
    DELETE_PTRVEC(_noniso_muons);
    DELETE_PTRVEC(_noniso_photons);



#undef DELETE_PTRVEC
*/
  }




  void FastSim::Baseline_Response() {
    // this function has the baseline response 
    //

    // if the baseline was not defined in the json file then we need to specified default ones
    PProperties *new_pprop;
    if (_electron == NULL) {
      new_pprop = new PProperties;
      new_pprop->_pid = 11;
      new_pprop->_min_pt = 25.0; // 25.0 GeV
      new_pprop->_min_eta = -2.5;
      new_pprop->_max_eta = 2.5;
      new_pprop->_iso = 20000;
      new_pprop->_test_acceptance = false;

      _electron = new_pprop;
    }
    if (_muon == NULL) {
      new_pprop = new PProperties;
      new_pprop->_pid = 13;
      new_pprop->_min_pt = 20.0; // 20.0 GeV
      new_pprop->_min_eta = -2.5;
      new_pprop->_max_eta = 2.5;
      new_pprop->_iso = 20000;
      new_pprop->_test_acceptance = false;

      _muon = new_pprop;
    }
    if (_photon == NULL) {
      new_pprop = new PProperties;
      new_pprop->_pid = 22;
      new_pprop->_min_pt = 20.0; // 20.0 GeV
      new_pprop->_min_eta = -3.2;
      new_pprop->_max_eta = 3.2;
      new_pprop->_iso = 20000;
      new_pprop->_test_acceptance = false;

      _photon = new_pprop;
    }
    if (_bjet == NULL) {
      new_pprop = new PProperties;
      new_pprop->_pid = 5;
      new_pprop->_min_pt = 10.0; // 20.0 GeV
      new_pprop->_min_eta = -2.5;
      new_pprop->_max_eta = 2.5;
      new_pprop->_iso = 20000;
      new_pprop->_test_acceptance = false;

      _bjet = new_pprop;
    }
    if (_jet == NULL) {
      new_pprop = new PProperties;
      new_pprop->_pid = 0;
      new_pprop->_min_pt = 10.0; // 20.0 GeV
      new_pprop->_min_eta = -4.0;
      new_pprop->_max_eta = 4.0;
      new_pprop->_iso = 20000;
      new_pprop->_test_acceptance = false;

      _jet = new_pprop;
    }
    if (_rest == NULL) {
      new_pprop = new PProperties;
      new_pprop->_pid = 1;
      new_pprop->_min_pt = 0.4; // 400 MeV
      new_pprop->_min_eta = -4.0;
      new_pprop->_max_eta = 4.0;
      new_pprop->_iso = 20000;
      new_pprop->_test_acceptance = false;

      _rest = new_pprop;
    }
    if (_tau == NULL) {
      new_pprop = new PProperties;
      new_pprop->_pid = 15;
      new_pprop->_min_pt = 20.0; // 25.0 GeV
      new_pprop->_min_eta = -2.5;
      new_pprop->_max_eta = 2.5;
      new_pprop->_iso = 20000;
      new_pprop->_test_acceptance = false;

      _tau = new_pprop;
    }

  }


  void FastSim::init(DetectorType which,int debug_level) {


    log_inst.set_verbosity(debug_level);
    _simtype = which;
    switch(which) {
      case NOMINAL:
        // initialise the cuts and calorimeter limits and granularity for  the NOMINAL
        LOG_INFO("FastSim initialised with the nominal Detector");

        /*
        _min_muon_pt = 20.0;  // GeV
        _min_ele_pt  = 25.0;  // GeV
        _min_jet_pt  = 10.0;  //GeV
        _min_bjet_pt  = 10.0;  //GeV
        _min_photon_pt = 10; //GeV
        _min_track_pt = 0.4; //GeV

        _minEt_isol_muon = 4.0;
        _minEt_isol_electron = 4.0;
        _minEt_isol_photon = 4.0; // GeV


        _max_jet_eta = 4.0;
        _max_ele_eta = 2.5;
        _max_muon_eta = 2.5;
        _max_photon_eta = 3.2;
        _max_tauhad_eta = 2.5;
        _min_jet_eta = -4.0;
        _min_ele_eta = -2.5;
        _min_muon_eta = -2.5;
        _min_photon_eta = -3.2;
        _min_tauhad_eta = -2.5;
        */

        _calo_etamax = 5.0;  //GeV
        _calo_dphi =  0.01;
        _calo_deta = 0.01;
        _calo_transition = 3.2;
        _calo_etthresh = 0.0;

        _calo_neta = int(2*_calo_etamax/_calo_deta);
        _calo_nphi = int(2*3.2/_calo_dphi);

        _et_min = 0.2;//GeV
        _et_seedmin = 1.5;// GeV
        _cluster_rcone = 0.4;
        _cluster_etmin = 5.0;

        _fastjet = true; // use fastjet antikt or cones (clustering)
        _count = 0;

        Baseline_Response(); 
        break;

      default:
        /// @todo Throw exception here?
        ;
    }
  }

  void FastSim::init(std::string init_filename,int debug_level) {

    // this constructor uses the json datacard to initialise the properties of the detector
    // and the physics object whose response will be simulated

    log_inst.set_verbosity(debug_level);
    if (FastSim_Reader(init_filename) < 0) {

      LOG_ERR("Error reading the FastSim param_card");
    }
    LOG_INFO("FastSim initialised with:",init_filename);

    for (std::vector<PProperties*>::iterator it = _detector_perf.begin();it != _detector_perf.end(); it++) {

      if ((*it)->_level == "baseline") {

        switch ((*it)->_pid) {
        
          case 11: _electron = (*it); break;
          case 13: _muon =  (*it);break;
          case 15: _tau =  (*it);break;
          case 0: _jet =  (*it);break;
          case 5: _bjet =  (*it);break;
          case 22: _photon =  (*it);break;
          case 1: _rest =  (*it);break;
          default: ;
        }
      }
//      else {
//
//     LOG_ERR(" error on ",(*it)->_level);
//      }

    }

    // initialize anything that was missed by the initialization
    Baseline_Response(); 

  }


  void FastSim::clear() {


#define DELETE_PTRVEC(vec) for (size_t i = 0; i < vec.size();i++) delete vec[i]; vec.clear()

    DELETE_PTRVEC(_stable_interacting_particles);
    DELETE_PTRVEC(_weakly_interacting);

    _prompt_electrons.clear();
    _prompt_muons.clear();
    _prompt_photons.clear();
    _jets.clear();
    _tauhads.clear();
    _bquarks.clear();
    _bjets.clear();
    _chargedhads.clear();

#undef DELETE_PTRVEC

   /*
    _chargedhads.clear();
    _stable_interacting_particles.clear();
    _prompt_electrons.clear();
    _prompt_muons.clear();
    _prompt_photons.clear();
    _iso_electrons.clear();
    _iso_muons.clear();
    _iso_photons.clear();
    _noniso_electrons.clear();
    _noniso_muons.clear();
    _noniso_photons.clear();

    _weakly_interacting.clear();
    _bquarks.clear();
    _tauhads.clear();
    _jets.clear();
    _bjets.clear();
    */

    _cellietph.clear();
    _cellmom.clear();
    _cellhits.clear();
    _celleta.clear();
    _cellphi.clear();
    _cellswitch.clear();
    _clusncell.clear();
    _clusnhits.clear();
    _clusswitch.clear();
    _cluseta.clear();
    _clusphi.clear();
    _clusweta.clear();
    _cluswphi.clear();
    _clusmom.clear();

  }



  void FastSim::setParticles(vector<Particle*> electrons, vector<Particle*> muons,
      vector<Particle*> photons, vector<Particle*> nonprompt_leptons,
      vector<Particle*>charged_hadrons,
      vector<Particle*> bquarks, vector<Particle*> tauhads, vector<Particle*> weaklyint ) {

    clear();

    setElectrons(electrons);
    setMuons(muons);
    setPhotons(photons);

    setBQuarks(bquarks);
    setChargedHadrons(charged_hadrons);
    setNonPromptLeptons(nonprompt_leptons);
    setTauHads(tauhads);
    setWeaklyInteracting(weaklyint);

  }


  void FastSim::setElectrons(vector<Particle*> particles) {
    for (size_t i = 0; i < particles.size(); ++i) {
      if (abs(particles[i]->pid()) != 11) {

        LOG_WARN("Particles with PID ",particles[i]->pid()," found in the electron list!!");
        continue;
      }

      if (_rest == NULL)
        std::cout << "_rest is NULL" << std::endl;
      if ((particles[i]->eta() < _rest->_min_eta) or (particles[i]->eta() > _rest->_max_eta)) continue;
      if (particles[i]->pT() < _rest->_min_pt) continue;

      Particle* chosen = new Particle(particles[i]);
      _stable_interacting_particles.push_back(chosen);

      // _electron represents the baseline reconstruction
      if ((particles[i]->eta() < _electron->_min_eta) or (particles[i]->eta() > _electron->_max_eta)) continue;

      if (particles[i]->pT() < _electron->_min_pt) continue;
      LOG_DEBUG1("Chosen Electron ",i,particles[i]->pT(), particles[i]->eta()," cuts applied are ", _electron->_min_pt, _electron->_min_eta,_electron->_max_eta);
      chosen = new Particle(particles[i]);
      chosen->set_prompt(true);
      _prompt_electrons.push_back(chosen);
    }
  }


  void FastSim::setMuons(vector<Particle*> particles) {
    for (size_t i = 0; i < particles.size(); ++i) {
      if (abs(particles[i]->pid()) != 13) {
        LOG_WARN("Particles with PID ",particles[i]->pid()," found in the muon list!!");
        continue;
      }

      if ((particles[i]->eta() < _rest->_min_eta) or (particles[i]->eta() > _rest->_max_eta)) continue;
      if (particles[i]->pT() < _rest->_min_pt) continue;

      Particle* chosen = new Particle(particles[i]);
      _stable_interacting_particles.push_back(chosen);

      LOG_DEBUG1("Set muon ",i,particles[i]->pT(), particles[i]->eta()," cuts applied are ", _muon->_min_pt, _muon->_min_eta,_muon->_max_eta);
      if ((particles[i]->eta() < _muon->_min_eta) or (particles[i]->eta() > _muon->_max_eta)) continue;


      if (particles[i]->pT() < _muon->_min_pt) continue;

      LOG_DEBUG1("Chosen muon ",i,particles[i]->pT(), particles[i]->eta()," cuts applied are ", _muon->_min_pt, _muon->_min_eta,_muon->_max_eta);
      chosen = new Particle(particles[i]);
      chosen->set_prompt(true);
      _prompt_muons.push_back(chosen);
    }
  }


  void FastSim::setPhotons(vector<Particle*> particles) {
    for (size_t i = 0; i < particles.size(); ++i) {
      if (abs(particles[i]->pid()) != 22) {
        LOG_WARN("Particles with PID ",particles[i]->pid()," found in the photon list!!");
        continue;
      }

      if ((particles[i]->eta() < _rest->_min_eta) or (particles[i]->eta() > _rest->_max_eta)) continue;
      if (particles[i]->pT() < _rest->_min_pt) continue;

      Particle* chosen = new Particle(particles[i]);
      _stable_interacting_particles.push_back(chosen); //< @todo Build later?

      //LOG_DEBUG1("Set photons ",i,particles[i]->pT(), particles[i]->eta()," cuts applied are ", _photon->_min_pt, _photon->_min_eta,_photon->_max_eta);
      if ((particles[i]->eta() < _photon->_min_eta) or (particles[i]->eta() > _photon->_max_eta)) continue;

      if (particles[i]->pT() < _photon->_min_pt) continue;
      //LOG_DEBUG1("Chosen");
      chosen = new Particle(particles[i]);
      chosen->set_prompt(true);
      _prompt_photons.push_back(chosen);

      //cout << " photon pt " <<  particles[i]->pT() << " " << _min_photon_pt << " eta " << particles[i]->eta() << " phi " << particles[i]->phi() <<  endl;
    }

    //  printPhotons();
  }


  void FastSim::setBQuarks(vector<Particle*> particles) {
    for (size_t i = 0; i < particles.size(); ++i) {
      if (abs(particles[i]->pid()) != 5) {
        LOG_WARN("Particles with PID ",particles[i]->pid()," found b-quark list!!");
        continue;
      }

      //if ((particles[i]->eta() < _rest->_min_eta) or (particles[i]->eta() > _rest->_max_eta)) continue;
      //if (particles[i]->pT() < _rest->_min_pt) continue;

      //Particle* chosen = new Particle(particles[i]);
      //_stable_interacting_particles.push_back(chosen); //< @todo Build later?


      //if ((particles[i]->eta() < _bjet->_min_eta) or (particles[i]->eta() > _bjet->_max_eta)) continue;
      //if (particles[i]->pT() < _bjet->_min_pt) continue;

      // just keep the b-quark, we will match it to a jet later on
      Particle* chosen = new Particle(particles[i]);
      _bquarks.push_back(chosen);
    }
  }


  void FastSim::setTauHads(vector<Particle*> particles) {
    for (size_t i = 0; i < particles.size(); ++i) {
      if (abs(particles[i]->pid()) != 15) {
        LOG_WARN("Particles with PID ",particles[i]->pid()," found in the hadronic tau list!!");
        continue;
      }

      //if ((particles[i]->eta() < _min_tauhad_eta) or (particles[i]->eta() > _max_tauhad_eta)) continue;

      //Particle* chosen = new Particle(particles[i]);
      //_stable_interacting_particles.push_back(chosen); //< @todo Build later?

      if ((particles[i]->eta() < _tau->_min_eta) or (particles[i]->eta() > _tau->_max_eta)) continue;
      if (particles[i]->pT() < _tau->_min_pt) continue;

      //chosen = new Particle(particles[i]);
      //_tauhads.push_back(chosen);
    }
  }


  void FastSim::setChargedHadrons(vector<Particle*> particles) {
    for (size_t i = 0; i < particles.size(); ++i) {
      if (abs(particles[i]->pid()) == 11 || abs(particles[i]->pid()) == 13) {
        LOG_WARN("leptons ",particles[i]->pid()," found in the charged hadron list!!");
        continue;
      }

      if ((particles[i]->eta() < _rest->_min_eta) or (particles[i]->eta() > _rest->_max_eta)) continue;
      if (particles[i]->pT() < _rest->_min_pt) continue;

      Particle* chosen = new Particle(particles[i]);
      _stable_interacting_particles.push_back(chosen);
      _chargedhads.push_back(chosen);
    }
  }

  void FastSim::setNonPromptLeptons(vector<Particle*> particles) {
    
    for (size_t i = 0; i < particles.size(); ++i) {
      if ((abs(particles[i]->pid()) != 11) && (abs(particles[i]->pid()) != 13)) {
        LOG_WARN("non leptons ",particles[i]->pid(),particles[i]->pT(),particles[i]->eta()," found in non-prompt lepton list!!");
        continue;
      }

      if ((particles[i]->eta() < _rest->_min_eta) or (particles[i]->eta() > _rest->_max_eta)) continue;
      if (particles[i]->pT() < _rest->_min_pt) continue;

      Particle* chosen = new Particle(particles[i]);
      _stable_interacting_particles.push_back(chosen);

      LOG_DEBUG1("non-prompt lepton ",i,particles[i]->pid(),particles[i]->pT(),particles[i]->eta(),particles[i]->phi());
    }
  }


  void FastSim::setWeaklyInteracting(vector<Particle*> particles) {
    // its up to the user to fill the weakly interacting particles 
    // correctly 

    Particle* chosen;
    for (size_t i = 0; i < particles.size(); ++i) {

      switch(abs(particles[i]->pid())) {
        case 12:
        case 14:
        case 16:
          chosen = new Particle(particles[i]);
          _weakly_interacting.push_back(chosen);
          break;
        default:
          chosen = new Particle(particles[i]);
          _weakly_interacting.push_back(chosen);
          LOG_WARN("Particles with PID ",particles[i]->pid()," found in the weakly interacting list!!");
       }
    }
  }


  /// @brief CalculateMET determines the transerve missing energy associated with the event
  ///
  /// It just adds the x and y components of the weakly interacting particles
  /// @todo Rename to calcMET()
  void FastSim::calcMET_truth() {
    double totalx = 0.0;
    double totaly = 0.0;
    for (size_t i = 0; i <_weakly_interacting.size();i++) {
      totalx += _weakly_interacting[i]->mom().px();
      totaly += _weakly_interacting[i]->mom().py();
    }
    /// @todo Need to apply the detector resolution to this value
    _metx_truth = totalx;
    _mety_truth = totaly;
  }

  void FastSim::calcMET_CaloSum() {
    double totalx = 0.0;
    double totaly = 0.0;
    int ncell = _cellswitch.size();
    for (int i = 0; i < ncell; i++) {

      totalx += _cellmom[i]*cos(_cellphi[i]);
      totaly += _cellmom[i]*sin(_cellphi[i]);
    }

    /// @todo Need to apply the detector resolution to this value
    _metx = totalx;
    _mety = totaly;
  }




  // Return the missing ET associated with the event
  // it only includes the weakly interacting particles
  void FastSim::MET_truth(double &met, double &phi) {
    met = -99.0;
    phi = -99.0;
    if (_metx_truth < 0) // the MET has not been calculated for the event
      calcMET_truth();

    met = sqrt(_metx_truth*_metx_truth + _mety_truth*_mety_truth);
    phi = atan2(_metx_truth,_mety_truth);
    
  }

  // Return the missing ET associated with the event
  // it uses the sum of the  calorimeter cells
  void FastSim::MET(double &met, double &phi) {
    met = -99.0;
    phi = -99.0;
    if (_metx < 0) // the MET has not been calculated for the event
      calcMET_CaloSum();

    met = sqrt(_metx*_metx + _mety*_mety);
    phi = atan2(_metx,_mety);
    
  }

  // Return the y-component of the event's missing ET
  double FastSim::METx() {
    if (_metx < 0) // the MET has not been calculated for the event
      calcMET_CaloSum();
    return _metx;
  }


  // Return the x-component of the event's missing ET
  double FastSim::METy() {
    if (_metx < 0) // the MET has not been calculated for the event
      calcMET_CaloSum();
    return _mety;
  }


  // Return the phi-component of the event's MET vector
  double FastSim::METphi() {
    if (_metx < 0) // the MET has not been calculated for the event
      calcMET_CaloSum();
    return atan2(_metx,_mety);
  }

  // Return the y-component of the event's missing ET
  double FastSim::METx_truth() {
    if (_metx < 0) // the MET has not been calculated for the event
      calcMET_truth();
    return _metx_truth;
  }


  // Return the x-component of the event's missing ET
  double FastSim::METy_truth() {
    if (_metx < 0) // the MET has not been calculated for the event
      calcMET_truth();
    return _mety_truth;
  }


  // Return the phi-component of the event's MET vector
  double FastSim::METphi_truth() {
    if (_metx < 0) // the MET has not been calculated for the event
      calcMET_truth();
    return atan2(_metx_truth,_mety_truth);
  }




  /// @todo Rename
  void FastSim::doDetectorResponse() {
    // this function runs the fast simulation for the particles list provided
    //


    //cout << "detector response " << endl;
    // determine which cells which have been traversed by the interacting particles
    FindCells();


    // for the individual particles determine the response of detector on their momentum
    ElectronResponse();
    PhotonResponse();
    MuonResponse();

    // add all the neutrinos and determine the maginute of the
    // missing momentum and phi
    LOG_DEBUG1("Calculating MET");
    calcMET_CaloSum();
    calcMET_truth();


    // now for the jets, only clustering for now, fast jet will be incorporated soon
    switch(_simtype)
    {
      case NOMINAL:
        if (_fastjet) {
          //cout << " doing jet " << endl;
          
          LOG_DEBUG1("Using FastJet");
          JetResponse();
        }
        else {

          LOG_DEBUG1("Using Clustering");
          Clustering();
        }
        break;
      default:;
    }


    // determine which leptons and photons are isolated
    // so far it only checks whether they overlap which jets
    // need to check consitency, the jets are electrons
    //
    //cout << "isolation " << endl;
    AppliedIsolation();
    // find bjets using the bquark list and the reconstructed jets
    LOG_DEBUG1("Selecting Bjets");
    ChooseBJets();

    // the taus to be implemented

    //  printParticles();
    //  printSummary();

  }

  void FastSim::selectParticles(vector<Particle*> stable_particles, PProperties *cuts, vector<Particle*> &chosen_particles) {
    // this function will select the particles according to the specified cuts
    // assume that particles have passed the baseline cuts

    if (not cuts->_test_acceptance) {
      // there is no acceptance information so it is assume it is perfect, return all the particles

      if (stable_particles.size() > 0)
        LOG_WARN("No Acceptance found for ",stable_particles[0]->pid()," ideal response assumed");
      else
        LOG_WARN("No Acceptance found - ideal response is assumed");

      for (size_t i=0;i<stable_particles.size();i++) {
        Particle *temp_p = new Particle(stable_particles[i]);
        chosen_particles.push_back(temp_p);
      }
      return;
    }

    for (size_t i=0;i<stable_particles.size();i++) {

      bool passed_cuts = true;
      for (size_t cc= 0l; cc < cuts->_response.size();cc++) {

        switch(cuts->_response[cc]._type) {
          case PT: //1d map
            passed_cuts = isParticleMeasured(cuts->_response[cc],stable_particles[i]->pT());
            break;
          case ETA: //1d map
            passed_cuts = isParticleMeasured(cuts->_response[cc],stable_particles[i]->eta());
            break;
          case ISO: //1d map
            passed_cuts = isParticleMeasured(cuts->_response[cc],stable_particles[i]->isol());
            break;
          case PT_ETA: //2d map
            passed_cuts = isParticleMeasured(cuts->_response[cc],stable_particles[i]->pT(),stable_particles[i]->eta());
            break;
          default:
            LOG_WARN("Response Variable Not Implemented Yet");
        }
        if (!passed_cuts) 
          break;  // do not test the rest of the responses, it has failed
      }
      if (passed_cuts) {
        Particle *temp_p = new Particle(stable_particles[i]);
        chosen_particles.push_back(temp_p);
      }
    }
    return;
  }

  void FastSim::selectJets(vector<Jet*> jets, PProperties *cuts, vector<Jet*> &measured_jets) {
    // this function will select the jets according to the specified cuts
    // assume the jets have passed the baseline selection

    if (not cuts->_test_acceptance) {
      // there is no acceptance information so it is assume it is perfect, return all the particles

        LOG_WARN("No Acceptance found for jets, ideal response assumed for all jets: ", jets.size());

      for (size_t i=0;i<jets.size();i++) {
        Jet *temp_j = new Jet(jets[i]->mom(),jets[i]->isBJet());
        measured_jets.push_back(temp_j);
      }
      return;
    }

    for (size_t i=0;i<jets.size();i++) {

      bool passed_cuts = true;
      for (size_t cc= 0l; cc < cuts->_response.size();cc++) {

        switch(cuts->_response[cc]._type) {
          case PT: //1d map
            passed_cuts = isParticleMeasured(cuts->_response[cc],jets[i]->pT());
            break;
          case ETA: //1d map
            passed_cuts = isParticleMeasured(cuts->_response[cc],jets[i]->eta());
            break;
          case PT_ETA: //2d map
            passed_cuts = isParticleMeasured(cuts->_response[cc],jets[i]->pT(),jets[i]->eta());
            break;
          default:
            LOG_WARN("Response Variable Not Implemented Yet");
        }
        if (!passed_cuts) 
          break;  // do not test the rest of the responses, it has failed
      }
      if (passed_cuts) {
        Jet *temp_j = new Jet(jets[i]->mom(),jets[i]->isBJet());
        measured_jets.push_back(temp_j);
      }
    }
    return;
  }




  bool FastSim::isParticleMeasured(Acceptance acceptance, double test_value) {
    // this function determines whether the particle passes the cut

    if (acceptance._ndim != 1) {
      LOG_ERR("Error the number of dimensions of the acceptance",acceptance._name,acceptance._ndim," should be 1");
      return false; // the number of dimensions should be equal to 1
    }

    int chosen_bin=-1;
    double dice_roll;

    //std::cout << " name of acceptance " << acceptance._name << " test val " << test_value << std::endl;
    // check if the value of the test_value is greater than the highest x-limit
    if ((test_value >= acceptance._bin_edges_x.back()) || ( test_value < acceptance._bin_edges_x[0])) {

      //chosen_bin = (int)acceptance._bin_edges_x.size()-2; // bins are counted from 0
      //        std::cout << " pass the limit  " << acceptance._name << " test val " << test_value << " chosen bin " << chosen_bin << " " << acceptance._binvals.back()<< std::endl;
      chosen_bin = -1;

    }
    else {

      for (size_t k=0; k<acceptance._bin_edges_x.size()-1;k++) {

        if ((test_value >= acceptance._bin_edges_x[k]) && (test_value < acceptance._bin_edges_x[k+1])) {
          //          std::cout << " value is located in bin " << k << " between " << acceptance._bin_edges_x[k] << " " << acceptance._bin_edges_x[k+1]<<  " bin val " << acceptance._binvals[k] <<  std::endl;
          // we found the bin that our test_value corresponds to
          chosen_bin = (int)k;
          break;
        }
      }
    }

    if (chosen_bin == -1) { 
      // the particle 
      return true;
    }

    //dice_roll =  hep_simple_lib::closed_interval_rand(0.0,1.0);
    dice_roll = gsl_rng_uniform(_random_num);
    if (dice_roll >= acceptance._binvals[chosen_bin]) {
      //std::cout << dice_roll << " failed " << acceptance._binvals[chosen_bin] << std::endl;
      return false; // particle not measured
    }
    else {
      //std::cout << dice_roll << " passed " << acceptance._binvals[chosen_bin] << std::endl;
      return true; // particle measured
    }

  }

  bool FastSim::isParticleMeasured(Acceptance acceptance, double test_valuex, double test_valuey) {
    // this function determines whether the particle passes the cut

    if (acceptance._ndim != 2) {
      LOG_ERR("Error the number of dimensions of the acceptance",acceptance._name,acceptance._ndim," should be 2");
      return false; // the number of dimensions should be equal to 1
    }

    int chosen_bin_y=-1;
    int chosen_bin_x=-1;
    double dice_roll;

    //std::cout << " name of acceptance " << acceptance._name << " test val " << test_value << std::endl;
    // check if the value of the test_value is greater than the highest x-limit
    if ((test_valuex >= acceptance._bin_edges_x.back()) || ( test_valuex < acceptance._bin_edges_x[0]) || (test_valuey >= acceptance._bin_edges_y.back()) || ( test_valuey < acceptance._bin_edges_y[0])) {

      chosen_bin_y = -1;
      chosen_bin_x = -1;

   //      chosen_bin = (int)acceptance._bin_edges_x.size()-2; // bins are counted from 0
   //        std::cout << " pass the limit  " << acceptance._name << " test val " << test_value << " chosen bin " << chosen_bin << " " << acceptance._binvals.back()<< std::endl;
    }
    else {

      for (size_t k=0; k<acceptance._bin_edges_x.size()-1;k++) {

        if ((test_valuex >= acceptance._bin_edges_x[k]) && (test_valuex < acceptance._bin_edges_x[k+1])) {
          chosen_bin_x = (int)k;
          break;
        }
      }

      for (size_t l=0; l<acceptance._bin_edges_y.size()-1;l++) {

        if ((test_valuey >= acceptance._bin_edges_y[l]) && (test_valuey < acceptance._bin_edges_y[l+1])) {
          //          std::cout << " value is located in bin " << k << " between " << acceptance._bin_edges_x[k] << " " << acceptance._bin_edges_x[k+1]<<  " bin val " << acceptance._binvals[k] <<  std::endl;
          // we found the bin that our test_value corresponds to
          chosen_bin_y = (int)l;
          break;
        }
      }
    }

    if ((chosen_bin_x == -1) || ( chosen_bin_y == -1)) { 
      // the particle 
      return true;
    }

    //dice_roll =  hep_simple_lib::closed_interval_rand(0.0,1.0);
    dice_roll = gsl_rng_uniform(_random_num);

    int index = chosen_bin_x + (acceptance._bin_edges_y.size()-1)*chosen_bin_y;
    LOG_DEBUG1("IsParticle 2D valuex ",test_valuex,"bin_x ",chosen_bin_x," valuey ",test_valuey," bin_y",chosen_bin_y," 2d map index ", index,acceptance._binvals[index]);

    if (dice_roll >= acceptance._binvals[index]) {
      LOG_DEBUG1("IsParticle 2D Failed - random result ",dice_roll);
      return false; // particle not measured
    }
    else {
      LOG_DEBUG1("IsParticle 2D Pass - random result ",dice_roll);
      return true; // particle measured
    }

  }


  void FastSim::getRecoEvent(Event &event) {
    // this function fills the gambit event interface with the reconstructed particles

    // we will copy the prompt particles only - the isolated ones

    //cout << " the size of physics electrons " << _prompt_electrons.size() << endl;
    //printElectrons();
    event.add_particles(_prompt_electrons);
    event.add_particles(_prompt_muons);
    event.add_particles(_prompt_photons);

    //  for (size_t i=0;i<_iso_electrons.size();i++)
    //    event.addParticle(_iso_electrons[i]);

    //  for (size_t i=0;i<_iso_muons.size();i++)
    //    event.addParticle(_iso_muons[i]);

    //  for (size_t i=0;i<_iso_photons.size();i++)
    //    event.addParticle(_iso_photons[i]);

    // the MET
    event.set_missingmom(P4::mkXYZM(METx(), METy(), 0., 0.));
    event.set_missingmom_truth(P4::mkXYZM(METx_truth(), METy_truth(), 0., 0.));

    for (size_t i=0;i<_jets.size();i++)
      event.addJet(_jets[i]);



    /* still need to add the jets, bjets and the taus
       if(candidate->BTag)
       recoJet = new Jet(momentum.Px(), momentum.Py(), momentum.Pz(),
       momentum.E(), true);
       else
       recoJet = new Jet(momentum.Px(), momentum.Py(), momentum.Pz(),
       momentum.E(), false);
       event.addJet(recoJet);
       */


  }

  void FastSim::getRecoEvent(Event &event, std::string electron_category, float electron_isolation,
                                           std::string muon_category, float muon_isolation,
                                           std::string btag_category) {
    // this function fills the gambit event interface with the reconstructed particles

    // we will copy the prompt particles only - the isolated ones

    //cout << " the size of physics electrons " << _prompt_electrons.size() << endl;
    //printElectrons();

    std::vector<Particle*> measured_electrons, measured_muons;
    std::vector<Jet*> measured_bjets;

    // search for the 


    for (std::vector<PProperties*>::iterator it = _detector_perf.begin();it != _detector_perf.end(); it++) {

      switch ((*it)->_pid) {

        case 11: // the electron
        case -11:
                 if ((*it)->_level == electron_category) {
                   // an acceptance with a matching level has been found
                   selectParticles(_prompt_electrons,(*it),measured_electrons);
                 }
                 break;
        case 13: // the electron
        case -13:
                 if ((*it)->_level == muon_category) {
                   // an acceptance with a matching level has been found
                   selectParticles(_prompt_muons,(*it),measured_muons);
                 }
                 break;
        case -5:
        case 5: // b-tagging
                 if ((*it)->_level == btag_category) {
                   // an acceptance with a matching level has been found
                   selectJets(_bjets,(*it),measured_bjets);
                 }
                 break;
        default: ;
      }
    }

    //std::cout << " size of measured electrons is " << measured_electrons.size() << std::endl;
    LOG_DEBUG1("The number of measured prompt electrons is:",measured_electrons.size(),"muons:",measured_muons.size(),"bjets",measured_bjets.size());
    event.add_particles(measured_electrons);
    event.add_particles(measured_muons);
    event.add_particles(_prompt_photons);

    //  for (size_t i=0;i<_iso_electrons.size();i++)
    //    event.addParticle(_iso_electrons[i]);

    //  for (size_t i=0;i<_iso_muons.size();i++)
    //    event.addParticle(_iso_muons[i]);

    //  for (size_t i=0;i<_iso_photons.size();i++)
    //    event.addParticle(_iso_photons[i]);

    // the MET
    event.set_missingmom(P4::mkXYZM(METx(), METy(), 0., 0.));
    event.set_missingmom_truth(P4::mkXYZM(METx_truth(), METy_truth(), 0., 0.));

    // b-tagging and measured jets
    for (size_t i=0;i<_jets.size();i++) {

      bool btagged = false;
      for (size_t b=0;b<_bjets.size();b++) {

        if (_jets[i]->mom().deltaR_rap(measured_bjets[b]->mom()) < 0.2) { // this jet is actuall b-tagged
          btagged = true;
          break;
        }
      }
      _jets[i]->setBJet(btagged);
      event.addJet(_jets[i]);
    }
  }

  /// @todo Rename
  void FastSim::MuonResponse() {

    //std::vector<Particle*>::iterator newEnd;
    //cout << " the number of muons is " << _prompt_muons.size() << endl;
    for (size_t j = 0; j < _prompt_muons.size(); j++) {
      switch(_simtype)
      {
        case NOMINAL:
          _nodetector.MuonResponse(*_prompt_muons[j]);
          // Clustering();
          break;

        case ATLAS:
          _atlas_simple_response.MuonResponse(*_prompt_muons[j]);
          //      _prompt_muons.erase(std::remove(_prompt_muons[j],_prompt_muons[j],1);
          //      if not _atlas_simple_response.muonEfficiency(*_stable_muon[j]) {

          //        newEnd = std::remove(_prompt_muons.begin(), _prompt_muons.end(),_prompt_muons[j]);
          //        _prompt_muons.erase(newEnd);
          break;

        default:
          /// @todo Exception
          ;
      }
      }

  }

  void FastSim::ElectronResponse() {

    //std::cout << " electron " << std::endl;
    for (size_t j = 0; j <_prompt_electrons.size(); j++) {

      switch(_simtype)
      {
        case NOMINAL:
          _nodetector.ElectronResponse(*_prompt_electrons[j]);
          // Clustering();
          break;

        case ATLAS:
          _atlas_simple_response.ElectronResponse(*_prompt_electrons[j]);
          break;

        default:
          /// @todo Exception
          ;
      }
    }
  }

  /// @todo Rename
  void FastSim::PhotonResponse() {

    //std::cout << " photon " << std::endl;
    for (size_t j = 0; j < _prompt_photons.size(); j++) {
      switch(_simtype) {
        case NOMINAL:
          _nodetector.PhotonResponse(*_prompt_photons[j]);
          // Clustering();
          break;

        case ATLAS:
          _atlas_simple_response.PhotonResponse(*_prompt_photons[j]);
          break;

        default:
          /// @todo Exception
          ;
      }
    }

  }


  /// @todo Rename
  void FastSim::JetResponse() {
    // this functions needs implementing

    Jet *chosen;
    vector<fastjet::PseudoJet> calo_cells;
    int ncell = _cellswitch.size();
    double px,py,pz,E;

    // save the hit cell map

    //  char filename[50];
    //  sprintf(filename,"cell_map_%d.txt",_count);
    //  FILE* f = fopen(filename,"w");
    //cout << " number of cells is " << _cellswitch.size() << " " << _min_jet_pt << endl;
    for (int i = 0; i < ncell; i++) {

      px = _cellmom[i]*cos(_cellphi[i]);
      py = _cellmom[i]*sin(_cellphi[i]);
      pz = _cellmom[i]*sinh(_celleta[i]);
      E = sqrt(_cellmom[i]*_cellmom[i] + pz*pz);

      //if (E > 10.0)
      //  cout << "calo cell " << i << " pt " <<  _cellmom[i] << " E " << E << " phi " << _cellphi[i] << " eta " << _celleta[i] << " hits " << _cellhits[i] << endl;

      // an event with three particles:   px    py  pz      E
      calo_cells.push_back( fastjet::PseudoJet(px,py,pz,E));

      //cout << " calo cell size " << calo_cells.size() << endl;
    }
    //  fclose(f);
    //  _count++;

    // choose a jet definition
    double R = 0.4;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

    // run the clustering, extract the jets
    fastjet::ClusterSequence cs(calo_cells, jet_def);
    vector<fastjet::PseudoJet> alljets = sorted_by_pt(cs.inclusive_jets(_jet->_min_pt));


    fastjet::Selector sel = fastjet::SelectorEtaRange(-_jet->_max_eta,_jet->_max_eta);
    vector<fastjet::PseudoJet> jets = sel(alljets);

    //cout << " size of jets " << jets.size() << " all jets " << alljets.size() << endl;
    for (unsigned i = 0; i < jets.size(); i++) {

      P4 jet_vector = P4::mkXYZE(jets[i].px(),jets[i].py(),jets[i].pz(),jets[i].E());
      chosen = new Jet(jet_vector,false);
      _jets.push_back(chosen);


      //    vector<PseudoJet> constituents = jets[i].constituents();
      //    printf(" jet %d: %.2f %.2f %.2f n %d\n",i,(float)jets[i].perp(),(float)jets[i].eta(),(float)jets[i].phi(),(int)constituents.size());
      //    printf(" _jet %d: %.2f %.2f %.2f\n",i,(float)_jets[i]->pT(),(float)_jets[i]->eta(),(float)_jets[i]->phi());

      //    cout << "jet " << i << ": "<< jets[i].perp() << " "
      //                   << jets[i].eta() << " " << jets[i].phi() << endl;
      //    vector<PseudoJet> constituents = jets[i].constituents();
      //    for (unsigned j = 0; j < constituents.size(); j++) {
      //      cout << "    constituent " << j << "'s pt: " << constituents[j].perp()
      //          << endl;
      //   }
    }
    LOG_DEBUG1(" size of jets ",_jets.size());
  }


  void FastSim::ChooseBJets() {
    // this function will find the overlap between the jets and b quarks
    // jets which overlap will be categorised as b-jets

    // lets find which jets overlap with the bquarks
    for (size_t b = 0; b <_bquarks.size(); b++) {
      for (size_t j = 0; j <_jets.size(); j++) {

        LOG_DEBUG2("Choosing BJets - Bquark",b,_bquarks[b]->pT(),_bquarks[b]->eta(),_bquarks[b]->phi()," Jet ",j,_jets[j]->pT(),_jets[j]->eta(),_jets[j]->phi());
        if (_jets[j]->mom().deltaR_rap(_bquarks[b]->mom()) < 0.2) { // we found a b-quark and a reconstructed jet overlappling
          Jet *bjet = new Jet(_jets[j]->mom(),true);  // copy the jet and tagged as a b jet
          _bjets.push_back(bjet);
          break;
        }
      }
    }
  }


  /// @todo Rename
  void FastSim::FindCells() {
    // find the cells hit
    //

    /* for checking momentum conservation
       float sumx=0.0;
       float sumy=0.0;

       float sumz=0.0;
       float sume=0.0;
       */


    double eta,phi,pt;

    double sum_x = 0.0;
    double sum_y = 0.0;
    //  printf(" filling the calorimeter:\n");
    for (size_t i = 0; i < _stable_interacting_particles.size();i++) {

      //    if (_stable_interacting_particles[i]->pid() == 13)
      //      continue;   // we are not including them

      eta = _stable_interacting_particles[i]->mom().eta();// can modify particles eta position (move the jet)
      phi = _stable_interacting_particles[i]->mom().phi();
      pt = _stable_interacting_particles[i]->mom().pT();

      sum_x += _stable_interacting_particles[i]->mom().px();
      sum_y += _stable_interacting_particles[i]->mom().py();

      //cout << "id " << _stable_interacting_particles[i]->pid() << " pt " << pt << " eta " << eta << " phi " << phi << endl;
      fillcellvector(pt,eta,phi);

      /*
      double totalx = 0.0;
      double totaly = 0.0;
      for (size_t k = 0; k < _cellmom.size(); k++) { 
        totalx += _cellmom[k]*cos(_cellphi[k]);
        totaly += _cellmom[k]*sin(_cellphi[k]);

        cout << k <<  " x " << _cellmom[k]*cos(_cellphi[k]) << " y " <<_cellmom[k]*sin(_cellphi[k]) << " totalx " << totalx  << " totaly " << totaly << endl;



        cout << "_cellmon[" << k << "] " << _cellmom[k] << " _cellhits[" << k << "] = " << _cellhits[k] << " phi " << _cellphi[k] << " " << phi << endl;
      }
      cout << " size " << _cellmom.size() << " totalx " << totalx << " totaly " << totaly << endl;
      */

      // sumx += pt*cos(phi);
      // sumy += pt*sin(phi);

      // sumz += _stable_interacting_particles[i]->mom.pz();
      // sume += _stable_interacting_particles[i]->mom.E();
    }

    //cout << " the sum is " << sum_x << " " << sum_y << " " << sqrt(sum_x*sum_x + sum_y*sum_y) << endl;

    // for debugging
    // printf("Direct from event list:            px %10.2f  py %10.2f pz %10.2f E %10.2f\n", sumx, sumy, sumz, sume);

    // sumx = 0;
    // sumy = 0;
    // int ncell = cellswitch.size();

    // // check with the cells
    // printf("Cells Formed { \n");

    // for (int i = 0; i<ncell; i++) {

    // sumx += cellmom[i]*cos(cellphi[i]);
    // sumy += cellmom[i]*sin(cellphi[i]);

    // //printf("Cell %4d pt %10.2f %10.2f %10.2f %5d \n",i,cellmom[i],celleta[i],cellphi[i],(int)cellswitch[i]);
    // }
  }


  /// @todo Rename
  bool FastSim::CheckOverlap(Particle *p1, Particle *p2) {
    // determine if the particles overlap
    double dr = p1->mom().deltaR_eta(p2->mom());
    return (dr < _min_dr);
  }


  /// @todo Rename
  void FastSim::AppliedIsolation() {
    // this function will determine which are muons, electron and photons are isolated
    // this function should be called last once jets, and cells are populated
    // At the moment it should check how energy is in the cells surrounding
    // the particles.
    // What is below should be changed

    double ConeEt;
    // electrons
    for (size_t ee=0;ee<_prompt_electrons.size();ee++) {

      //    cout << "Electron eta " << _prompt_electrons[ee]->mom().eta() << " phi " <<_prompt_electrons[ee]->mom().phi() << endl;
      ConeEt = calcIsoEt(_prompt_electrons[ee]->mom().eta() ,_prompt_electrons[ee]->mom().phi());

      _prompt_electrons[ee]->set_isol(ConeEt);

      /*

      //cout << "E coneEt " << ConeEt << endl;
      if (ConeEt > _minEt_isol_electron) {
        _prompt_electrons[ee]->set_prompt(false);
//        Particle *temp_p = new Particle(_prompt_electrons[ee]);
//        _noniso_electrons.push_back(temp_p);
        isolated = false;
      }
      if (isolated) {
        cout << "Isolated Electron eta " << _prompt_electrons[ee]->mom().eta() << " phi " <<_prompt_electrons[ee]->mom().phi() << endl;
        _prompt_electrons[ee]->set_prompt(true);
        Particle *temp_p = new Particle(_prompt_electrons[ee]);
        _iso_electrons.push_back(temp_p);
      }
      */
    }


    // muon
    for (size_t mu=0;mu<_prompt_muons.size();mu++) {

      //    cout << "Muon eta " << _prompt_muons[mu]->mom().eta() << " phi " <<_prompt_muons[mu]->mom().phi() << endl;
      ConeEt = calcIsoEt(_prompt_muons[mu]->mom().eta() ,_prompt_muons[mu]->mom().phi());

      _prompt_muons[mu]->set_isol(ConeEt);

      /*
      if (ConeEt > _minEt_isol_muon) {
        _prompt_muons[mu]->set_prompt(false);
        Particle *temp_p = new Particle(_prompt_muons[mu]);
        _noniso_muons.push_back(temp_p);
        isolated = false;
      }
      if (isolated) {
        _prompt_muons[mu]->set_prompt(true);
        Particle *temp_p = new Particle(_prompt_muons[mu]);
        _iso_muons.push_back(temp_p);
      }
      */
    }

    // photon
    for (size_t ph=0;ph<_prompt_photons.size();ph++) {
      //    cout << "photon eta " << _prompt_photons[ph]->mom().eta() << " phi " <<_prompt_photons[ph]->mom().phi() << endl;
      ConeEt = calcIsoEt(_prompt_photons[ph]->mom().eta() ,_prompt_photons[ph]->mom().phi());
      //cout << "photon coneEt " << ConeEt << endl;

      _prompt_photons[ph]->set_isol(ConeEt);

      /*
      if (ConeEt > _minEt_isol_photon) {

        _prompt_photons[ph]->set_prompt(false);
        Particle *temp_p = new Particle(_prompt_photons[ph]);
        _noniso_photons.push_back(temp_p);
        isolated = false;
      }
      if (isolated) {

        cout << "Isolated photon eta " << _prompt_photons[ph]->mom().eta() << " phi " <<_prompt_photons[ph]->mom().phi() << endl;
        //
        _prompt_photons[ph]->set_prompt(true);
        Particle *temp_p = new Particle(_prompt_photons[ph]);
        _iso_photons.push_back(temp_p);
      }
      */
    }



    // electrons and jets
    /*
       bool isolated;
       for (int ee=0;ee<_prompt_electrons.size();ee++) {

       isolated = true;
       for (int jj=0;jj<_jets.size();jj++) {

       if (CheckOverlap(_prompt_electrons[ee],_jets[jj])) {
       _noniso_electrons.push_back(_prompt_electrons[ee]);
       isolated = false;
       continue;
       }
       }
       if (isolated)
       _iso_electrons.push_back(_prompt_electrons[ee]);
       }

    // muon and jets
    for (int mu=0;mu<_prompt_muons.size();mu++) {

    isolated = true;
    for (int jj=0;jj<_jets.size();jj++) {

    if (CheckOverlap(_prompt_muons[mu],_jets[jj])) {
    _noniso_muons.push_back(_prompt_muons[mu]);
    isolated = false;
    continue;
    }
    }
    if (isolated)
    _iso_muons.push_back(_prompt_muons[mu]);
    }

    // photon and jets
    for (int ph=0;ph<_prompt_photons.size();ph++) {

    isolated = true;
    for (int jj=0;jj<_jets.size();jj++) {

    if (CheckOverlap(_prompt_photons[ph],_jets[jj])) {
    _noniso_photons.push_back(_prompt_photons[ph]);
    isolated = false;
    continue;
    }
    }
    if (isolated)
    _iso_photons.push_back(_prompt_photons[ph]);
    }
    */

    // now lets check against the cells of the calorimeter
  }


  /// @todo Rename
  void FastSim::fillcellvector(double pt, double part_eta, double part_phi) {
    // this method digitises the position of the particle within
    // the calorimeter
    // a number of vectors are filled
    //

    int ieta, iphi;
    double cell_eta,cell_phi;

    //  cout << " phi " << phi << endl;
    // calculate the corresponding cell index of the particle using its eta and phi
    if (abs(part_eta) < _calo_transition) {
      ieta = int((part_eta + _calo_etamax) /(2 * _calo_etamax) * _calo_neta) + 1;
      iphi = int((part_phi + M_PI)/(2 * M_PI) * _calo_nphi) + 1;
    }else{ // in the barrel region the granularity is twice the barrel region granularity
      ieta = 2*int((part_eta + _calo_etamax)/(4 * _calo_etamax) * _calo_neta) + 1;
      iphi = 2*int((part_phi + M_PI)/(4 * M_PI) * _calo_nphi) + 1;
    }

    // calculate the corresponding index of the cell using their eta and phi position and the number phi
    // channels, it is a cell id
    int   ietph = _calo_nphi * ieta + iphi;

    // check if the cell already exists
    bool exists = false;
    for (int i=0; i<(int)_cellietph.size(); i++) {
      if (_cellietph[i] == ietph) {  // it exists add the momentum and increment the hit on the cell
        _cellmom[i] += pt;
        _cellhits[i] += 1;
        exists = true;
      }
    }

    //printf("%5d %6.3f %6.3f ", ietph,  eta, phi);

    // calculate the digitized position (eta,phi) of the particle
    if (abs(part_eta) < _calo_transition) {
      cell_eta = 2.0 * _calo_etamax * ((ieta) - 1.0 + 0.5)/_calo_neta  - _calo_etamax;
      cell_phi = 2.0 * M_PI * ((iphi)-1.0 + 0.5)/_calo_nphi - M_PI;
    }
    else{
      cell_eta = 2.0 * _calo_etamax * ((ieta) - 1.0 + 1.0)/_calo_neta  - _calo_etamax;
      cell_phi = 2.0 * M_PI * ((iphi)-1.0 + 1.0)/_calo_nphi - M_PI;
    }

    //  printf("particle %6.3f %6.3f  calo %6.3f %6.2f\n", part_eta, part_phi,eta,phi);

    // first time cell has been hit
    if (not exists) {
      _cellietph.push_back(ietph); // the cell id
      _cellmom.push_back(pt);      // transverse momentum
      _cellhits.push_back(1);      // number of hits
      _celleta.push_back(cell_eta);     // digitized eta
      _cellphi.push_back(cell_phi);     // digitized phi
      _cellswitch.push_back(2);    // for tracking purposes
      //    printf("cell eta diff %.4f phi diff %.4f %.2f\n",cell_eta-part_eta,cell_phi-part_phi,pt);
    }
  }


  /// @todo Rename
  void FastSim::Clustering() {
    // this function takes the cells, creates clusters from seed cells

    int ndum, ncell, imax;
    double etmax, phi, eta, etac, phic, dr, etad, phid; //dphia

    _clusncell.clear();
    _clusnhits.clear();
    _clusswitch.clear();
    _cluseta.clear();
    _clusphi.clear();
    _clusweta.clear();
    _cluswphi.clear();
    _clusmom.clear();

    int done = 0;
    ncell = _cellswitch.size();
    while(done == 0) {


      ndum = 0;
      imax=-1;
      etmax = 0;
      eta=-999;
      phi=-999;
      //find unused cell with highest energy
      //seed for the cluster
      for (int i = 0; i < ncell; i++) {
        if ((_cellswitch[i]==2)&&(_cellmom[i]>etmax)) {
          if (imax != -1) _cellswitch[imax] =2;
          imax = i;
          eta = _celleta[i];
          phi = _cellphi[i];
          etmax = _cellmom[i];
          _cellswitch[i] = 1;
          //printf(" seed cells %5d %6.2f %6.2f %6.2f\n",i,etmax,eta,phi);
        }
      }
      //printf("etmax %.2f  etini %.2f %d %d %d\n",etmax,etini,imax,ncell,cellswitch.size());
      //printf("cell switch size %d\n", cellswitch.size());
      //  if chosen cell does meet the ET requirement use as SEED
      if (etmax >= _et_seedmin) {

        if (((int)_cellswitch.size() > imax) && (imax >= 0)) {

          // initialize the cluster
          int  temp_clusncell = 1;
          int  temp_clusnhits = 0;
          double  temp_clusweta  = 0.0;
          double  temp_cluswphip  = 0.0;
          double  temp_cluswphin  = 0.0;
          double  temp_cluseta = eta;
          double  temp_clusphi = phi;
          double  temp_clusmomp   = 0.0;
          double  temp_clusmomn   = 0.0;


          // Add unused cells within distance of SEED
          for (int i = 0; i < ncell; i++) {
            if (_cellswitch[i]!=0) {
              phic  = _cellphi[i];
              etac  = _celleta[i];
              etad = etac-temp_cluseta;
              phid = hep_simple_lib::delta_phi(phic, temp_clusphi);

              dr = sqrt(pow(etad,2)+ pow(phid,2));

              if (dr < _cluster_rcone) {
                _cellswitch[i]   = -_cellswitch[i];
                temp_clusncell += 1;
                temp_clusnhits += _cellhits[i];
                temp_clusweta  += (etac) * _cellmom[i]; //
                if (phic>0) {
                  temp_cluswphip  += (phic+10.0) * _cellmom[i];
                  temp_clusmomp   += _cellmom[i];

                }
                if (phic<0) {
                  temp_cluswphin  += (phic+10.0) * _cellmom[i];
                  temp_clusmomn   += _cellmom[i];
                }
              }
            }
          }
          //printf("cluster momentum %6.2f\n",temp_clusmom);
          // does the cluster have the minimum energy required
          if (temp_clusmomp + temp_clusmomn>_cluster_etmin) {
            // save the cluster
            ndum = ndum + 1;
            _clusncell.push_back(temp_clusncell);
            _clusnhits.push_back(temp_clusnhits);
            _clusswitch.push_back(1);
            _cluseta.push_back(temp_cluseta);
            _clusphi.push_back(temp_clusphi);
            if (temp_clusmomp + temp_clusmomn > 0) {
              _clusweta.push_back(temp_clusweta/(temp_clusmomp + temp_clusmomn));

              double smartp = temp_cluswphip/temp_clusmomp;
              double smartn = temp_cluswphin/temp_clusmomn;

              double smart = smartp + (smartp+ smartn)*temp_clusmomp/(temp_clusmomp+temp_clusmomn);
              // if (smart>pi) smart += -2* pi;

              _cluswphi.push_back(smart);
            }
            else {
              _clusweta.push_back(-99.0);
              _cluswphi.push_back(-99.0);
            }


            _clusmom.push_back(temp_clusmomp+temp_clusmomn);

            for (int i = 0; i < ncell; i++) {
              if (_cellswitch[i]<0) _cellswitch[i] = 0;
            }
          }
          else {
            for (int i = 0; i < ncell; i++) {
              if (_cellswitch[i]<0) _cellswitch[i] = -_cellswitch[i];
            }
          }
        }
        else {
          cout << "imax " << imax << " size of cellswitch " << (int)_cellswitch.size() << " " << ncell << endl;
          return;
        }
      }
      else
      {
        done = 1;
      }
    }

    // now that we have a list of the clusters, they will be our jets

    Jet *chosen;
    double px,py,pz;
    for (size_t i = 0; i <_clusmom.size();i++) {

      // if the jet pt is less than the minimum discarded
      if  (_clusmom[i] <= _jet->_min_pt)
        continue;

      if  (fabs(_clusweta[i]) <= _jet->_max_eta)
        continue;

      px = _clusmom[i]*cos(_cluswphi[i]);
      py = _clusmom[i]*sin(_cluswphi[i]);
      pz = _clusmom[i]*sinh(_clusweta[i]);

      P4 jet_vector = P4::mkEtaPhiME(_clusweta[i], _cluswphi[i],0.0,sqrt(px*px+py*py+pz*pz));
      chosen = new Jet(jet_vector,false);

      _jets.push_back(chosen);

    }
  }

  double FastSim::calcIsoEt(double or_eta, double or_phi) {
    /* add the transverse energy on the cells within a deltaR of eta and phi of the particle
     * and return the value
     */

    int ncell;
    double dphi,deta,dR,SumEt,eta,phi;
    SumEt = 0;
    ncell = _cellswitch.size();
    for(int i = 0; i<ncell; i++){
      eta = _celleta[i];
      phi = _cellphi[i];


      dphi = hep_simple_lib::delta_phi(or_phi,phi);
      deta = fabs(or_eta - eta);
      dR = sqrt(dphi*dphi + deta*deta);

      if ((dR > 0.05) && (dR < 0.2))
        SumEt += _cellmom[i];
    }

    return SumEt;
  }


  void FastSim::printParticles() {

    printMuons();
    printElectrons();
    printPhotons();
  }


  /// @todo Rename
  void FastSim::printMuons() {
    LOG_INFO("Number of muons is ",_prompt_muons.size());
    for (int j=0;j<(int)_prompt_muons.size();j++) {
      LOG_INFO("Muons ",j,_prompt_muons[j]->pid()," P: ", _prompt_muons[j]->mom().px(),_prompt_muons[j]->mom().py(),_prompt_muons[j]->mom().pz(),
        _prompt_muons[j]->mom().E()," Eta: ",_prompt_muons[j]->mom().eta()," Phi: ",_prompt_muons[j]->mom().phi());
    }
  }


  /// @todo Rename
  void FastSim::printElectrons() {

    LOG_INFO("Number of electrons is ",_prompt_electrons.size());
    for (int j=0;j<(int)_prompt_electrons.size();j++) {
      LOG_INFO("Electrons ",j,_prompt_electrons[j]->pid()," P: ", _prompt_electrons[j]->mom().px(),_prompt_electrons[j]->mom().py(),_prompt_electrons[j]->mom().pz(),
        _prompt_electrons[j]->mom().E()," Eta: ",_prompt_electrons[j]->mom().eta()," Phi: ",_prompt_electrons[j]->mom().phi());
    }
  }


  /// @todo Rename
  void FastSim::printPhotons() {
    LOG_INFO("Number of photons is ",_prompt_photons.size());
    for (int j=0;j<(int)_prompt_photons.size();j++) {
      LOG_INFO("Photons ",j,_prompt_photons[j]->pid()," P: ", _prompt_photons[j]->mom().px(),_prompt_photons[j]->mom().py(),_prompt_photons[j]->mom().pz(),
        _prompt_photons[j]->mom().E()," Eta: ",_prompt_photons[j]->mom().eta()," Phi: ",_prompt_photons[j]->mom().phi());
    }
  }

  void FastSim::printSummary() {

    LOG_INFO("Number of prompt electrons is ",_prompt_electrons.size() );
    LOG_INFO("Number of prompt muons is ",_prompt_muons.size() );
    LOG_INFO("Number of prompt photons is ",_prompt_photons.size() );

  }
  }
