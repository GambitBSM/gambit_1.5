//  GAMBIT: Global and Modular BSM Inference Tool
//  //  ********************************************
//  //
//  //  Class Definition for FastSim
//  //
//  //  ********************************************
//  //
//  //  Authors
//  //  =======
//  //
//  //  (add name and date if you modify)
//  //
//  //  Aldo Saavedra
//  //  2013 March
//  //  2013 June 13,20
//  //  2014 Feb
//  //
//  //  ********************************************
//
//

#ifndef __FASTSIM_HPP__
#define __FASTSIM_HPP__

#include "simple_hep_lib/Particle.hpp"
#include "simple_hep_lib/Jet.hpp"
#include "simple_hep_lib/Event.hpp"                 // Gambit's event interface
#include "DetectorResponse.hpp"
#include "json/json.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

namespace fast_sim {

  enum DetectorType { NOMINAL,ATLAS, CMS, EXPERIMENTAL};
  // NOMINAL the energy response is smeared and the calorimeter parameters are the same as the ACERDET
  // ACERDET the energy of the different particles is smeared by the respective functions as used in the ACERDET paper and the calorimeter parameters are
  //         the same as the ACERDET paper
  // the values below are not implemented yet
  // FASTATLAS  will be defined with a simple smear functions but the calorimeter parameters will reflect the ATLAS detector
  // FASTCMS    will be defined with a simple smear functions but the calorimeter parameters will reflect the CMS detector
  // ATLAS2011  will try to reproduce the performance of the ATLAS detector
  // CMS2011    will try to reproduce the performance of the CMS detector

  // predefined responses that can be used
  enum Response_Type { PT, ETA, PT_ETA, ISO };

/*
  typedef struct {

    std::vector<hep_simple_lib::Jet*> _jets;
    std::vector<hep_simple_lib::Jet*> _bjets;

    std::vector<hep_simple_lib::Particle*> _iso_electrons;
    std::vector<hep_simple_lib::Particle*> _iso_muons;
    std::vector<hep_simple_lib::Particle*> _iso_photons;
    std::vector<hep_simple_lib::Particle*> _noniso_electrons;
    std::vector<hep_simple_lib::Particle*> _noniso_muons;
    std::vector<hep_simple_lib::Particle*> _noniso_photons;

    double _MET;
    double _METPhi;

  } FastSimEvent;
  */

  struct Acceptance {
    bool _init;
    Response_Type _type;
    std::string _name;
    int _nbinsx;
    int _nbinsy;
    int _ndim;
    std::vector<double> _bin_edges_x;
    std::vector<double> _bin_edges_y;
    std::vector<double> _binvals;
  };

  struct PProperties {

    int _pid;
    std::string _level;
    double _min_pt;
    double _min_eta;
    double _max_eta;
    double _iso;
    bool _test_acceptance; // determines whether the particle has an acceptance or not
    std::vector<Acceptance> _response;
  };


  class FastSim {
    public:
      FastSim();
      ~FastSim();

      
      void init(DetectorType which,int debug_level);
      void init(std::string filename,int debug_level);


      void Baseline_Response();

      
      int FastSim_Reader(std::string filename);
      int FastSim_ObjectReader(Json::Value phys_object, PProperties &particle_props);


      int FastSim_AcceptanceItemReader(const Json::Value reco_object, Acceptance &efficiency);
//      int FastSim_AcceptanceItemReader(const Json::Value reco_object, Acceptance_2D &efficiency);
      int FastSim_InfoReader(const Json::Value histo_object, int &ndim, int &nbinsx, int &nbinsy,
      std::vector<double> &bin_edgesx, std::vector<double> &bin_edgesy, std::vector<double> &histo_values);
      int FastSim_ArrayReader(const Json::Value array, std::vector<double> &values );
      


      void doDetectorResponse();
      void FindCells();

      // this functions set the particle list for each type
      void setParticles(std::vector<hep_simple_lib::Particle*> electrons, std::vector<hep_simple_lib::Particle*> muons,
          std::vector<hep_simple_lib::Particle*> photons,std::vector<hep_simple_lib::Particle*> nonprompt_leptons,
          std::vector<hep_simple_lib::Particle*>charged_hadrons,
          std::vector<hep_simple_lib::Particle*> bjets, std::vector<hep_simple_lib::Particle*> tauhads, std::vector<hep_simple_lib::Particle*> weaklyint );
      void setBQuarks(std::vector<hep_simple_lib::Particle*> particles);
      void setElectrons(std::vector<hep_simple_lib::Particle*> particles);
      void setMuons(std::vector<hep_simple_lib::Particle*> particles);
      void setPhotons(std::vector<hep_simple_lib::Particle*> particles);
      void setTauHads(std::vector<hep_simple_lib::Particle*> particles);
      void setChargedHadrons(std::vector<hep_simple_lib::Particle*> particles);
      void setNonPromptLeptons(std::vector<hep_simple_lib::Particle*> particles);
      void setWeaklyInteracting(std::vector<hep_simple_lib::Particle*> particles);

      // acceptance methods
      void selectParticles(std::vector<hep_simple_lib::Particle*> stable_particles, PProperties *cuts, std::vector<hep_simple_lib::Particle*> &chosen_particles);

      void selectJets(std::vector<hep_simple_lib::Jet*> jets, PProperties *cuts, std::vector<hep_simple_lib::Jet*> &measured_jets);
      bool isParticleMeasured(Acceptance acceptance, double test_value_x,double test_value_y);
      bool isParticleMeasured(Acceptance acceptance, double test_value_x);

      void clear();

      void Clustering();

      void ElectronResponse();
      void MuonResponse();
      void PhotonResponse();
      void JetResponse();
      void ChooseBJets();
      void AppliedIsolation();

      double calcIsoEt(double or_eta, double or_phi);

      bool CheckOverlap(hep_simple_lib::Particle *p1, hep_simple_lib::Particle *p2);

      void calcMET_CaloSum();
      void calcMET_truth();
      void MET_truth(double &met, double &phi);
      void MET(double &met, double &phi);


      double METx();
      double METy();
      double METphi();
      double METx_truth();
      double METy_truth();
      double METphi_truth();



      void printMuons();
      void printElectrons();
      void printPhotons();
      void printSummary();
      void printParticles();

      void fillcellvector(double pt, double eta, double phi);
      int NElectrons(){ return _prompt_electrons.size(); }
      int NMuons(){ return _prompt_muons.size(); }


      void getRecoEvent(hep_simple_lib::Event &event);
      void getRecoEvent(hep_simple_lib::Event &event, std::string electron_category, float electron_isolation, std::string muon_category, float muon_isolation,
                        std::string btag_category);
      
    private:

      // calorimeter range and granularity
      double _calo_dphi;  // the delta phi of the calorimeter in the barrel region
      double _calo_deta; // the delta eta of the calorimeter in the barrel region
      double _calo_transition;  // the eta position of the transition region
      double _calo_etthresh; // the energy threshold of a calorimeter cell
      double _calo_bfield_ptmin; // the minimum pt that the b field
      double _calo_etamax; // the total coverage of the calorimeter
      int    _calo_neta; // number of eta channels in the calorimeter
      int    _calo_nphi; // number of phi channels in the calorimeter

      double _et_min; // the minimum energy of a cell
      double _et_seedmin; // the minimum transverse energy a cell can have to be a seed for a cluster
      double _cluster_rcone;
      double _cluster_etmin; // the minimum energy that a cluster can have otherwise discarded
      double _min_dr;

      bool _fastjet;
      int _count;

      // the particles
      /*
      double _min_muon_pt;
      double _min_ele_pt;
      double _min_photon_pt;
      double _min_bjet_pt;
      double _min_jet_pt;
      double _min_tauhad_pt;
      double _min_dr;

      double _min_track_pt; //GeV

      double _max_jet_eta;
      double _max_bjet_eta;
      double _max_ele_eta;
      double _max_muon_eta;
      double _max_photon_eta;
      double _max_tauhad_eta;
      double _min_jet_eta;
      double _min_bjet_eta;
      double _min_ele_eta;
      double _min_muon_eta;
      double _min_photon_eta;
      double _min_tauhad_eta; 

      // the isolation
      double _minEt_isol_muon; 
      double _minEt_isol_electron; 
      double _minEt_isol_photon; */


      std::vector<PProperties*> _detector_perf;

      // the loose properties, the baseline cuts
      PProperties* _electron;
      PProperties* _muon;
      PProperties* _tau;
      PProperties* _photon;
      PProperties* _jet;
      PProperties* _bjet;
      PProperties* _rest;

      
      // the particles
      std::vector<hep_simple_lib::Particle*> _chargedhads;
      std::vector<hep_simple_lib::Particle*> _stable_interacting_particles;
      std::vector<hep_simple_lib::Particle*> _prompt_electrons;
      std::vector<hep_simple_lib::Particle*> _prompt_muons;
      std::vector<hep_simple_lib::Particle*> _prompt_photons;
      //std::vector<hep_simple_lib::Particle*> _iso_electrons;
      //std::vector<hep_simple_lib::Particle*> _iso_muons;
      //std::vector<hep_simple_lib::Particle*> _iso_photons;
      //std::vector<hep_simple_lib::Particle*> _noniso_electrons;
      //std::vector<hep_simple_lib::Particle*> _noniso_muons;
      //std::vector<hep_simple_lib::Particle*> _noniso_photons;
      std::vector<hep_simple_lib::Particle*> _weakly_interacting;
      std::vector<hep_simple_lib::Particle*> _bquarks;
      std::vector<hep_simple_lib::Particle*> _tauhads;
      std::vector<hep_simple_lib::Jet*> _jets;
      std::vector<hep_simple_lib::Jet*> _bjets;
      
      // measurements
      double _metx;
      double _mety;
      double _metx_truth;
      double _mety_truth;


      
      DetectorType _simtype;
      DetectorResponse _nodetector;
      ATLAS_Simple_Response _atlas_simple_response;
      

      /*
         double etacel = 5.0; //rapidity coverage
         double ptmin = 0.5; // min pt for b field
         double etthr = 0.0; //ein et for cell
         double caloth = 3.2; // eta transition in cells granularity
         double dbeta = 0.01; // granularity of eta,  twice outside
         double dbphi = 0.01; // granulrity of phi, twice outside
         */

      // list of the cells hit
      //  vector<double> cellietph(0), cellmom(0), cellhits(0), celleta(0), cellphi(0), cellswitch(0);
      std::vector<double> _cellietph,_cellmom,_cellhits,_celleta,_cellphi,_cellswitch;
      std::vector<double> _clusncell,_clusnhits,_clusswitch,_cluseta,_clusphi,_clusweta,_cluswphi,_clusmom;

      // for the random numbers using the gnu scientific library
      gsl_rng *_random_num;

  };

//  }

}

#endif
