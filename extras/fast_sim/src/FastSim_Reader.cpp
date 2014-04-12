//  GAMBIT: Global and Modular BSM Inference Tool
//  //  ********************************************
//  //
//  //  FastSim_Reader
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
//  //  2014 March
//  //  
//  //  This reads the json file and initialises the detector
//  //  ********************************************
//
//



#include "FastSim.hpp"
#include "json/json.h"

#include <algorithm>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>


namespace fast_sim {

  int FastSim::FastSim_Reader(std::string init_filename) {
    // this functuon reads the jason file and initialises the class

    std::ifstream fastsim_datacard(init_filename.c_str());
    if (not fastsim_datacard.good()) {
      // could not open the file
      return -1;
    }

    std::stringstream buffer;
    buffer << fastsim_datacard.rdbuf();

    Json::Value root;   // will contains the root value after parsing.
    Json::Reader reader;
    bool parsingSuccessful = reader.parse( buffer.str(), root );
    if ( !parsingSuccessful )
    {
      // report to the user the failure and their locations in the document.
      std::cout  << "Failed to parse configuration\n"
        << reader.getFormattedErrorMessages();
      return -1;
    }

    // lets start reading the datacard, if detector field empty use ATLAS
    std::string type_str,level_str;
    double min_pt,min_eta,max_eta,iso_cut;
    const Json::Value detector = root["detector"];

    std::string which_detector = detector.get("name","ATLAS").asString();
    std::transform(which_detector.begin(), which_detector.end(),which_detector.begin(), ::tolower);
    if (which_detector == "nominal")
      _simtype = NOMINAL;
    else if (which_detector == "atlas")
      _simtype = ATLAS;
    else if (which_detector == "cms")
      _simtype = CMS;
    else if (which_detector == "x")
      _simtype = EXPERIMENTAL;
    else
      _simtype = ATLAS;

    _calo_etamax = detector.get("calo_eta_max",5.0).asDouble();
    _calo_dphi = detector.get("calo_dphi",0.01).asDouble();
    _calo_deta = detector.get("calo_deta",0.01).asDouble();
    _calo_transition = detector.get("calo_transition",3.2).asDouble();
    _et_min = detector.get("et_min",0.2).asDouble();
    _et_seedmin = detector.get("et_seed_min",1.5).asDouble();
    _cluster_rcone = detector.get("cluster_rcone",0.4).asDouble();
    _cluster_etmin = detector.get("cluster_etmin",5.0).asDouble();
    _min_dr = detector.get("dr_min",0.2).asDouble();

    const Json::Value perf_objects  = root["Performance"];

    PProperties *particle_props;
    std::cout << "size of the list is " << perf_objects.size() << std::endl;
    for (Json::ValueIterator itr = perf_objects.begin(); itr != perf_objects.end(); itr++) {
//    for ( int index = 0; index < (int)perf_objects.size(); ++index ) { // Iterates over the sequence elements.

      particle_props = new PProperties;

      FastSim_ObjectReader(*itr,*particle_props);
      _detector_perf.push_back(particle_props);
    }

      /*
      type_str = phys_objects[index].get("type","jet").asString();

      min_pt = phys_objects[index].get("min_pt","10.0").asDouble();
      min_eta = phys_objects[index].get("min_eta","-4.0").asDouble();
      max_eta = phys_objects[index].get("max_eta","4.0").asDouble();

      if (type_str == "hadrons") { // hadrons only defined the minimum track pt
        _min_track_pt  = min_pt;  // GeV
        break;
      }

      level_str = phys_objects[index].get("level","loose").asString();
      if ((type_str == "electron") || (level_str == "loose")) {

        iso_cut = phys_objects[index].get("iso","0.4").asDouble();
        _min_ele_pt  = min_pt;  // GeV
        _minEt_isol_electron = iso_cut;
        _min_ele_eta = min_eta;
        _max_ele_eta = max_eta;
      }
      else if ((type_str == "muon") || (level_str == "loose")) {

        iso_cut = phys_objects[index].get("iso","0.4").asDouble();
        _minEt_isol_muon = iso_cut;
        _min_muon_pt  = min_pt;  // GeV
        _min_muon_eta = min_eta;
        _max_muon_eta = max_eta;
      }
      else if ((type_str == "photon") || (level_str == "loose")) {

        iso_cut = phys_objects[index].get("iso","0.4").asDouble();
        _minEt_isol_photon = iso_cut;
        _min_photon_pt  = min_pt;  // GeV
        _min_photon_eta = min_eta;
        _max_photon_eta = max_eta;
      }
      else if ((type_str == "jet") || (level_str == "loose")) {

        _min_jet_pt  = min_pt;  // GeV
        _min_jet_eta = min_eta;
        _max_jet_eta = max_eta;
      }
      else if ((type_str == "bjet") || (level_str == "loose")) {
        _min_bjet_pt  = min_pt;  // GeV
        _min_bjet_eta = min_eta;
        _max_bjet_eta = max_eta;
      }
      else if ((type_str == "tauhad") || (level_str == "loose")) {
        _min_tauhad_pt  = min_pt;  // GeV
        _min_tauhad_eta = min_eta;
        _max_tauhad_eta = max_eta;
      }
    }
    */

    _calo_neta = int(2*_calo_etamax/_calo_deta);
    _calo_nphi = int(2*3.2/_calo_dphi);



    // everything is fine
    return 0;

  }


  int FastSim::FastSim_ObjectReader(const Json::Value phys_objects,PProperties &particle_props) {

    // need to add warnings if the fields are missing, when we get the logger
     std::cout << "type " << std::endl;

    if (not phys_objects["type"].empty()) 
      particle_props._pid = phys_objects.get("type","13").asInt();

     std::cout << "level " << std::endl;
    if (not phys_objects["level"].empty()) 
      particle_props._level = phys_objects.get("level","loose").asString();

    if (not phys_objects["min_pt"].empty()) 
      particle_props._min_pt = phys_objects.get("min_pt","0.0").asDouble();
    std::cout << "min_pt " << std::endl;
    if (not phys_objects["min_eta"].empty()) 
      particle_props._min_eta = phys_objects.get("min_eta","-20000.0").asDouble();
    std::cout << "min_eta " << std::endl;
    if (not phys_objects["max_eta"].empty()) 
      particle_props._max_eta = phys_objects.get("max_eta","20000.0").asDouble();
    if (not phys_objects["iso"].empty()) 
      particle_props._iso = phys_objects.get("iso","20000").asDouble();
    std::cout << "iso " << std::endl;

    std::cout << phys_objects["acceptance_reco"].empty() << std::endl;
    if (not phys_objects["acceptance_reco"].empty()) {
      const Json::Value accept_objects = phys_objects["acceptance_reco"];

      particle_props._test_acceptance = true; 
      for (Json::ValueIterator itr = accept_objects.begin(); itr != accept_objects.end(); itr++) {

        if (((*itr)["name"].empty()) || ((*itr)["info"].empty())) {
          std::cout << " name or the histogram of acceptance is missing " << std::endl;
          continue;
        }
          
        std::cout << "reading the acceptance " << std::endl;
        if ("eta" == (*itr).get("name","eta").asString()) {

          std::cout << "eta " << (*itr).get("name","eta").asString() << std::endl;
          FastSim_AcceptanceItemReader(*itr,particle_props._eta);
        }
        else
        {
          std::cout << "not eta " <<  (*itr).get("name","eta").asString() << std::endl;
          FastSim_AcceptanceItemReader(*itr,particle_props._pt);
        }

      }
    }
    else {// the acceptance was null was it should be perfect
        std::cout << "acceptance was null" << std::endl;

     particle_props._test_acceptance = false;
    }


    return 0;
}


  int FastSim::FastSim_AcceptanceItemReader(const Json::Value reco_object,Acceptance &efficiency) {

    std::vector<double> bin_edgesx,bin_edgesy,histo_values;
    int ndim,nbinsx,nbinsy;
    
    efficiency._name = reco_object.get("name","eta").asString();
    FastSim_InfoReader(reco_object["info"],ndim,nbinsx,nbinsy, bin_edgesx,bin_edgesy,histo_values);

    return 0;
  }

  int FastSim::FastSim_InfoReader(const Json::Value histo_object, int &ndim, int &nbinsx, int &nbinsy,
      std::vector<double> &bin_edgesx, std::vector<double> &bin_edgesy, std::vector<double> &histo_values) {

    if ((histo_object["ndim"].empty()) || (histo_object["nbinx"].empty()) || (histo_object["axisx"].empty())) {
      std::cout << " missing info fields " << std::endl;
      return -1;
    }

    ndim = histo_object.get("ndim","1").asInt();
    nbinsx = histo_object.get("nbinx","50").asInt();

    // now get the bin edges
    const Json::Value axisx  = histo_object["axisx"];
    FastSim_ArrayReader(axisx,bin_edgesx);
    std::cout << " counted " << bin_edgesx.size() << " nbinsx " << nbinsx << std::endl;
    if ((int)bin_edgesx.size() != (nbinsx+1)) {
      std::cout << "error number of x bins and the number of y bin edges are wrong: " << bin_edgesx.size() << " " << nbinsx << std::endl;
      return -1;
    }

    if (ndim == 1) {
      const Json::Value hvalues  = histo_object["val"];
      FastSim_ArrayReader(hvalues,histo_values);
      if ((int)histo_values.size() != nbinsx) {
        std::cout << "1d error number of bins and number of values are different: " << histo_values.size() << " " << nbinsx << std::endl;
        return -1;
      }
    }
    else if (ndim == 2) {// we need to read the axisy

      nbinsy = histo_object.get("nbiny","50").asInt();

      const Json::Value axisy  = histo_object["axisy"];
      FastSim_ArrayReader(axisy,bin_edgesy);
      if ((int)bin_edgesy.size() != nbinsy+1) {
        std::cout << "error number of y bins and the number of y bin edges are wrong: " << bin_edgesy.size() << " " << nbinsx << std::endl;
        return -1;
      }
      const Json::Value hvalues  = histo_object["val"];
      FastSim_ArrayReader(hvalues,histo_values);
      if ((int)hvalues.size() != (nbinsx*nbinsy)) {
        std::cout << "2d error number of bins (binx*biny) and number of values are different: " << hvalues.size() << " " << nbinsx*nbinsy << std::endl;
        return -1;
      }
    }
    else {
      std::cout << " the number of dimesions is not supported " << ndim << std::endl;
    }

    return 0;
  }

  int FastSim::FastSim_ArrayReader(const Json::Value array, std::vector<double> &values ) {

    int size = array.size();

    for ( int index = 0; index < size; ++index ) { // Iterates over the binedges
      values.push_back(array[index].asDouble());
    }

    return 0;

  }
}
