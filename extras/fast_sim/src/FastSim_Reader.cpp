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
#include "FastSim_Logger.hpp"
#include "json/json.h"

#include <algorithm>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>


extern logging::logger log_inst(0);

namespace fast_sim {

  int FastSim::FastSim_Reader(std::string init_filename) {
    // this functuon reads the jason file and initialises the class

    std::ifstream fastsim_datacard(init_filename.c_str());
    if (not fastsim_datacard.good()) {

      LOG_ERR("JSonReader Could not open file:",init_filename);
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
      LOG_WARN("JSonReader Failed to parse configuration:",reader.getFormattedErrorMessages());
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
    std::string temp = detector.get("jet_reco","fastjet").asString();
    if (temp == "fastjet")
      _fastjet = true;
    else
      _fastjet = false;

    const Json::Value perf_objects  = root["Performance"];

    PProperties *particle_props;

    LOG_DEBUG1("JSonReader Number of Performance objects:",perf_objects.size());
    for (Json::ValueIterator itr = perf_objects.begin(); itr != perf_objects.end(); itr++) {

      particle_props = new PProperties;

      FastSim_ObjectReader(*itr,*particle_props);
      _detector_perf.push_back(particle_props);
    }

    // a hack for now
    _calo_neta = int(2*_calo_etamax/_calo_deta);
    _calo_nphi = int(2*3.2/_calo_dphi);

    // everything is fine
    return 0;

  }


  int FastSim::FastSim_ObjectReader(const Json::Value phys_objects,PProperties &particle_props) {

    int val; // return value

    // need to add warnings if the fields are missing, when we get the logger
    if (not phys_objects["type"].empty()) 
      particle_props._pid = phys_objects.get("type","13").asInt();
    LOG_DEBUG1("JSonReader Particle Type ",particle_props._pid);

    if (not phys_objects["level"].empty()) 
      particle_props._level = phys_objects.get("level","loose").asString();
    LOG_DEBUG1("JSonReader Particle Category ",particle_props._level);

    if (not phys_objects["min_pt"].empty()) 
      particle_props._min_pt = phys_objects.get("min_pt","0.0").asDouble();
    LOG_DEBUG1("JSonReader Particle Min_Pt Cut ",particle_props._min_pt);

    if (not phys_objects["min_eta"].empty()) 
      particle_props._min_eta = phys_objects.get("min_eta","-20000.0").asDouble();
    LOG_DEBUG1("JSonReader Particle Min_Eta Cut ",particle_props._min_eta);

    if (not phys_objects["max_eta"].empty()) 
      particle_props._max_eta = phys_objects.get("max_eta","20000.0").asDouble();
    LOG_DEBUG1("JSonReader Particle Max_Eta Cut ",particle_props._max_eta);

    if (not phys_objects["iso"].empty()) 
      particle_props._iso = phys_objects.get("iso","20000").asDouble();
    LOG_DEBUG1("JSonReader Particle Isolation Cut ",particle_props._iso);

    LOG_DEBUG1("JSonReader Particle Acceptance Reco ");
    if (not phys_objects["acceptance_reco"].empty()) {
      const Json::Value accept_objects = phys_objects["acceptance_reco"];

      particle_props._test_acceptance = true; 
      for (Json::ValueIterator itr = accept_objects.begin(); itr != accept_objects.end(); itr++) {

        if (((*itr)["name"].empty()) || ((*itr)["info"].empty())) {
          LOG_WARN("JSonReader Particle name or histogram of acceptance is missing ");
          continue;
        }

        Acceptance new_acceptance;

        std::string accept_name = (*itr).get("name","xxx").asString();
        if ("eta" == accept_name) {
          LOG_DEBUG1("JSonReader Particle Acceptance Name ", (*itr).get("name","eta").asString());
          val=FastSim_AcceptanceItemReader(*itr,new_acceptance);
          if (val != 0)
            return val;
          new_acceptance._init = true;
          new_acceptance._type = ETA;
          particle_props._response.push_back(new_acceptance);
          LOG_DEBUG1("JSonReader Particle Acceptance Ok ", new_acceptance._name);
        }
        else if ("pt" == accept_name) 
        {
          LOG_DEBUG1("JSonReader Particle Acceptance Name ", (*itr).get("name","pt").asString());
          val=FastSim_AcceptanceItemReader(*itr,new_acceptance);
          if (val != 0)
            return val;
          new_acceptance._init = true;
          new_acceptance._type = PT;
          particle_props._response.push_back(new_acceptance);
          LOG_DEBUG1("JSonReader Particle Acceptance Ok ", new_acceptance._name);
        }
        else if ("pt_eta" ==  accept_name) 
        {
          LOG_DEBUG1("JSonReader Particle Acceptance Name ", (*itr).get("name","pt_eta").asString());
          val=FastSim_AcceptanceItemReader(*itr,new_acceptance);
          if (val != 0)
            return val;
          new_acceptance._init = true;
          new_acceptance._type = PT_ETA;
          particle_props._response.push_back(new_acceptance);
          LOG_DEBUG1("JSonReader Particle Acceptance Read ", new_acceptance._name);
        }
        else if ("iso" == accept_name) 
        {
          LOG_DEBUG1("JSonReader Particle Acceptance Name ", (*itr).get("name","iso").asString());
          val=FastSim_AcceptanceItemReader(*itr,new_acceptance);
          if (val != 0)
            return val;
          new_acceptance._init = true;
          new_acceptance._type = ISO;
          particle_props._response.push_back(new_acceptance);
          LOG_DEBUG1("JSonReader Particle Acceptance Read ", new_acceptance._name);
        }
        else {

          LOG_ERR("JSonReader Response Object not defined ", accept_name);
          return -1;
        }
      }
    }
    else {// the acceptance was null was it should be perfect

      LOG_WARN("JSonReader Acceptance was Null, efficiency is set to 1 for the complete range ");
     particle_props._test_acceptance = false;
    }


    return 0;
}


  int FastSim::FastSim_AcceptanceItemReader(const Json::Value reco_object,Acceptance &efficiency) {

    std::vector<double> bin_edgesx,bin_edgesy,histo_values;
    int ndim,nbinsx,nbinsy,error;
    
    efficiency._name = reco_object.get("name","eta").asString();

    LOG_DEBUG1("JSonReader Acceptance name ",efficiency._name);
    error = FastSim_InfoReader(reco_object["info"],ndim,nbinsx,nbinsy, bin_edgesx,bin_edgesy,histo_values);
    if (error != 0) // check whether that is something went wrong when reading the Info
      return error;

    efficiency._ndim = ndim;
    efficiency._nbinsx = nbinsx;

    for (size_t kk=0;kk<bin_edgesx.size();kk++)
      efficiency._bin_edges_x.push_back(bin_edgesx[kk]);

    if (ndim == 2) {

      efficiency._nbinsy = nbinsy;
      for (size_t kk=0;kk<bin_edgesy.size();kk++)
        efficiency._bin_edges_y.push_back(bin_edgesy[kk]);
    }

    for (size_t kk=0;kk<histo_values.size();kk++)
      efficiency._binvals.push_back(histo_values[kk]);

    return 0;
  }

/*
  int FastSim::FastSim_AcceptanceItemReader(const Json::Value reco_object,Acceptance_2D &efficiency) {

    std::vector<double> bin_edgesx,bin_edgesy,histo_values;
    int ndim,nbinsx,nbinsy;
    
    efficiency._name = reco_object.get("name","eta").asString();
    FastSim_InfoReader(reco_object["info"],ndim,nbinsx,nbinsy, bin_edgesx,bin_edgesy,histo_values);

    efficiency._ndim = ndim;
    efficiency._nbinsx = nbinsx;
    efficiency._nbinsy = nbinsy;

    for (size_t kk=0;kk<bin_edgesx.size();kk++)
      efficiency._bin_edges_x.push_back(bin_edgesx[kk]);

    for (size_t kk=0;kk<bin_edgesy.size();kk++)
      efficiency._bin_edges_y.push_back(bin_edgesy[kk]);

    for (size_t kk=0;kk<histo_values.size();kk++)
      efficiency._binvals.push_back(histo_values[kk]);

    std::cout << "Acceptance 2D name " << efficiency._name << std::endl;
    for (size_t kk=0;kk<histo_values.size();kk++)
      std::cout << " bin edges " << efficiency._bin_edges_x[kk] << std::endl;;
    for (size_t kk=0;kk<histo_values.size();kk++)
      std::cout << " bin edges " << efficiency._bin_edges_x[kk] << std::endl;;

    return 0;
  }
*/



  int FastSim::FastSim_InfoReader(const Json::Value histo_object, int &ndim, int &nbinsx, int &nbinsy,
      std::vector<double> &bin_edgesx, std::vector<double> &bin_edgesy, std::vector<double> &histo_values) {

    if ((histo_object["ndim"].empty()) || (histo_object["nbinx"].empty()) || (histo_object["axisx"].empty())) {
      LOG_ERR("JSonReader missing info fields - either ndim, nbinx axisx nbiny axisy are not defined in file ");
      return -1;
    }

    ndim = histo_object.get("ndim","1").asInt();
    LOG_DEBUG1("JSonReader - Number of dimensions is",ndim);
    switch(ndim) {
      case 1: { // read only the x-axis
                nbinsx = histo_object.get("nbinx","50").asInt();
                const Json::Value axisx  = histo_object["axisx"];
                FastSim_ArrayReader(axisx,bin_edgesx);
                if ((int)bin_edgesx.size() != (nbinsx+1)) {
                  LOG_ERR("JSonReader Error the number of bins specified in the x-axis and the number of edges are incompatible",bin_edgesx.size(),nbinsx);
                  return -1;
                }

                const Json::Value hvalues  = histo_object["val"];
                FastSim_ArrayReader(hvalues,histo_values);
                if ((int)histo_values.size() != nbinsx) {
                  LOG_ERR("JSonReader Error the number of bins and the number bin values are incompatible",histo_values.size(),nbinsx);
                  return -1;
                }
              }
              break;
      case 2: {// we need to read the y-axis as well as the x-axis

                nbinsx = histo_object.get("nbinx","50").asInt();
                const Json::Value axisx  = histo_object["axisx"];
                FastSim_ArrayReader(axisx,bin_edgesx);
                if ((int)bin_edgesx.size() != (nbinsx+1)) {
                  LOG_ERR("JSonReader Error the number of bins specified in the x-axis and the number of edges are incompatible",bin_edgesx.size(),nbinsx);
                  return -1;
                }

                nbinsy = histo_object.get("nbiny","50").asInt();
                const Json::Value axisy  = histo_object["axisy"];
                FastSim_ArrayReader(axisy,bin_edgesy);
                if ((int)bin_edgesy.size() != nbinsy+1) {
                  LOG_ERR("JSonReader Error the number of bins specified in the y-axis and the number bin values are incompatible",histo_values.size(),nbinsx);
                  return -1;
                }
                const Json::Value hvalues  = histo_object["val"];
                FastSim_ArrayReader(hvalues,histo_values);
                if ((int)hvalues.size() != (nbinsx*nbinsy)) {
                  LOG_ERR("JSonReader Error the number of bins (nbinx*nbiny)",nbinsx*nbinsy,"and the number bin values are incompatible",histo_values.size());
                  return -1;
                }
              }
              break;
      default:
              LOG_ERR("JSonReader Error the number of dimensions is not supported",ndim);
              return -1;
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
