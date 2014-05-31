#include "simple_hep_lib/Particle.hpp"

#include <gsl/gsl_rng.h>
#include <time.h>


namespace fast_sim {


// the base class of the Detector Response
class DetectorResponse
{
  protected:

    double _muon_resolution;
    double _electron_resolution;
    double _photon_resolution;
    double _jet_resolution;


    gsl_rng *_random_num;

  public:
    DetectorResponse();
    ~DetectorResponse();
    virtual void MuonResponse(hep_simple_lib::Particle& muon) { };
    virtual void PhotonResponse(hep_simple_lib::Particle& ph) {};
    virtual void ElectronResponse(hep_simple_lib::Particle& ele) {};
    virtual void JetResponse(hep_simple_lib::Particle& jet) {};
};

class ATLAS_Simple_Response: public DetectorResponse
{

  public:
    ATLAS_Simple_Response();
    void MuonResponse(hep_simple_lib::Particle& muon);
    void PhotonResponse(hep_simple_lib::Particle& ph);
    void ElectronResponse(hep_simple_lib::Particle& ele);
    void JetResponse(hep_simple_lib::Particle& jet);
};


//RESHAD============================================

/*
double RESHAD(double e, double eta, double Caloth){

  if(SMEAR ==1){
  double distr = rndm.Gaus(0,1);

  double res = -999.99;

  //Check which region of the detector.
  if(eta < Caloth){
    res = 0.5;
  }else if(eta > Caloth){
    res = 1.0;
  }

  double sig = distr * res *1/sqrt(e);

  return sig;
  }
  else return 0;
}

//}
*/

}
