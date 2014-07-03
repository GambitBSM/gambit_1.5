#include "FastSim.hpp"



#include <iostream>

using namespace std;

fast_sim::FastSim _my_fast_sim;


int FastSim_Init(int which);
int FastSim_ReadEvents( int &events);


int FastSim_Init(int which)
{
  switch(which) {
    default:
    case 0:_my_fast_sim.init(fast_sim::NOMINAL, 0);break; // added '0' in default spot so that gambit would compile
    case 1:_my_fast_sim.init(fast_sim::ATLAS, 0);break;
  }

  return 0;
}

int FastSim_ReadEvents( int &events) {

  vector<hep_simple_lib::Particle*> _electrons;
  vector<hep_simple_lib::Particle*> _muons;
  vector<hep_simple_lib::Particle*> _photons;
  vector<hep_simple_lib::Particle*> _nonprompt_leptons; //GM:  add so that gambit will compile
  vector<hep_simple_lib::Particle*> _bjets;
  vector<hep_simple_lib::Particle*> _tauhads;
  vector<hep_simple_lib::Particle*> _chargedhads;
  vector<hep_simple_lib::Particle*> _weakly_interacting; // stdm neutrinos, susy neutralinos

  hep_simple_lib::Particle *chosen;
  chosen = new hep_simple_lib::Particle(40.5, -32.6, 0.5, 51.992884131, 13);
  _electrons.push_back(chosen);
 
  hep_simple_lib::Event reco_event;
  _my_fast_sim.setParticles(_electrons, _muons, _photons, _nonprompt_leptons, _chargedhads, _bjets, _tauhads, _weakly_interacting);

  //cout << " do detector response " << std::endl;
  _my_fast_sim.doDetectorResponse();
  _my_fast_sim.getRecoEvent(reco_event);

  if ((int)reco_event.electrons().size() > 0)
    cout << " electron detected " << endl;
  else
    cout << " electron not detected " << endl;

  events = 1;

  for (size_t kk = 0; kk < _electrons.size();kk++) 
    delete _electrons[kk]; 
  _electrons.clear();



  return 0;
}


