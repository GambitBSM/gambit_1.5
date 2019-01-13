//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit event loop functions returning
///  detector simulations.
///
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Abram Krislock
///          (a.m.b.krislock@fys.uio.no)
///
///  \author Aldo Saavedra
///
///  \author Andy Buckley
///
///  \author Chris Rogan
///          (crogan@cern.ch)
///  \date 2014 Aug
///  \date 2015 May
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2015 Jul
///  \date 2018 Jan
///  \date 2019 Jan
///
///  \author Anders Kvellestad
///          (anders.kvellestad@fys.uio.no)
///  \date   2017 March
///  \date   2018 Jan
///  \date   2018 May
///
///  *********************************************

#include <memory>

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "gambit/ColliderBit/CMSEfficiencies.hpp"
#include "gambit/ColliderBit/detectors/BuckFast_definitions.hpp"

#include "HEPUtils/FastJet.h"

// #define COLLIDERBIT_DEBUG

namespace Gambit
{

  namespace ColliderBit
  {

    /// Get a BuckFast detector simulation
    template<typename EventT>
    BaseDetector<EventT>* getBuckFast(const str& detname,
                                      const MCLoopInfo& RunMC,
                                      bool use_effs,
                                      int iteration,
                                      const Options& runOptions)
    {
      static std::vector<bool> partonOnly;
      static std::vector<double> antiktR;

      // Where the real action is
      static std::unique_ptr<BuckFast<EventT>[]> bucky(new BuckFast<EventT>[omp_get_max_threads()]);
      int mine = omp_get_thread_num();

      if (iteration == BASE_INIT)
      {
        std::vector<bool> default_partonOnly(RunMC.collider_names.size(), false);
        partonOnly = runOptions.getValueOrDef<std::vector<bool> >(default_partonOnly, "partonOnly");
        CHECK_EQUAL_VECTOR_LENGTH(partonOnly,RunMC.collider_names)

        std::vector<double> default_antiktR(RunMC.collider_names.size(), 0.4);
        antiktR = runOptions.getValueOrDef<std::vector<double> >(default_antiktR, "antiktR");
        CHECK_EQUAL_VECTOR_LENGTH(antiktR,RunMC.collider_names)
      }

      if (iteration == START_SUBPROCESS)
      {
        // Each thread gets its own copy of the detector sim, so it is initialised *after* COLLIDER_INIT, within omp parallel.
        bucky[mine].init(partonOnly[RunMC.current_collider_index], antiktR[RunMC.current_collider_index]);
        // Assign detector functions
        if (detname == "ATLAS")
        {
          bucky[mine].smearElectronEnergy = &ATLAS::smearElectronEnergy;
          bucky[mine].smearMuonMomentum   = &ATLAS::smearMuonMomentum;
          bucky[mine].smearTaus           = &ATLAS::smearTaus;
          bucky[mine].smearJets           = &ATLAS::smearJets;
          if (use_effs)
          {
            bucky[mine].applyElectronEff  = &ATLAS::applyElectronEff;
            bucky[mine].applyMuonEff      = &ATLAS::applyMuonEff;
          }
        }
        else if (detname == "CMS")
        {
          bucky[mine].smearElectronEnergy = &CMS::smearElectronEnergy;
          bucky[mine].smearMuonMomentum   = &CMS::smearMuonMomentum;
          bucky[mine].smearTaus           = &CMS::smearTaus;
          bucky[mine].smearJets           = &CMS::smearJets;
          if (use_effs)
          {
            bucky[mine].applyElectronEff  = &CMS::applyElectronEff;
            bucky[mine].applyMuonEff      = &CMS::applyMuonEff;
          }
        }
        else if (detname == "Identity") { /* relax */ }
        else
        {
          ColliderBit_error().raise(LOCAL_INFO, "Unrecognised detector name.");
        }
      }

      // Paper-bag it
      return &bucky[mine];

    }

    /// Retrieve a BuckFast sim of CAT(EXPERIMENT,SUFFIX), that accepts type EVENT
    #define GET_BUCKFAST_AS_BASE_DETECTOR(NAME, EVENT, EXPERIMENT, SUFFIX)      \
    void NAME(BaseDetector<EVENT>* &result)                                     \
    {                                                                           \
      using namespace Pipes::NAME;                                              \
      result = getBuckFast<EVENT>(#EXPERIMENT, *Dep::RunMC,                     \
       IF_ELSE_EMPTY(SUFFIX,true,false), *Loop::iteration, *runOptions);        \
    }

  }

}
