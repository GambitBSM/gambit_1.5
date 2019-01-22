//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Initialisation functions for ColliderBit
///  analyses.
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

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"

// #define COLLIDERBIT_DEBUG

namespace Gambit
{

  namespace ColliderBit
  {

    /// Retrieve an analysis container for a specific detector
    void getAnalysisContainer(HEPUtilsAnalysisContainer& result,
                              const str& detname,
                              const MCLoopInfo& RunMC,
                              const BaseCollider& HardScatteringSim,
                              int iteration)
    {
      if (RunMC.analyses.empty() or iteration == BASE_INIT) return;
      if (not RunMC.current_analyses_exist_for(detname)) return;

      if (iteration == START_SUBPROCESS)
      {
        // Register analysis container
        result.register_thread(detname+"AnalysisContainer");

        // Set current collider
        result.set_current_collider(RunMC.current_collider());

        // Initialize analysis container or reset all the contained analyses
        if (!result.has_analyses())
        {
          try
          {
            result.init(RunMC.current_analyses_for(detname));
          }
          catch (std::runtime_error& e)
          {
            piped_errors.request(LOCAL_INFO, e.what());
          }
        }
        else result.reset();
      }

      if (iteration == END_SUBPROCESS && RunMC.event_generation_began)
      {
        if (not RunMC.exceeded_maxFailedEvents)
        {
          const double xs_fb = HardScatteringSim.xsec_pb() * 1000.;
          const double xserr_fb = HardScatteringSim.xsecErr_pb() * 1000.;
          result.add_xsec(xs_fb, xserr_fb);

          #ifdef COLLIDERBIT_DEBUG
          cout << debug_prefix() << "xs_fb = " << xs_fb << " +/- " << xserr_fb << endl;
          #endif
        }
      }

      if (iteration == COLLIDER_FINALIZE)
      {
        result.collect_and_add_signal();
        result.collect_and_improve_xsec();
        result.scale();
      }

    }

    /// Retrieve a container for analyses with EXPERIMENT
    #define GET_ANALYSIS_CONTAINER(NAME, EXPERIMENT)               \
    void NAME(HEPUtilsAnalysisContainer& result)                   \
    {                                                              \
      using namespace Pipes::NAME;                                 \
      getAnalysisContainer(result, #EXPERIMENT, *Dep::RunMC,       \
       *(*Dep::HardScatteringSim), *Loop::iteration);              \
    }

    GET_ANALYSIS_CONTAINER(getATLASAnalysisContainer, ATLAS)
    GET_ANALYSIS_CONTAINER(getCMSAnalysisContainer, CMS)
    GET_ANALYSIS_CONTAINER(getIdentityAnalysisContainer, Identity)


  }
}
