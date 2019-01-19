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
                              int iteration,
                              const Options& runOptions)
    {
      static std::map<str, std::vector<str> > analyses;
      static bool first = true;

      if (iteration == BASE_INIT)
      {
        // Only run this once
        if (first)
        {
          // Loop over colliders
          for (auto& collider : RunMC.collider_names)
          {
            // Read analysis names from the yaml file
            if (runOptions.hasKey(collider))
            {
              YAML::Node colNode = runOptions.getValue<YAML::Node>(collider);
              analyses[collider] = colNode[collider].as<std::vector<str>>();
            }
            else analyses[collider] = std::vector<str>();

            // Check that the analysis names listed in the yaml file all correspond to actual ColliderBit analyses
            for (str& analysis_name : analyses[collider])
            {
              if (!checkAnalysis(analysis_name))
              {
                str errmsg = "The analysis " + analysis_name + " is not a known ColliderBit analysis.";
                ColliderBit_error().raise(LOCAL_INFO, errmsg);
              }
            }
          }
          first = false;
        }
      }

      if (analyses.empty()) return;

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
            result.init(analyses.at(RunMC.current_collider()));
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
       *(*Dep::HardScatteringSim), *Loop::iteration, *runOptions); \
    }

    GET_ANALYSIS_CONTAINER(getATLASAnalysisContainer, ATLAS)
    GET_ANALYSIS_CONTAINER(getATLASmultieffAnalysisContainer, ATLASmultieff)
    GET_ANALYSIS_CONTAINER(getCMSAnalysisContainer, CMS)
    GET_ANALYSIS_CONTAINER(getCMSmultieffAnalysisContainer, CMSmultieff)
    GET_ANALYSIS_CONTAINER(getIdentityAnalysisContainer, Identity)


  }
}
