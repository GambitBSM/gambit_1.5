//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Type for holding event loop information.
///
///  *********************************************
///
///  Authors (add name if you modify):
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2019 Jan
///
///  *********************************************


#include "gambit/ColliderBit/MCLoopInfo.hpp"
#include "gambit/Utils/standalone_error_handlers.hpp"


namespace Gambit
{

  namespace ColliderBit
  {

    int MCLoopInfo::current_seed_base() const
    {
      auto it = seed_base.find(current_collider);
      if (it == seed_base.end())
      {
        str msg = "Current collider \"" + current_collider + "\" not found in MCLoopInfo::seed_base map.";
        utils_error().raise(LOCAL_INFO, msg);
      }
      return it->second;
    }

    /// Number of events generated for the current collider
    int MCLoopInfo::current_event_count() const
    {
      auto it = event_count.find(current_collider);
      if (it == event_count.end())
      {
        str msg = "Current collider \"" + current_collider + "\" not found in MCLoopInfo::current_event_count map.";
        utils_error().raise(LOCAL_INFO, msg);
      }
      return it->second;
    }

    /// Reinitialise
    void MCLoopInfo::clear()
    {
      event_generation_began = false;
      maxFailedEvents = -1;
      exceeded_maxFailedEvents = false;
      current_collider_index = 0;
      current_collider = "";
      collider_names.clear();
      seed_base.clear();
      event_count.clear();
    }

  }

}
