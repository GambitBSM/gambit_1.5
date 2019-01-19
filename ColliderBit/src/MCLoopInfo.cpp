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

    /// Reinitialise
    void MCLoopInfo::clear()
    {
      event_generation_began = false;
      exceeded_maxFailedEvents = false;
      collider_names.clear();
      _current_collider = "";
      maxFailedEvents.clear();
      seed_base.clear();
      event_count.clear();
    }

    /// Set the current collider
    void MCLoopInfo::set_current_collider(str& col)
    {
      // Save the current collider
      _current_collider = col;

      // Save the current maxFailedEvents
      auto it = maxFailedEvents.find(_current_collider);
      if (it == maxFailedEvents.end())
      {
        str msg = "Current collider \"" + _current_collider + "\" not found in MCLoopInfo::maxFailedEvents map.";
        utils_error().raise(LOCAL_INFO, msg);
      }
      _current_maxFailedEvents_it = it;

      // Save the current seed base
      it = seed_base.find(_current_collider);
      if (it == seed_base.end())
      {
        str msg = "Current collider \"" + _current_collider + "\" not found in MCLoopInfo::seed_base map.";
        utils_error().raise(LOCAL_INFO, msg);
      }
      _current_seed_base_it = it;

      // Save the number of events generated for the current collider
      it = event_count.find(_current_collider);
      if (it == event_count.end())
      {
        str msg = "Current collider \"" + _current_collider + "\" not found in MCLoopInfo::current_event_count map.";
        utils_error().raise(LOCAL_INFO, msg);
      }
      _current_event_count_it = it;

    }

    const str& MCLoopInfo::current_collider() const { return _current_collider; }

    int& MCLoopInfo::current_maxFailedEvents() const { return _current_maxFailedEvents_it->second; }

    int& MCLoopInfo::current_seed_base() const { return _current_seed_base_it->second; }

    int& MCLoopInfo::current_event_count() const { return _current_event_count_it->second; }

  }

}
