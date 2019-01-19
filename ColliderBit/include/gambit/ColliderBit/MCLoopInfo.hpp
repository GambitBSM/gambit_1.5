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


#pragma once

#include <vector>
#include "gambit/Utils/util_types.hpp"
#include "gambit/ColliderBit/MC_convergence.hpp"

namespace Gambit
{

  namespace ColliderBit
  {

    /// @brief Container for event loop status data and settings
    struct MCLoopInfo
    {
      /// Event generation has started
      bool event_generation_began;

      /// Maximum allowed number of failed events has been reached and MC loop terminated
      bool exceeded_maxFailedEvents;

      /// The names of all colliders
      std::vector<str> collider_names;

      /// Maximum allowable number of failed events before MC loop is terminated for each collider
      std::map<str,int> maxFailedEvents;

      /// The random seed base of each collider
      std::map<str,int> seed_base;

      /// Number of events generated for each collider
      std::map<str,int> event_count;

      /// Convergence options for each collider
      std::map<str,convergence_settings> convergence_options;

      /// Set the current collider
      void set_current_collider(str&);

      /// Get the current collider
      const str& current_collider() const;

      /// Get/set maximum allowable number of failed events before MC loop is terminated for the current collider
      int& current_maxFailedEvents() const;

      /// Get/set the random seed base of the current collider
      int& current_seed_base() const;

      /// Get/set the number of events generated for the current collider
      int& current_event_count() const;

      /// Reinitialise
      void clear();

      private:

        /// The name of the current collider
        str _current_collider;

        /// Iterator to the current maxFailedEvents
        std::map<str,int>::iterator _current_maxFailedEvents_it;

        /// Iterator to the current seed base
        std::map<str,int>::iterator _current_seed_base_it;

        /// Iterator to the current event count
        std::map<str,int>::iterator _current_event_count_it;

    };

  }
}

