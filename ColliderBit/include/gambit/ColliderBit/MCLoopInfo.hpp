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


#ifndef __MCLoopInfo_hpp__
#define __MCLoopInfo_hpp__

#include <vector>
#include "gambit/Utils/util_types.hpp"

namespace Gambit
{

  namespace ColliderBit
  {

    /// @brief Container for event loop status data and settings
    struct MCLoopInfo
    {
      /// Event generation has started
      bool event_generation_began;

      /// Maximum allowable number of failed events before MC loop is terminated
      int maxFailedEvents;

      /// Maximum allowed number of failed events has been reached and MC loop terminated
      bool exceeded_maxFailedEvents;

      /// The index of the current collider
      unsigned int current_collider_index;

      /// The name of the current collider
      str current_collider;

      /// The names of all colliders
      std::vector<str> collider_names;

      /// The random seed bases of all colliders
      std::map<str,int> seed_base;

      /// Number of events generated for all colliders
      std::map<str,int> event_count;

      /// The random seed base of the current collider
      int current_seed_base() const;

      /// Number of events generated for the current collider
      int current_event_count() const;

      /// Reinitialise
      void clear();
    };

  }
}



#endif /* defined __MCLoopInfo_hpp__ */
