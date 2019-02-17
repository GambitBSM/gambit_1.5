#pragma once
//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  The BaseDetector class.

#include <string>
#include <vector>
#include <exception>
#include <memory>

#include "gambit/Elements/shared_types.hpp"

namespace Gambit
{

  namespace ColliderBit
  {

    /// An abstract base class for detector simulators within ColliderBit.
    template<typename EventT>
    class BaseDetector
    {

      public:

        /// Constructor
        BaseDetector() {}
        /// Destructor
        virtual ~BaseDetector() {}
        /// Reset this instance for reuse, avoiding the need for "new" or "delete".
        virtual void clear() {}

        /// Perform the actual simulation on the next collider event.
        virtual void processEvent(const EventT&, HEPUtils::Event&) const = 0;

        /// @name (Re-)Initialization functions
        ///@{
        /// Settings parsing and initialization for each sub-class.
        virtual void init(const std::vector<std::string>&) {}
        /// General init for any collider of this type - no settings version.
        virtual void init() {};
        ///@}

    };

  }
}
