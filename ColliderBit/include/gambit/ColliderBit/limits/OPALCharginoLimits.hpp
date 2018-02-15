#pragma once
#include "gambit/ColliderBit/limits/BaseLimitContainer.hpp"

namespace Gambit
{
  namespace ColliderBit 
  {
    /// @brief A class to contain the limit data from OPAL, hep-ex/0210043, figure 5a (in colour)
    class OPALCharginoMassSplitLimitAt208GeV : public BaseLimitContainer {

      public:

        /// @read off the csv file containting the data
        std::vector<std::pair<double,double> > dataForLimit(double);
        /// @brief Check to see if the point is within the exclusion region
        bool isWithinExclusionRegion(double x, double y, double) const;

        /// @name Construction, initializing with all necessary data from the plot
        OPALCharginoMassSplitLimitAt208GeV();
    };

  }
}
