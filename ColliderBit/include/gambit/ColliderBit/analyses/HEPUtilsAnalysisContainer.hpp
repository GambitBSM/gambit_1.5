#pragma once
#include <string>
#include <stdexcept>
#include <vector>

// Forward declarations, to avoid header-chaining into CB_types.hpp
namespace HEPUtils { class Event; }
namespace Gambit {
  namespace ColliderBit {
    template <typename EventT>
    class BaseAnalysis;
    using HEPUtilsAnalysis = BaseAnalysis<HEPUtils::Event>;
  }
}


namespace Gambit {
  namespace ColliderBit {


    /// Create a new analysis based on a name string
    /// @note The caller is responsible for deleting the returned analysis object.
    /// @todo Move to a separate file
    HEPUtilsAnalysis* mkAnalysis(const std::string& name);


    /// More sophisticated than just std::vector<HEPUtilsAnalysis*>
    struct HEPUtilsAnalysisContainer
    {
      /// @todo Really need to make our minds up about wtf this type is to be called -- it has 3 different names.
      ///       Just pick one and delete the typedefs, they make things messier and harder, not easier.
      typedef HEPUtilsAnalysis Analysis;

      std::vector<HEPUtilsAnalysis*> analyses;
      bool ready;

      /// @name Construction, destruction, and recycling:
      //@{
      HEPUtilsAnalysisContainer() : ready(false) { }
      ~HEPUtilsAnalysisContainer() { clear(); }
      /// Delete and clear the analyses contained within this instance.
      void clear();
      //@}

      // /// @name Access contents info:
      // //@{
      // /// List of contained analysis names
      // /// Number of contained analyses
      // //@}

      /// @name (Re-)Initialization functions:
      //@{
      /// Initialize analyses by their names.
      void init(const std::vector<std::string>& analysisNames);
      /// Re-initialize/reset already-instantiated analyses.
      void reset();
      //@}

      /// @name Event analysis
      //@{
      /// Analyze an event.
      /// @todo Analyze is const?! What about internal counter updates?
      void analyze(const HEPUtils::Event&) const;
      //@}

      /// @name Analysis post-processing (combination, scaling, etc.):
      //@{
      /// Add cross-sections and errors for two different process types
      void add_xsec(double, double);
      /// Add cross-sections and errors for two different process types
      void add_xsec(const HEPUtilsAnalysisContainer& e) { add_xsec(&e); }
      /// Add cross-sections and errors for two different process types
      void add_xsec(const HEPUtilsAnalysisContainer*);
      /// Combine cross-sections and errors for the same process type
      void improve_xsec(double, double);
      /// Combine cross-sections and errors for the same process type
      void improve_xsec(const HEPUtilsAnalysisContainer& e) { improve_xsec(&e); }
      /// Combine cross-sections and errors for the same process type
      void improve_xsec(const HEPUtilsAnalysisContainer*);
      /// Add the results of all analyses from this instance to the given one.
      void add(const HEPUtilsAnalysisContainer& e) { add(&e); }
      /// Add the results of all analyses from this instance to the given one.
      void add(const HEPUtilsAnalysisContainer*);
      /// Set cross-sections and errors for each analysis.
      /// @note If factor is negative (the default), each analysis uses its own nEvents, xsec, and luminosity to scale
      void scale(double factor=-1);
      //@}

    };


  }
}
