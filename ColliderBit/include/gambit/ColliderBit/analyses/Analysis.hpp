//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Class for ColliderBit analyses.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Abram Krislock
///          (a.m.b.krislock@fys.uio.no)
///
///  \author Andy Buckley
///          (mostlikelytobefound@facebook.com)
///
///  \author Anders Kvellestad
///          (anders.kvellestad@fys.uio.no)
///  \date often
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2019 Feb
///
///  *********************************************

#pragma once

#include <string>
#include "HEPUtils/Event.h"
#include "gambit/ColliderBit/analyses/AnalysisData.hpp"

namespace Gambit
{
  typedef std::string str;

  namespace ColliderBit
  {

    /// A class for collider analyses within ColliderBit.
    class Analysis
    {

      public:

        /// @name Construction, Destruction, and Recycling:
        ///@{

        Analysis();
        virtual ~Analysis() { }

        /// @brief Public method to reset this instance for reuse, avoiding the need for "new" or "delete".
        /// @note This method calls _clear() to reset the base class variables, and then the overridden
        /// clear() to reset analysis-specific variables.
        /// @todo For v2.0: Simplify this reset scheme.
        void reset();

        /// @name Event analysis, event number, and cross section functions:
        ///@{
        /// Analyze the event (accessed by reference).
        void do_analysis(const HEPUtils::Event&);
        /// Analyze the event (accessed by pointer).
        void do_analysis(const HEPUtils::Event*);

        /// Return the total number of events seen so far.
        double num_events() const;
        /// Return the cross-section (in pb).
        double xsec() const;
        /// Return the cross-section error (in pb).
        double xsec_err() const;
        /// Return the cross-section relative error.
        double xsec_relerr() const;
        /// Return the cross-section per event seen (in pb).
        double xsec_per_event() const;
        /// Return the integrated luminosity (in inverse pb).
        void set_xsec(double, double);
        /// Set the integrated luminosity (in inverse pb).
        double luminosity() const;
        /// Set the cross-section and its error (in pb).
        void set_luminosity(double);
        /// Set the analysis name
        void set_analysis_name(str);
        /// Get the analysis name
        str analysis_name();

        /// Get the collection of SignalRegionData for likelihood computation.
        const AnalysisData& get_results();
        /// An overload of get_results() with some additional consistency checks.
        const AnalysisData& get_results(str&);
        /// Get a pointer to _results.
        const AnalysisData* get_results_ptr();
        /// Get a pointer to _results.
        const AnalysisData* get_results_ptr(str&);
        ///@}

        /// @name (Re-)initialization functions
        ///@{
        /// General init for any analysis of this type.
        virtual void init(const std::vector<std::string>&) {}
        /// General init for any collider of this type - no settings version.
        virtual void init() {}
        /// Scale by number of input events and xsec.
        virtual void scale(double factor=-1);
        ///@}

        /// @name Analysis combination operations
        ///@{
        /// An operator to do xsec-weighted combination of analysis runs.
        virtual void add(Analysis* other);
        /// Add cross-sections and errors for two different process types.
        void add_xsec(double xs, double xserr);
        /// Combine cross-sections and errors for the same process type, assuming uncorrelated errors.
        void improve_xsec(double xs, double xserr);
        ///@}

      protected:

        /// @brief Overloadable method to reset the analysis-specific variables.
        //
        /// @note Internally, reset() is called clear() -- we avoid this
        /// externally because of confusion with std::vector::clear(), esp. on
        /// AnalysisContainer.
        virtual void clear() = 0;

        /// @name Collection functions
        ///@{
        /// Analyze the event (accessed by pointer).
        /// @note Needs to be called from Derived::analyze().
        virtual void analyze(const HEPUtils::Event*);
        /// Add the given result to the internal results list.
        void add_result(const SignalRegionData& sr);
        /// Set the covariance matrix, expressing SR correlations
        void set_covariance(const Eigen::MatrixXd& srcov);
        /// A convenience function for setting the SR covariance from a nested vector/initialiser list
        void set_covariance(const std::vector<std::vector<double>>&);
        /// Gather together the info for likelihood calculation.
        virtual void collect_results() = 0;
        ///@}

      private:

        double _ntot, _xsec, _xsecerr, _luminosity;
        bool _xsec_is_set, _luminosity_is_set, _is_scaled, _needs_collection;
        AnalysisData _results;
        std::string _analysis_name;

        /// @brief Reset the private base class variables.
        /// @todo For v2.0: Avoid this 'duplication' of reset/clear methods.
        void _clear();
        ///@}

    };


    /// For analysis factory function definition
    #define DEFINE_ANALYSIS_FACTORY(ANAME)                                     \
      Analysis* create_Analysis_ ## ANAME()                                    \
      {                                                                        \
        return new Analysis_ ## ANAME();                                       \
      }                                                                        \
      std::string getDetector_ ## ANAME()                                      \
      {                                                                        \
        return std::string(Analysis_ ## ANAME::detector);                      \
      }


  }
}
