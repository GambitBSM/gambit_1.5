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

        /// Public method to reset this instance for reuse, avoiding the need for "new" or "delete".
        void reset();

        /// @name Event analysis, event number, and cross section functions:
        ///@{
        /// Analyze the event (accessed by reference).
        void analyze(const HEPUtils::Event&);
        /// Analyze the event (accessed by pointer).
        void analyze(const HEPUtils::Event*);

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
        /// Set the cross-section and its error (in pb).
        void set_xsec(double, double);
        /// Return the integrated luminosity (in inverse pb).
        double luminosity() const;
        /// Set the integrated luminosity (in inverse pb).
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
        /// General init for any analysis of this type - no settings version.
        virtual void init() {}
        /// Scale by number of input events and xsec.
        void scale(double factor=-1);
        ///@}

        /// @name Analysis combination operations
        ///@{
        /// Add the results of another analysis to this one. Argument is not const, because the other needs to be able to gather its results if necessary.
        void add(Analysis* other);
        /// Add the analysis-specific variables of another analysis to this one.
        virtual void combine(const Analysis* other) = 0;
        /// Add cross-sections and errors for two different process types.
        void add_xsec(double xs, double xserr);
        /// Combine cross-sections and errors for the same process type, assuming uncorrelated errors.
        void improve_xsec(double xs, double xserr);
        ///@}

      protected:

        /// Reset the analysis-specific variables.
        virtual void analysis_specific_reset() = 0;

        /// @name Collection functions
        ///@{
        /// Run the analysis.
        virtual void run(const HEPUtils::Event*) = 0;
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

        double _ntot;
        double _xsec;
        double _xsecerr;
        double _luminosity;
        bool _xsec_is_set;
        bool _luminosity_is_set;
        bool _is_scaled;
        bool _needs_collection;
        AnalysisData _results;
        std::string _analysis_name;

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
