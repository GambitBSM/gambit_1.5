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

#include <vector>

#include "HEPUtils/MathUtils.h"
#include "HEPUtils/Event.h"

#include "gambit/ColliderBit/analyses/Analysis.hpp"

namespace Gambit
{
  namespace ColliderBit
  {

    Analysis::Analysis() : _ntot(0)
                         , _xsec(0)
                         , _xsecerr(0)
                         , _luminosity(0)
                         , _xsec_is_set(false)
                         , _luminosity_is_set(false)
                         , _is_scaled(false)
                         , _needs_collection(true)
                         { }

    /// Public method to reset this instance for reuse, avoiding the need for "new" or "delete".
    void Analysis::reset()
    {
      _ntot = 0;
      _xsec = 0;
      _xsecerr = 0;
      _xsec_is_set = false;
      _is_scaled = false;
      _needs_collection = true;
      _results.clear();
      analysis_specific_reset();
    }

    /// Analyze the event (accessed by reference).
    void Analysis::analyze(const HEPUtils::Event& e) { analyze(&e); }

    /// Analyze the event (accessed by pointer).
    void Analysis::analyze(const HEPUtils::Event* e)
    {
      _needs_collection = true;
      _ntot += 1;
      run(e);
    }

    /// Return the total number of events seen so far.
    double Analysis::num_events() const { return _ntot; }

    /// Return the cross-section (in pb).
    double Analysis::xsec() const { return _xsec; }

    /// Return the cross-section error (in pb).
    double Analysis::xsec_err() const { return _xsecerr; }

    /// Return the cross-section relative error.
    double Analysis::xsec_relerr() const { return xsec() > 0 ? xsec_err()/xsec() : 0; }

    /// Return the cross-section per event seen (in pb).
    double Analysis::xsec_per_event() const { return (xsec() >= 0 && num_events() > 0) ? xsec()/num_events() : 0; }

    /// Return the integrated luminosity (in inverse pb).
    void Analysis::set_xsec(double xs, double xserr) { _xsec_is_set = true; _xsec = xs; _xsecerr = xserr; }

    /// Set the integrated luminosity (in inverse pb).
    double Analysis::luminosity() const { return _luminosity; }

    /// Set the cross-section and its error (in pb).
    void Analysis::set_luminosity(double lumi) { _luminosity_is_set = true; _luminosity = lumi; }

    /// Set the analysis name
    void Analysis::set_analysis_name(str aname)
    {
      _analysis_name = aname;
      _results.analysis_name = _analysis_name;
    }

    /// Get the analysis name
    str Analysis::analysis_name() { return _analysis_name; }

    /// Get the collection of SignalRegionData for likelihood computation.
    const AnalysisData& Analysis::get_results()
    {
      if (_needs_collection)
      {
        collect_results();
        _needs_collection = false;
      }
      return _results;
    }

    /// An overload of get_results() with some additional consistency checks.
    const AnalysisData& Analysis::get_results(str& warning)
    {
      warning = "";
      if (not _xsec_is_set)
        warning += "Cross section has not been set for analysis " + _analysis_name + "." ;
      if (not _luminosity_is_set)
        warning += "Luminosity has not been set for analysis " + _analysis_name + ".";
      if (not _is_scaled)
        warning += "Results have not been scaled for analysis " + _analysis_name + ".";
      if (_ntot < 1)
        warning += "No events have been analyzed for analysis " + _analysis_name + ".";

      /// @todo We need to shift the 'analysis_name' property from class SignalRegionData
      ///       to this class. Then we can add the class name to this error message.
      // warning = "Ooops! In analysis " + analysis_name + ": " + warning

      return get_results();
    }

    /// Get a pointer to _results.
    const AnalysisData* Analysis::get_results_ptr()
    {
      return &get_results();
    }

    /// Get a pointer to _results.
    const AnalysisData* Analysis::get_results_ptr(str& warning)
    {
      return &get_results(warning);
    }

    /// Add the given result to the internal results list.
    void Analysis::add_result(const SignalRegionData& sr) { _results.add(sr); }

    /// Set the covariance matrix, expressing SR correlations
    void Analysis::set_covariance(const Eigen::MatrixXd& srcov) { _results.srcov = srcov; }

    /// A convenience function for setting the SR covariance from a nested vector/initialiser list
    void Analysis::set_covariance(const std::vector<std::vector<double>>& srcov)
    {
      Eigen::MatrixXd cov(srcov.size(), srcov.front().size());
      for (size_t i = 0; i < srcov.size(); ++i)
      {
        for (size_t j = 0; j < srcov.front().size(); ++j)
        {
          cov(i,j) = srcov[i][j];
        }
      }
      set_covariance(cov);
    }

    /// Scale by number of input events and xsec.
    void Analysis::scale(double factor)
    {
      if (factor < 0)
      {
        factor = (num_events() == 0 ? 0 : (luminosity() * xsec()) / num_events());
      }
      assert(factor >= 0);
      for (SignalRegionData& sr : _results)
      {
        sr.n_signal_at_lumi = factor * sr.n_signal;
      }
      _is_scaled = true;
    }

    /// Add the results of another analysis to this one. Argument is not const, because the other needs to be able to gather its results if necessary.
    void Analysis::add(Analysis* other)
    {
      if (_results.empty()) collect_results();
      const AnalysisData otherResults = other->get_results();
      /// @todo Access by name, including merging disjoint region sets?
      assert(otherResults.size() == _results.size());
      for (size_t i = 0; i < _results.size(); ++i)
      {
        _results[i].n_signal += otherResults[i].n_signal;
      }
      _ntot += other->num_events();
      combine(other);
    }

    /// Add cross-sections and errors for two different process types.
    void Analysis::add_xsec(double xs, double xserr)
    {
      if (xs > 0)
      {
        if (xsec() <= 0)
        {
          set_xsec(xs, xserr);
        }
        else
        {
          _xsec += xs;
          _xsecerr = HEPUtils::add_quad(xsec_err(), xserr);
        }
      }
    }

    /// Combine cross-sections and errors for the same process type, assuming uncorrelated errors.
    void Analysis::improve_xsec(double xs, double xserr)
    {
      if (xs > 0)
      {
        if (xsec() <= 0)
        {
          set_xsec(xs, xserr);
        }
        else
        {
          /// @todo Probably shouldn't be combined with equal weight?!?
          _xsec = _xsec/2.0 + xs/2.0;
          _xsecerr = HEPUtils::add_quad(xsec_err(), xserr) / 2.0;
        }
      }
    }

  }
}
