#pragma once
//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  The SignalRegionData and AnalysisData structs.

#include "gambit/ColliderBit/ColliderBit_macros.hpp"
#include "gambit/ColliderBit/Utils.hpp"

#include "Eigen/Core"

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <limits>
#include <memory>
#include <iomanip>
#include <algorithm>

namespace Gambit {
  namespace ColliderBit {


    /// A simple container for the result of one signal region from one analysis.
    struct SignalRegionData
    {

      /// Constructor with {n,nsys} pair args
      SignalRegionData(const std::string& name, const std::string& sr,
                       double nobs, const std::pair<double,double>& nsig, const std::pair<double,double>& nbkg,
                       double nsigatlumi=-1)
       : SignalRegionData(name, sr, nobs, nsig.first, nbkg.first, nsig.second, nbkg.second, nsigatlumi)
      {}

      /// Constructor with separate n & nsys args
      SignalRegionData(const std::string& name, const std::string& sr,
                       double nobs, double nsig, double nbkg,
                       double syssig, double sysbkg, double nsigatlumi=-1)
       : analysis_name(name), sr_label(sr),
         n_observed(nobs), n_signal(nsig), n_signal_at_lumi(nsigatlumi), n_background(nbkg),
         signal_sys(syssig), background_sys(sysbkg)
      {}

      /// Default constructor
      SignalRegionData() {}

      /// @name Analysis and signal region specification
      //@{
      std::string analysis_name; ///< The name of the analysis common to all signal regions
      std::string sr_label; ///< A label for the particular signal region of the analysis
      //@}

      /// @name Signal region data
      //@{
      double n_observed = 0; ///< The number of events passing selection for this signal region as reported by the experiment
      double n_signal = 0; ///< The number of simulated model events passing selection for this signal region
      double n_signal_at_lumi = -1; ///< n_signal, scaled to the experimental luminosity
      double n_background = 0; ///< The number of standard model events expected to pass the selection for this signal region, as reported by the experiment.
      double signal_sys = 0; ///< The absolute systematic error of n_signal
      double background_sys = 0; ///< The absolute systematic error of n_background
      //@}

    };


    /// A container for the result of an analysis, potentially with many signal regions and correlations
    ///
    /// @todo Access by name?
    /// @todo Guarantee ordering?
    /// @todo How to combine covariance matrices -- require?
    struct AnalysisData
    {

      /// Default constructor
      AnalysisData() { clear(); }

      /// @brief Constructor from a list of SignalRegionData and an optional correlation (or covariance?) matrix
      ///
      /// If corrs is a null matrix (the default), this AnalysisData is to be interpreted as having no correlation
      /// information, and hence the likelihood calculation should use the single best-expected-limit SR.
      AnalysisData(const std::vector<SignalRegionData>& srds, const Eigen::MatrixXd& cov=Eigen::MatrixXd())
        : srdata(srds), srcov(cov)
      {
        _checkConsistency();
      }

      /// Clear the list of SignalRegionData, and nullify the covariance matrix
      /// @todo It'd be good to *not* have to re-enter most of the SRData and the covariance on every point: they don't change
      void clear()
      {
        srdata.clear();
        srcov = Eigen::MatrixXd();
      }

      /// Number of analyses
      size_t size() const
      {
        _checkConsistency();
        return srdata.size();
      }

      /// Is this container empty of signal regions?
      bool empty() const { return size() == 0; }

      /// Is there non-null correlation data?
      bool hasCorrs() const
      {
        _checkConsistency();
        return srcov.rows() == 0;
      }

      /// @brief Add a SignalRegionData
      /// @todo Allow naming the SRs?
      void add(const SignalRegionData& srd)
      {
        std::string key = srd.analysis_name + srd.sr_label;
        auto loc = srdata_identifiers.find(key);
        if (loc == srdata_identifiers.end())
        {
          // If the signal region doesn't exist in this object yet, add it
          srdata.push_back(srd);
          srdata_identifiers[key] = srdata.size();
        }
        else
        {
          // If it does, just update the signal count in the existing SignalRegionData object
          srdata[loc->second].n_signal = srd.n_signal;
        }
      }

      /// Access the i'th signal region's data
      SignalRegionData& operator[] (size_t i) { return srdata[i]; }
      /// Access the i'th signal region's data (const)
      const SignalRegionData& operator[] (size_t i) const { return srdata[i]; }

      /// Iterators (sugar for direct access to this->srdata)
      std::vector<SignalRegionData>::iterator begin() { return srdata.begin(); }
      std::vector<SignalRegionData>::const_iterator begin() const { return srdata.begin(); }
      std::vector<SignalRegionData>::iterator end() { return srdata.end(); }
      std::vector<SignalRegionData>::const_iterator end() const { return srdata.end(); }

      /// List of signal regions' data summaries
      std::vector<SignalRegionData> srdata;

      /// Map of names and indices of all entries in srdata, for easy lookup
      std::map<std::string, int> srdata_identifiers;

      /// Optional covariance matrix between SRs (0x0 null matrix = no correlation info)
      Eigen::MatrixXd srcov;

      /// Check that the size of the SRData list and the covariance matrix are consistent
      void _checkConsistency() const
      {
        assert(srcov.rows() == 0 || srcov.rows() == (int) srdata.size());
      }

    };


  }
}
