#pragma once

#include <string>
#include <stdexcept>
#include <vector>
#include <map>

// Forward declarations, to avoid header-chaining into CB_types.hpp
namespace HEPUtils { class Event; }
namespace Gambit
{
  typedef std::string str;
  namespace ColliderBit
  {
    template <typename EventT>
    class BaseAnalysis;
    using HEPUtilsAnalysis = BaseAnalysis<HEPUtils::Event>;
  }
}


namespace Gambit
{

  namespace ColliderBit
  {

    /// Create a new analysis based on a name string
    /// @note The caller is responsible for deleting the returned analysis object.
    /// @todo Move to a separate file
    HEPUtilsAnalysis* mkAnalysis(const str& name);

    /// Return the detector to be used for a given analysis name, checking that the analysis exists.
    str getDetector(const str& name);


    /// A class for managing collections of HEPUtilsAnalysis instances.
    class HEPUtilsAnalysisContainer
    {

      private:

        /// A map of maps of pointer-to-HEPUtilsAnalysis.
        /// First key is the collider name, second key is the analysis name.
        std::map<str,std::map<str,HEPUtilsAnalysis*> > analyses_map;

        /// String identifying the currently active collider
        str current_collider;

        /// Has this class instance been initialized?
        bool ready; //< @todo Currently not used for anything. Do we need it?

        /// Has this instance been registered in the instances_map?
        bool is_registered;

        /// OpenMP info
        int n_threads;

        /// Key for the instances_map
        str base_key;

        /// A vector with pointers to all instances of this class. The key is the OMP thread number.
        /// (There should only be one instance of this class per OMP thread.)
        static std::map<str,std::map<int,HEPUtilsAnalysisContainer*> > instances_map;


      public:

        /// Constructor
        HEPUtilsAnalysisContainer();

        /// Destructor
        ~HEPUtilsAnalysisContainer();

        /// Add container to instances map
        void register_thread(str);

        /// Delete and clear the analyses contained within this instance.
        void clear();

        /// Set name of current collider
        void set_current_collider(str);

        /// Get name of current collider
        str get_current_collider() const;

        /// Does this instance contain analyses for the given collider
        bool has_analyses(str) const;
        /// Does this instance contain analyses for the current collider
        bool has_analyses() const;

        /// Initialize analyses (by names) for a specified collider
        void init(const std::vector<str>&, str);
        /// Initialize analyses (by names) for the current collider
        void init(const std::vector<str>&);

        /// Reset specific analysis
        void reset(str, str);
        /// Reset all analyses for given collider
        void reset(str);
        /// Reset all analyses for the current collider
        void reset();
        /// Reset all analyses for all colliders
        void reset_all();

        /// Get pointer to specific analysis
        const HEPUtilsAnalysis* get_analysis_pointer(str, str) const;
        /// Get analyses map for a specific collider
        const std::map<str,HEPUtilsAnalysis*>& get_collider_analyses_map(str) const;
        /// Get analyses map for the current collider
        const std::map<str,HEPUtilsAnalysis*>& get_current_analyses_map() const;
        /// Get the full analyses map
        const std::map<str,std::map<str,HEPUtilsAnalysis*> >& get_full_analyses_map() const;

        /// Pass event through specific analysis
        void analyze(const HEPUtils::Event&, str, str) const;
        /// Pass event through all analysis for a specific collider
        void analyze(const HEPUtils::Event&, str) const;
        /// Pass event through all analysis for the current collider
        void analyze(const HEPUtils::Event&) const;

        /// Add cross-sections and errors for two different processes,
        /// for a specific analysis
        void add_xsec(double, double, str, str);
        /// Add cross-sections and errors for two different processes,
        /// for all analyses for a given collider
        void add_xsec(double, double, str);
        /// Add cross-sections and errors for two different processes,
        /// for all analyses for the current collider
        void add_xsec(double, double);

        /// Weighted combination of xsecs and errors for the same process,
        /// for a specific analysis
        void improve_xsec(double, double, str, str);
        /// Weighted combination of xsecs and errors for the same process,
        /// for all analyses for a given collider
        void improve_xsec(double, double, str);
        /// Weighted combination of xsecs and errors for the same process,
        /// for all analyses for the current collider
        void improve_xsec(double, double);

        /// Collect signal predictions from other threads and add to this one,
        /// for specific analysis
        void collect_and_add_signal(str, str);
        /// Collect signal predictions from other threads and add to this one,
        /// for all analyses for given collider
        void collect_and_add_signal(str);
        /// Collect signal predictions from other threads and add to this one,
        /// for all analyses for the current collider
        void collect_and_add_signal();

        /// Collect xsec predictions from other threads and do a weighted combination,
        /// for specific analysis
        void collect_and_improve_xsec(str, str);
        /// Collect xsec predictions from other threads and do a weighted combination,
        /// for all analyses for given collider
        void collect_and_improve_xsec(str);
        /// Collect xsec predictions from other threads and do a weighted combination,
        /// for all analyses for the current collider
        void collect_and_improve_xsec();

        /// Scale results for specific analysis
        void scale(str, str, double factor=-1);
        /// Scale results for all analyses for given collider
        void scale(str, double factor=-1);
        /// Scale results for all analyses for the current collider
        void scale(double factor=-1);
        /// Scale results for all analyses across all colliders
        void scale_all(double factor=-1);

    };

  }
}
