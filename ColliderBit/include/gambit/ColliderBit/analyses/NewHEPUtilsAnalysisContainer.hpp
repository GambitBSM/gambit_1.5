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


    /// A class for managing collections of HEPUtilsAnalysis instances.
    class NewHEPUtilsAnalysisContainer
    {

      private:

        /// A map of maps of pointer-to-HEPUtilsAnalysis. 
        /// First key is the collider name, second key is the analysis name.
        std::map<string,std::map<string,HEPUtilsAnalysis*> > analyses_map;

        /// String identifying the currently active collider
        string current_collider;

        // /// A map indicating if a group of analyses has been properly initialized.
        // /// The key is the collider name.
        // std::map<string,bool> is_ready;

        /// Has this class instance been initialized?
        bool ready; //< @todo Currently not used for anything. Do we need it?

        /// OpenMP info
        int n_threads;
        int my_thread;

        /// A vector with pointers to all instances of this class. The key is the OMP thread number.
        /// (There should only be one instance of this class per OMP thread.)
        static std::map<int,NewHEPUtilsAnalysisContainer*> instances_map;


      public:

        /// Constructor
        NewHEPUtilsAnalysisContainer();

        /// Destructor
        ~NewHEPUtilsAnalysisContainer();

        /// Delete and clear the analyses contained within this instance.
        void clear();

        /// Set name of current collider
        void set_current_collider(string);

        /// Get name of current collider
        string get_current_collider() const;

        /// Initialize analyses (by names) for a specified collider
        void init(const std::vector<std::string>&, string);
        /// Initialize analyses (by names) for the current collider
        void init(const std::vector<std::string>&);

        /// Reset specific analysis
        void reset(string, string);
        /// Reset all analyses for given collider
        void reset(string);
        /// Reset all analyses for the current collider
        void reset();
        /// Reset all analyses for all colliders
        void reset_all();

        /// Get pointer to specific analysis
        HEPUtilsAnalysis* get_analysis_pointer(string, string) const;
        /// Get analyses map for a specific collider
        std::map<string,HEPUtilsAnalysis*>& get_collider_analyses_map(string) const;
        /// Get analyses map for the current collider
        std::map<string,HEPUtilsAnalysis*>& get_current_analyses_map() const;
        /// Get the full analyses map
        std::map<string,std::map<string,HEPUtilsAnalysis*> >& get_full_analyses_map() const;

        /// Pass event through specific analysis
        void analyze(const HEPUtils::Event&, string, string);
        /// Pass event through all analysis for a specific collider
        void analyze(const HEPUtils::Event&, string);
        /// Pass event through all analysis for the current collider
        void analyze(const HEPUtils::Event&);

        /// Add cross-sections and errors for two different processes, 
        /// for a specific analysis
        void add_xsec(double, double, string collider_name);
        /// Add cross-sections and errors for two different processes, 
        /// for all analyses for a given collider
        void add_xsec(double, double, string);
        /// Add cross-sections and errors for two different processes,
        /// for all analyses for the current collider
        void add_xsec(double, double);

        /// Weighted combination of xsecs and errors for the same process,
        /// for a specific analysis
        void improve_xsec(double, double, string collider_name);
        /// Weighted combination of xsecs and errors for the same process,
        /// for all analyses for a given collider
        void improve_xsec(double, double, string);
        /// Weighted combination of xsecs and errors for the same process,
        /// for all analyses for the current collider
        void improve_xsec(double, double);

        /// Collect signal predictions from other threads and add to this one,
        /// for specific analysis
        void collect_and_add_signal(string, string);
        /// Collect signal predictions from other threads and add to this one,
        /// for all analyses for given collider
        void collect_and_add_signal(string);
        /// Collect signal predictions from other threads and add to this one,
        /// for all analyses for the current collider
        void collect_and_add_signal();

        /// Collect xsec predictions from other threads and do a weighted combination,
        /// for specific analysis
        void collect_and_improve_xsec(string, string);
        /// Collect xsec predictions from other threads and do a weighted combination,
        /// for all analyses for given collider
        void collect_and_improve_xsec(string);
        /// Collect xsec predictions from other threads and do a weighted combination,
        /// for all analyses for the current collider
        void collect_and_improve_xsec();

        /// Scale results for specific analysis
        void scale(double factor=-1, string, string);
        /// Scale results for all analyses for given collider
        void scale(double factor=-1, string);
        /// Scale results for all analyses for the current collider
        void scale(double factor=-1);
        /// Scale results for all analyses across all colliders
        void scale_all(double factor=-1);

        // /// Add the results of all analyses from this instance to the given one.
        // void add(const HEPUtilsAnalysisContainer& e) { add(&e); }
        // /// Add the results of all analyses from this instance to the given one.
        // void add(const HEPUtilsAnalysisContainer*);

        // /// Collect results from other instances and improve xsec
        // void collect_results();

    };

  }
}
