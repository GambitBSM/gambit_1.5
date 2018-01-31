#include "gambit/ColliderBit/ColliderBit_macros.hpp"
#include "gambit/ColliderBit/analyses/NewHEPUtilsAnalysisContainer.hpp"
#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include <stdexcept>
// _Anders
#include <omp.h>
using namespace std;

namespace Gambit
{
  namespace ColliderBit
  {

    // Add analysis names here and only here (trick to avoid duplication)
    #define MAP_ANALYSES(F)                          \
      F(Minimum)                                     \
      F(Covariance)                                  \
      F(Perf)                                        \
      F(ATLAS_13TeV_0LEP_13invfb)                    \
      F(ATLAS_13TeV_MultiLEP_36invfb)                \
      F(ATLAS_13TeV_0LEPStop_20invfb)                \
      F(ATLAS_13TeV_1LEPStop_36invfb)                \
      F(ATLAS_13TeV_RJ3L_lowmass_36invfb)            \
      F(ATLAS_8TeV_0LEP_20invfb)                     \
      F(ATLAS_8TeV_0LEPStop_20invfb)                 \
      F(ATLAS_8TeV_1LEPStop_20invfb)                 \
      F(ATLAS_8TeV_2bStop_20invfb)                   \
      F(ATLAS_8TeV_2LEPEW_20invfb)                   \
      F(ATLAS_8TeV_2LEPStop_20invfb)                 \
      F(ATLAS_8TeV_3LEPEW_20invfb)                   \
      F(ATLAS_8TeV_1LEPbb_20invfb)                   \
      F(CMS_13TeV_0LEP_13invfb)                      \
      F(CMS_13TeV_1LEPbb_36invfb)                    \
      F(CMS_13TeV_MultiLEP_36invfb)                  \
      F(CMS_13TeV_2OSLEP_36invfb)                    \
      F(CMS_13TeV_2LEPsoft_36invfb)                  \
      F(CMS_8TeV_1LEPDMTOP_20invfb)                  \
      F(CMS_8TeV_2LEPDMTOP_20invfb)                  \
      F(CMS_8TeV_3LEPEW_20invfb)                     \
      F(CMS_8TeV_MONOJET_20invfb)
    


    /// Forward declarations using #DECLARE_ANALYSIS_FACTORY(ANAME)

    MAP_ANALYSES(DECLARE_ANALYSIS_FACTORY)


    // Factory definition
    HEPUtilsAnalysis* mkAnalysis(const std::string& name)
    {

      MAP_ANALYSES(IF_X_RTN_CREATE_ANA_X)

      throw std::runtime_error(name + " isn't a known collider analysis!");
      return nullptr;
    }


    /// A map with pointers to all instances of this class. The key is the thread number.
    std::map<int,NewHEPUtilsAnalysisContainer*> NewHEPUtilsAnalysisContainer::instances_map;

    /// Constructor
    NewHEPUtilsAnalysisContainer::NewHEPUtilsAnalysisContainer() : 
      ready(false),
      n_threads(omp_get_max_threads()),
      my_thread(omp_get_thread_num()),
      current_collider("")
    { 
      #pragma omp critical
      {
        // Check that no other instances exist for the current OMP thread
        if (instances_map.count(my_thread) > 0)
        {
          utils_error().raise(LOCAL_INFO, "Only one instance of NewHEPUtilsAnalysisContainer allowed per OpenMP thread.");
        }

        // Add this instance to the instances map
        instances_map[my_thread] = this;
      }
    }


    /// Destructor
    NewHEPUtilsAnalysisContainer::~NewHEPUtilsAnalysisContainer() 
    { 
      clear(); 
    }


    /// Delete and clear the analyses contained within this instance.
    void NewHEPUtilsAnalysisContainer::clear()
    {
      /// @todo Storing smart ptrs to Analysis would make this way easier
      // Loop through double map and delete the analysis pointers
      for(auto& collider_map_pair : analyses_map)
      {
        for(auto& analysis_pointer_pair : collider_map_pair.second)
        {
          delete analysis_pointer_pair.second;
          analysis_pointer_pair.second = nullptr;
        }
      }

      // Clear the double map
      analyses_map.clear();
      ready = false;
    }


    /// Set name of the current collider
    void NewHEPUtilsAnalysisContainer::set_current_collider(string collider_name)
    {
      current_collider = collider_name;
    }


    /// Get the name of the current collider
    string NewHEPUtilsAnalysisContainer::get_current_collider() const
    {
      return current_collider;
    }


    /// Initialize analyses (by names) for a specified collider
    void NewHEPUtilsAnalysisContainer::init(const std::vector<std::string>& analysis_names, string collider_name);
    {
      // If a map of analyses already exist for this collider, clear it
      if (analyses_map.count(collider_name) > 0)
      {
        analyses_map[collider_name].clear();
      }

      // Create analysis pointers and add to the map 
      for (auto& aname : analysis_names)
      {
        analyses_map[collider_name][aname] = mkAnalysis(aname);
      }

      // This instance is now ready
      ready = true;
    }

    /// Initialize analyses (by names) for the current collider
    void NewHEPUtilsAnalysisContainer::init(const std::vector<std::string>& analysis_names);
    {
      init(analysis_names, current_collider);
    }


    /// Reset specific analysis
    void NewHEPUtilsAnalysisContainer::reset(string collider_name, string analysis_name)
    {
      analyses_map[collider_name][analysis_name]->reset();
    }

    /// Reset all analyses for given collider
    void NewHEPUtilsAnalysisContainer::reset(string collider_name)
    {
      for (auto& analysis_pointer_pair : analyses_map[collider_name])
      {
        analysis_pointer_pair.second->reset();
      }
    }

    /// Reset all analyses for the current collider
    void NewHEPUtilsAnalysisContainer::reset()
    {
      reset_analyses(current_collider);
    }

    /// Reset all analyses for all colliders
    void NewHEPUtilsAnalysisContainer::reset_all()
    {
      for(auto& collider_map_pair : analyses_map)
      {
        for(auto& analysis_pointer_pair : collider_map_pair.second)
        {
          analysis_pointer_pair.second->reset();
        }
      }
    }


    /// Get pointer to specific analysis
    HEPUtilsAnalysis* NewHEPUtilsAnalysisContainer::get_analysis_pointer(string collider_name, string analysis_name) const
    {
      return analyses_map[collider_name][analysis_name];
    }

    /// Get analyses map for a specific collider
    std::map<string,HEPUtilsAnalysis*>& get_collider_analyses_map(string collider_name) const;
    {
      return analysis_name[collider_name];
    }

    /// Get analyses map for the current collider
    std::map<string,HEPUtilsAnalysis*>& get_current_analyses_map() const;
    {
      return analysis_name[current_collider];
    }

    /// Get the full analyses map
    std::map<string,std::map<string,HEPUtilsAnalysis*> >& NewHEPUtilsAnalysisContainer::get_full_analyses_map() const
    {
      return analyses_map;
    }


    /// Pass event through specific analysis
    void NewHEPUtilsAnalysisContainer::analyze(const HEPUtils::Event& event, string collider_name, string analysis_name)
    { 
      analyses_map[collider_name][analysis_name]->do_analysis(event);
    }

    /// Pass event through all analysis for a specific collider
    void NewHEPUtilsAnalysisContainer::analyze(const HEPUtils::Event& event, string collider_name)
    {
      for (auto& analysis_pointer_pair : analyses_map[collider_name])
      {
        analysis_pointer_pair.second->do_analysis(event);
      }
    }

    /// Pass event through all analysis for the current collider
    void NewHEPUtilsAnalysisContainer::analyze(const HEPUtils::Event& event)
    { 
      analyze(event, current_collider);
    }


    /// Add cross-sections and errors for two different processes,
    /// for specific analysis
    void NewHEPUtilsAnalysisContainer::add_xsec(double xs, double xserr, string collider_name, string analysis_name)
    {
      analyses_map[collider_name][analysis_name]->add_xsec(xs, xserr);
    }

    /// Add cross-sections and errors for two different processes,
    /// for all analyses for a given collider
    void NewHEPUtilsAnalysisContainer::add_xsec(double xs, double xserr, string collider_name)
    {
      for (auto& analysis_pointer_pair : analyses_map[collider_name])
      {
        analysis_pointer_pair.second->add_xsec(xs, xserr);
      }
    }

    /// Add cross-sections and errors for two different processes,
    /// for all analyses for the current collider
    void NewHEPUtilsAnalysisContainer::add_xsec(double xs, double xserr)
    {
      add_xsec(xs, xserr, current_collider);
    }


    /// Weighted combination of cross-sections and errors for the same process,
    /// for a specific analysis
    void NewHEPUtilsAnalysisContainer::improve_xsec(double xs, double xserr, string collider_name, string analysis_name)
    {
      analyses_map[collider_name][analysis_name]->improve_xsec(xs, xserr);
    }

    /// Weighted combination of cross-sections and errors for the same process,
    /// for all analyses for a given collider
    void NewHEPUtilsAnalysisContainer::improve_xsec(double xs, double xserr, string collider_name)
    {
      for (auto& analysis_pointer_pair : analyses_map[collider_name])
      {
        analysis_pointer_pair.second->improve_xsec(xs, xserr);
      }
    }

    /// Weighted combination of cross-sections and errors for the same process,
    /// for all analyses for the current collider
    void NewHEPUtilsAnalysisContainer::improve_xsec(double xs, double xserr)
    {
      improve_xsec(xs, xserr, current_collider);
    }


    // 
    // @todo Add the 'collect_and_add_signal' functions
    // 

    /// Collect signal predictions from other threads and add to this one,
    /// for specific analysis
    void NewHEPUtilsAnalysisContainer::collect_and_add_signal(string collider_name, string analysis_name)
    {
      for (auto& thread_container_pair : instances_map)
      {
        if (thread == my_thread) continue;

        NewHEPUtilsAnalysisContainer* other_container = thread_container_pair.second;
        HEPUtilsAnalysis* other_analysis = other_container->get_analysis_pointer(collider_name, analysis_name);

        analyses_map[collider_name][analysis_name]->add(other_analysis);
    }

    /// Collect signal predictions from other threads and add to this one,
    /// for all analyses for given collider
    void NewHEPUtilsAnalysisContainer::collect_and_add_signal(string collider_name)
    {
      for (auto& analysis_pointer_pair : analyses_map[collider_name])
      {
        string analysis_name = analysis_pointer_pair.first;
        collect_and_add_signal(collider_name, analysis_name);
      }
    }

    /// Collect signal predictions from other threads and add to this one,
    /// for all analyses for the current collider
    void NewHEPUtilsAnalysisContainer::collect_and_add_signal()
    {
      collect_and_add_signal(current_collider);
    }


    /// Collect xsec predictions from other threads and do a weighted combination,
    /// for specific analysis
    void NewHEPUtilsAnalysisContainer::collect_and_improve_xsec(string collider_name, string analysis_name)
    {
      for (auto& thread_container_pair : instances_map)
      {
        if (thread == my_thread) continue;

        NewHEPUtilsAnalysisContainer* other_container = thread_container_pair.second;
        HEPUtilsAnalysis* other_analysis = other_container->get_analysis_pointer(collider_name, analysis_name);

        double other_xsec = other_analysis->xsec();
        double other_xsec_err = other_analysis->xsec_err();

        improve_xsec(other_xsec, other_xsec_err, collider_name, analysis_name);
    }

    /// Collect xsec predictions from other threads and do a weighted combination,
    /// for all analyses for given collider
    void NewHEPUtilsAnalysisContainer::collect_and_improve_xsec(string collider_name)
    {
      for (auto& analysis_pointer_pair : analyses_map[collider_name])
      {
        string analysis_name = analysis_pointer_pair.first;
        collect_and_improve_xsec(collider_name, analysis_name);
      }
    }

    /// Collect xsec predictions from other threads and do a weighted combination,
    /// for all analyses for the current collider
    void NewHEPUtilsAnalysisContainer::collect_and_improve_xsec(string collider_name)
    {
      collect_and_improve_xsec(current_collider);
    }


    /// Scale results for specific analysis
    void NewHEPUtilsAnalysisContainer::scale(double factor, string collider_name, string analysis_name)
    {
      analyses_map[collider_name][analysis_name]->scale(factor);
    }

    /// Scale results for all analyses for given collider
    void NewHEPUtilsAnalysisContainer::scale(double factor, string collider_name)
    {
      for (auto& analysis_pointer_pair : analyses_map[collider_name])
      {
        string analysis_name = analysis_pointer_pair.first;
        scale(factor, collider_name, analysis_name);
      }      
    }

    /// Scale results for all analyses for the current collider
    void NewHEPUtilsAnalysisContainer::scale(double factor)
    {
      scale(factor, current_collider);
    }

    /// Scale results for all analyses across all colliders
    void NewHEPUtilsAnalysisContainer::scale_all(double factor)
    {
      for (auto& collider_map_pair : analyses_map)
      {
        string collider_name = collider_map_pair.first;
        scale(factor, collider_name);
      }
    }

  }
}
