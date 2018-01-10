#include "gambit/ColliderBit/ColliderBit_macros.hpp"
#include "gambit/ColliderBit/analyses/HEPUtilsAnalysisContainer.hpp"
#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include <stdexcept>
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


    void HEPUtilsAnalysisContainer::clear()
    {
      /// @todo Storing smart ptrs to Analysis would make this way easier
      for (Analysis* a : analyses) {
        delete a;
        a = nullptr;
      }
      //if (analyses.size() != 0)
      //{
      //  for (auto it = analyses.begin(); it != analyses.end(); ++it)
      //  {
      //    delete *it;
      //    *it = nullptr;
      //  }
      //  analyses.clear();
      //}
      analyses.clear();
      ready = false;
    }


    void HEPUtilsAnalysisContainer::init(const std::vector<std::string>& analysisNames)
    {
      // assert(!analysisNames.empty()); //< @todo Really necessary? It's immediately emptied!!
      clear();
      for (const string& a : analysisNames)
        analyses.push_back(mkAnalysis(a));
      ready = true;
    }


    void HEPUtilsAnalysisContainer::reset()
    {
      ready = false;
      for (Analysis* a : analyses) a->reset();
      ready = true;
    }


    void HEPUtilsAnalysisContainer::analyze(const HEPUtils::Event& event) const
    {
      // assert(!analyses.empty()); //< @todo Really necessary?
      assert(ready); //< @todo Really necessary?
      for (Analysis* a : analyses) a->analyze(event);
    }


    void HEPUtilsAnalysisContainer::add_xsec(double xs, double xserr)
    {
      // assert(!analyses.empty()); //< @todo Really necessary?
      assert(ready); //< @todo Really necessary?
      for (Analysis* a : analyses) a->add_xsec(xs, xserr);
    }


    void HEPUtilsAnalysisContainer::add_xsec(const HEPUtilsAnalysisContainer* other)
    {
      assert(other->analyses.size() != 0); //< @todo Really necessary?
      assert(ready); //< @todo Really necessary?
      const Analysis* otherana = other->analyses.front();
      add_xsec(otherana->xsec(), otherana->xsec_err());
    }


    void HEPUtilsAnalysisContainer::improve_xsec(double xs, double xserr)
    {
      assert(!analyses.empty()); //< @todo Really necessary?
      assert(ready); //< @todo Really necessary?
      for (Analysis* a : analyses) a->improve_xsec(xs, xserr);
    }


    void HEPUtilsAnalysisContainer::improve_xsec(const HEPUtilsAnalysisContainer* other)
    {
      assert(other->analyses.size() != 0); //< @todo Really necessary?
      assert(ready); //< @todo Really necessary?
      const Analysis* otherana = other->analyses.front();
      improve_xsec(otherana->xsec(), otherana->xsec_err());
    }


    void HEPUtilsAnalysisContainer::add(const HEPUtilsAnalysisContainer* other)
    {
      assert(other->analyses.size() != 0); //< @todo Really necessary?
      assert(analyses.size() == other->analyses.size());
      assert(ready); //< @todo Really necessary?
      for (size_t i = 0; i < analyses.size(); ++i)
        analyses[i]->add(other->analyses[i]);
      // auto myIter = analyses.begin();
      // auto otherIter = other->analyses.begin();
      // while (myIter != analyses.end())
      // {
      //   (*myIter++)->add(*otherIter++);
      // }
    }


    void HEPUtilsAnalysisContainer::scale(double factor)
    {
      assert(!analyses.empty()); //< @todo Really necessary?
      assert(ready); //< @todo Really necessary?
      for (Analysis* a : analyses)
        a->scale(factor);
    }

  }
}
