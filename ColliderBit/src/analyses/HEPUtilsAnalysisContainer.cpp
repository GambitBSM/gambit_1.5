#include "gambit/ColliderBit/ColliderBit_macros.hpp"
#include "gambit/ColliderBit/analyses/HEPUtilsAnalysisContainer.hpp"
#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include <stdexcept>
using namespace std;

namespace Gambit
{
  namespace ColliderBit
  {


    /// Forward declarations using #DECLARE_ANALYSIS_FACTORY(ANAME)
    /// @{
    DECLARE_ANALYSIS_FACTORY(Minimum);
    DECLARE_ANALYSIS_FACTORY(Covariance);
    DECLARE_ANALYSIS_FACTORY(Perf);
    DECLARE_ANALYSIS_FACTORY(ATLAS_13TeV_2LEPStop_36invfb);
    DECLARE_ANALYSIS_FACTORY(ATLAS_13TeV_0LEP_13invfb);
    DECLARE_ANALYSIS_FACTORY(ATLAS_13TeV_MultiLEP_36invfb);
    DECLARE_ANALYSIS_FACTORY(ATLAS_13TeV_0LEPStop_20invfb);
    DECLARE_ANALYSIS_FACTORY(ATLAS_8TeV_0LEP_20invfb);
    DECLARE_ANALYSIS_FACTORY(ATLAS_8TeV_0LEPStop_20invfb);
    DECLARE_ANALYSIS_FACTORY(ATLAS_8TeV_1LEPStop_20invfb);
    DECLARE_ANALYSIS_FACTORY(ATLAS_8TeV_2bStop_20invfb);
    DECLARE_ANALYSIS_FACTORY(ATLAS_8TeV_2LEPEW_20invfb);
    DECLARE_ANALYSIS_FACTORY(ATLAS_8TeV_2LEPStop_20invfb);
    DECLARE_ANALYSIS_FACTORY(ATLAS_8TeV_3LEPEW_20invfb);
    DECLARE_ANALYSIS_FACTORY(ATLAS_8TeV_1LEPbb_20invfb);
    DECLARE_ANALYSIS_FACTORY(CMS_13TeV_MONOJET_36invfb);
    DECLARE_ANALYSIS_FACTORY(CMS_13TeV_0LEP_13invfb);
    DECLARE_ANALYSIS_FACTORY(CMS_13TeV_1LEPbb_36invfb);
    DECLARE_ANALYSIS_FACTORY(CMS_13TeV_MultiLEP_36invfb);
    DECLARE_ANALYSIS_FACTORY(CMS_13TeV_2OSLEP_36invfb);
    DECLARE_ANALYSIS_FACTORY(CMS_13TeV_2LEPsoft_36invfb);
    DECLARE_ANALYSIS_FACTORY(CMS_13TeV_1LEPStop_36invfb);
    DECLARE_ANALYSIS_FACTORY(CMS_13TeV_2LEPStop_36invfb);
    DECLARE_ANALYSIS_FACTORY(CMS_8TeV_1LEPDMTOP_20invfb);
    DECLARE_ANALYSIS_FACTORY(CMS_8TeV_2LEPDMTOP_20invfb);
    DECLARE_ANALYSIS_FACTORY(CMS_8TeV_3LEPEW_20invfb);
    DECLARE_ANALYSIS_FACTORY(CMS_8TeV_MONOJET_20invfb);
    /// @}


    // Factory definition
    HEPUtilsAnalysis* mkAnalysis(const std::string& name)
    {
      IF_X_RTN_CREATE_ANA_X(Minimum);
      IF_X_RTN_CREATE_ANA_X(Covariance);
      IF_X_RTN_CREATE_ANA_X(Perf);
      IF_X_RTN_CREATE_ANA_X(ATLAS_13TeV_2LEPStop_36invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_13TeV_0LEP_13invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_13TeV_MultiLEP_36invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_13TeV_0LEPStop_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_8TeV_0LEP_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_8TeV_0LEPStop_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_8TeV_1LEPStop_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_8TeV_2bStop_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_8TeV_2LEPEW_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_8TeV_2LEPStop_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_8TeV_3LEPEW_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_8TeV_1LEPbb_20invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_13TeV_MONOJET_36invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_13TeV_0LEP_13invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_13TeV_1LEPbb_36invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_13TeV_MultiLEP_36invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_13TeV_2OSLEP_36invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_13TeV_2LEPsoft_36invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_13TeV_1LEPStop_36invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_13TeV_2LEPStop_36invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_8TeV_1LEPDMTOP_20invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_8TeV_2LEPDMTOP_20invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_8TeV_3LEPEW_20invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_8TeV_MONOJET_20invfb);
      throw std::runtime_error(name + " isn't a known collider analysis!");
      return nullptr;
    }


    void HEPUtilsAnalysisContainer::clear()
    {
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
      assert(!analysisNames.empty());
      clear();
      for (const string& a : analysisNames)
        analyses.push_back(mkAnalysis(a));
      ready = true;
    }


    void HEPUtilsAnalysisContainer::analyze(const HEPUtils::Event& event) const
    {
      assert(!analyses.empty());
      assert(ready);
      for (Analysis* a : analyses) a->do_analysis(event);
    }


    void HEPUtilsAnalysisContainer::add_xsec(double xs, double xserr)
    {
      assert(!analyses.empty());
      assert(ready);
      for (Analysis* a : analyses) a->add_xsec(xs, xserr);
    }


    void HEPUtilsAnalysisContainer::add_xsec(const HEPUtilsAnalysisContainer* other)
    {
      assert(other->analyses.size() != 0);
      assert(ready);
      const Analysis* otherana = other->analyses.front();
      add_xsec(otherana->xsec(), otherana->xsec_err());
    }


    void HEPUtilsAnalysisContainer::improve_xsec(double xs, double xserr)
    {
      assert(!analyses.empty());
      assert(ready);
      for (Analysis* a : analyses) a->improve_xsec(xs, xserr);
    }


    void HEPUtilsAnalysisContainer::improve_xsec(const HEPUtilsAnalysisContainer* other)
    {
      assert(other->analyses.size() != 0);
      assert(ready);
      const Analysis* otherana = other->analyses.front();
      improve_xsec(otherana->xsec(), otherana->xsec_err());
    }


    void HEPUtilsAnalysisContainer::add(const HEPUtilsAnalysisContainer* other)
    {
      assert(other->analyses.size() != 0);
      assert(analyses.size() == other->analyses.size());
      assert(ready);
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
      assert(!analyses.empty());
      assert(ready);
      for (Analysis* a : analyses) a->scale(factor);
    }

  }
}
