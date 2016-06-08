#include "gambit/ColliderBit/ColliderBit_macros.hpp"
#include "gambit/ColliderBit/analyses/HEPUtilsAnalysisContainer.hpp"
#include <stdexcept>
using namespace std;

namespace Gambit {
  namespace ColliderBit {

    /// @todo Move these to a separate file

    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(Minimum);
    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(ATLAS_0LEP_20invfb);
    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(ATLAS_0LEPStop_20invfb);
    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(ATLAS_1LEPStop_20invfb);
    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(ATLAS_2bStop_20invfb);
    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(ATLAS_2LEPEW_20invfb);
    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(ATLAS_2LEPStop_20invfb);
    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(ATLAS_3LEPEW_20invfb);
    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(CMS_1LEPDMTOP_20invfb);
    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(CMS_2LEPDMTOP_20invfb);
    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(CMS_3LEPEW_20invfb);
    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(CMS_MONOJET_20invfb);
    /// Forward declaration using #DECLARE_ANALYSIS_FACTORY(ANAME)
    DECLARE_ANALYSIS_FACTORY(Perf);

    // Factory definition
    /// @todo Move to a separate file
    HEPUtilsAnalysis* mkAnalysis(const string& name) {
      IF_X_RTN_CREATE_ANA_X(Minimum);
      IF_X_RTN_CREATE_ANA_X(ATLAS_0LEP_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_0LEPStop_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_1LEPStop_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_2bStop_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_2LEPEW_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_2LEPStop_20invfb);
      IF_X_RTN_CREATE_ANA_X(ATLAS_3LEPEW_20invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_1LEPDMTOP_20invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_2LEPDMTOP_20invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_3LEPEW_20invfb);
      IF_X_RTN_CREATE_ANA_X(CMS_MONOJET_20invfb);
      IF_X_RTN_CREATE_ANA_X(Perf);
      throw runtime_error(name + " isn't a known collider analysis, you fool of a Took!");
      return nullptr;
    }


    void HEPUtilsAnalysisContainer::clear() {
      for (Analysis* a : analyses) {
        delete a;
        a = nullptr;
      }
      analyses.clear();
      ready = false;
    }


    void HEPUtilsAnalysisContainer::init(const vector<string>& analysisNames) {
      assert(!analysisNames.empty());
      clear();
      for (const string& a : analysisNames)
        analyses.push_back(mkAnalysis(a));
      ready = true;
    }


    void HEPUtilsAnalysisContainer::analyze(const HEPUtils::Event& event) const {
      assert(!analyses.empty());
      assert(ready);
      for (Analysis* a : analyses) a->analyze(event);
    }


    void HEPUtilsAnalysisContainer::add_xsec(double xs, double xserr) {
      assert(!analyses.empty());
      assert(ready);
      for (Analysis* a : analyses) a->add_xsec(xs, xserr);
    }


    void HEPUtilsAnalysisContainer::add_xsec(const HEPUtilsAnalysisContainer* other) {
      assert(other->analyses.size() != 0);
      assert(ready);
      const Analysis* otherana = other->analyses.front();
      add_xsec(otherana->xsec(), otherana->xsec_err());
    }


    void HEPUtilsAnalysisContainer::improve_xsec(double xs, double xserr) {
      assert(!analyses.empty());
      assert(ready);
      for (Analysis* a : analyses) a->improve_xsec(xs, xserr);
    }


    void HEPUtilsAnalysisContainer::improve_xsec(const HEPUtilsAnalysisContainer* other) {
      assert(other->analyses.size() != 0);
      assert(ready);
      const Analysis* otherana = other->analyses.front();
      improve_xsec(otherana->xsec(), otherana->xsec_err());
    }


    void HEPUtilsAnalysisContainer::add(const HEPUtilsAnalysisContainer* other) {
      assert(other->analyses.size() != 0);
      assert(analyses.size() == other->analyses.size());
      assert(ready);
      for (size_t i = 0; i < analyses.size(); ++i)
        analyses[i]->add(other->analyses[i]);
    }


    void HEPUtilsAnalysisContainer::scale(double factor) {
      assert(!analyses.empty());
      assert(ready);
      for (Analysis* a : analyses) a->scale(factor);
    }

  }
}
