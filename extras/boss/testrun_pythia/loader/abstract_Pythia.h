#ifndef __ABSTRACT_PYTHIA_H__
#define __ABSTRACT_PYTHIA_H__

#include <iostream>  // FOR DEBUG: Allow virtual member functions to print a warning if executed.

#include "forward_decls_abstract_classes.h"
#include <string>
#include <istream>
#include "abstract_Event.h"
#include "abstract_Vec4.h"
#include <ostream>
#include "abstract_Info.h"
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wreturn-type"
namespace Pythia8
{
    class Abstract_Pythia
    {
        public:

            virtual Abstract_Event& process_ref_gambit() {std::cout << "Called virtual function" << std::endl;};

            virtual Abstract_Event& event_ref_gambit() {std::cout << "Called virtual function" << std::endl;};

            virtual Abstract_Info& info_ref_gambit() {std::cout << "Called virtual function" << std::endl;};
            // IGNORED: Field  -- Name: settings  -- XML id: _23802
            // IGNORED: Field  -- Name: particleData  -- XML id: _23803
            // IGNORED: Field  -- Name: rndm  -- XML id: _23804
            // IGNORED: Field  -- Name: couplings  -- XML id: _23805
            // IGNORED: Field  -- Name: couplingsPtr  -- XML id: _23806
            // IGNORED: Field  -- Name: slhaInterface  -- XML id: _23807
            // IGNORED: Field  -- Name: partonSystems  -- XML id: _23808
            // IGNORED: Field  -- Name: merging  -- XML id: _23809
            // IGNORED: Field  -- Name: mergingHooksPtr  -- XML id: _23810
        private:
            // IGNORED: Variable  -- Name: VERSIONNUMBERCODE  -- XML id: _23811
            // IGNORED: Variable  -- Name: NTRY  -- XML id: _23812
            // IGNORED: Variable  -- Name: SUBRUNDEFAULT  -- XML id: _23813
            // IGNORED: Field  -- Name: xmlPath  -- XML id: _23814
            // IGNORED: Field  -- Name: doProcessLevel  -- XML id: _23815
            // IGNORED: Field  -- Name: doPartonLevel  -- XML id: _23816
            // IGNORED: Field  -- Name: doHadronLevel  -- XML id: _23817
            // IGNORED: Field  -- Name: doDiffraction  -- XML id: _23818
            // IGNORED: Field  -- Name: doResDec  -- XML id: _23819
            // IGNORED: Field  -- Name: doFSRinRes  -- XML id: _23820
            // IGNORED: Field  -- Name: decayRHadrons  -- XML id: _23821
            // IGNORED: Field  -- Name: abortIfVeto  -- XML id: _23822
            // IGNORED: Field  -- Name: checkEvent  -- XML id: _23823
            // IGNORED: Field  -- Name: checkHistory  -- XML id: _23824
            // IGNORED: Field  -- Name: nErrList  -- XML id: _23825
            // IGNORED: Field  -- Name: epTolErr  -- XML id: _23826
            // IGNORED: Field  -- Name: epTolWarn  -- XML id: _23827
            // IGNORED: Field  -- Name: isConstructed  -- XML id: _23828
            // IGNORED: Field  -- Name: isInit  -- XML id: _23829
            // IGNORED: Field  -- Name: isUnresolvedA  -- XML id: _23830
            // IGNORED: Field  -- Name: isUnresolvedB  -- XML id: _23831
            // IGNORED: Field  -- Name: showSaV  -- XML id: _23832
            // IGNORED: Field  -- Name: showMaD  -- XML id: _23833
            // IGNORED: Field  -- Name: idA  -- XML id: _23834
            // IGNORED: Field  -- Name: idB  -- XML id: _23835
            // IGNORED: Field  -- Name: frameType  -- XML id: _23836
            // IGNORED: Field  -- Name: boostType  -- XML id: _23837
            // IGNORED: Field  -- Name: nCount  -- XML id: _23838
            // IGNORED: Field  -- Name: nShowLHA  -- XML id: _23839
            // IGNORED: Field  -- Name: nShowInfo  -- XML id: _23840
            // IGNORED: Field  -- Name: nShowProc  -- XML id: _23841
            // IGNORED: Field  -- Name: nShowEvt  -- XML id: _23842
            // IGNORED: Field  -- Name: mA  -- XML id: _23843
            // IGNORED: Field  -- Name: mB  -- XML id: _23844
            // IGNORED: Field  -- Name: pxA  -- XML id: _23845
            // IGNORED: Field  -- Name: pxB  -- XML id: _23846
            // IGNORED: Field  -- Name: pyA  -- XML id: _23847
            // IGNORED: Field  -- Name: pyB  -- XML id: _23848
            // IGNORED: Field  -- Name: pzA  -- XML id: _23849
            // IGNORED: Field  -- Name: pzB  -- XML id: _23850
            // IGNORED: Field  -- Name: eA  -- XML id: _23851
            // IGNORED: Field  -- Name: eB  -- XML id: _23852
            // IGNORED: Field  -- Name: pzAcm  -- XML id: _23853
            // IGNORED: Field  -- Name: pzBcm  -- XML id: _23854
            // IGNORED: Field  -- Name: eCM  -- XML id: _23855
            // IGNORED: Field  -- Name: betaZ  -- XML id: _23856
            // IGNORED: Field  -- Name: gammaZ  -- XML id: _23857
            // IGNORED: Field  -- Name: pAinit  -- XML id: _23858
            // IGNORED: Field  -- Name: pBinit  -- XML id: _23859
            // IGNORED: Field  -- Name: pAnow  -- XML id: _23860
            // IGNORED: Field  -- Name: pBnow  -- XML id: _23861
            // IGNORED: Field  -- Name: MfromCM  -- XML id: _23862
            // IGNORED: Field  -- Name: MtoCM  -- XML id: _23863
            // IGNORED: Field  -- Name: nErrEvent  -- XML id: _23864
            // IGNORED: Field  -- Name: iErrId  -- XML id: _23865
            // IGNORED: Field  -- Name: iErrCol  -- XML id: _23866
            // IGNORED: Field  -- Name: iErrEpm  -- XML id: _23867
            // IGNORED: Field  -- Name: iErrNan  -- XML id: _23868
            // IGNORED: Field  -- Name: iErrNanVtx  -- XML id: _23869
            // IGNORED: Field  -- Name: pdfAPtr  -- XML id: _23870
            // IGNORED: Field  -- Name: pdfBPtr  -- XML id: _23871
            // IGNORED: Field  -- Name: pdfHardAPtr  -- XML id: _23872
            // IGNORED: Field  -- Name: pdfHardBPtr  -- XML id: _23873
            // IGNORED: Field  -- Name: pdfPomAPtr  -- XML id: _23874
            // IGNORED: Field  -- Name: pdfPomBPtr  -- XML id: _23875
            // IGNORED: Field  -- Name: useNewPdfA  -- XML id: _23876
            // IGNORED: Field  -- Name: useNewPdfB  -- XML id: _23877
            // IGNORED: Field  -- Name: useNewPdfHard  -- XML id: _23878
            // IGNORED: Field  -- Name: useNewPdfPomA  -- XML id: _23879
            // IGNORED: Field  -- Name: useNewPdfPomB  -- XML id: _23880
            // IGNORED: Field  -- Name: beamA  -- XML id: _23881
            // IGNORED: Field  -- Name: beamB  -- XML id: _23882
            // IGNORED: Field  -- Name: beamPomA  -- XML id: _23883
            // IGNORED: Field  -- Name: beamPomB  -- XML id: _23884
            // IGNORED: Field  -- Name: doLHA  -- XML id: _23885
            // IGNORED: Field  -- Name: useNewLHA  -- XML id: _23886
            // IGNORED: Field  -- Name: lhaUpPtr  -- XML id: _23887
            // IGNORED: Field  -- Name: decayHandlePtr  -- XML id: _23888
            // IGNORED: Field  -- Name: handledParticles  -- XML id: _23889
            // IGNORED: Field  -- Name: userHooksPtr  -- XML id: _23890
            // IGNORED: Field  -- Name: hasUserHooks  -- XML id: _23891
            // IGNORED: Field  -- Name: doVetoProcess  -- XML id: _23892
            // IGNORED: Field  -- Name: doVetoPartons  -- XML id: _23893
            // IGNORED: Field  -- Name: retryPartonLevel  -- XML id: _23894
            // IGNORED: Field  -- Name: beamShapePtr  -- XML id: _23895
            // IGNORED: Field  -- Name: useNewBeamShape  -- XML id: _23896
            // IGNORED: Field  -- Name: doMomentumSpread  -- XML id: _23897
            // IGNORED: Field  -- Name: doVertexSpread  -- XML id: _23898
            // IGNORED: Field  -- Name: sigmaPtrs  -- XML id: _23899
            // IGNORED: Field  -- Name: resonancePtrs  -- XML id: _23900
            // IGNORED: Field  -- Name: timesDecPtr  -- XML id: _23901
            // IGNORED: Field  -- Name: timesPtr  -- XML id: _23902
            // IGNORED: Field  -- Name: spacePtr  -- XML id: _23903
            // IGNORED: Field  -- Name: useNewTimes  -- XML id: _23904
            // IGNORED: Field  -- Name: useNewSpace  -- XML id: _23905
            // IGNORED: Field  -- Name: processLevel  -- XML id: _23906
            // IGNORED: Field  -- Name: partonLevel  -- XML id: _23907
            // IGNORED: Field  -- Name: trialPartonLevel  -- XML id: _23908
            // IGNORED: Field  -- Name: hasMergingHooks  -- XML id: _23909
            // IGNORED: Field  -- Name: hasOwnMergingHooks  -- XML id: _23910
            // IGNORED: Field  -- Name: doMerging  -- XML id: _23911
            // IGNORED: Field  -- Name: hadronLevel  -- XML id: _23912
            // IGNORED: Field  -- Name: sigmaTot  -- XML id: _23913
            // IGNORED: Field  -- Name: rHadrons  -- XML id: _23914
        public:

            virtual bool readString(std::string line, bool warn) {std::cout << "Called virtual function" << std::endl;};

            virtual bool readFile(std::string fileName, bool warn, int subrun) {std::cout << "Called virtual function" << std::endl;};

            virtual bool readFile(std::string fileName, int subrun) {std::cout << "Called virtual function" << std::endl;};

            virtual bool readFile(std::istream& is, bool warn, int subrun) {std::cout << "Called virtual function" << std::endl;};

            virtual bool readFile(std::istream& is, int subrun) {std::cout << "Called virtual function" << std::endl;};

            virtual bool init() {std::cout << "Called virtual function" << std::endl;};

            virtual bool init(int idAin, int idBin, double eCMin) {std::cout << "Called virtual function" << std::endl;};

            virtual bool init(int idAin, int idBin, double eAin, double eBin) {std::cout << "Called virtual function" << std::endl;};

            virtual bool init(int idAin, int idBin, double pxAin, double pyAin, double pzAin, double pxBin, double pyBin, double pzBin) {std::cout << "Called virtual function" << std::endl;};

            virtual bool init(std::string LesHouchesEventFile, bool skipInit) {std::cout << "Called virtual function" << std::endl;};

            virtual bool next() {std::cout << "Called virtual function" << std::endl;};

            virtual int forceTimeShower(int iBeg, int iEnd, double pTmax, int nBranchMax) {std::cout << "Called virtual function" << std::endl;};

            virtual bool forceHadronLevel(bool findJunctions) {std::cout << "Called virtual function" << std::endl;};

            virtual bool moreDecays() {std::cout << "Called virtual function" << std::endl;};

            virtual bool forceRHadronDecays() {std::cout << "Called virtual function" << std::endl;};

            virtual void LHAeventList(std::ostream& os) {std::cout << "Called virtual function" << std::endl;};

            virtual bool LHAeventSkip(int nSkip) {std::cout << "Called virtual function" << std::endl;};

            virtual void stat() {std::cout << "Called virtual function" << std::endl;};

            virtual void statistics(bool all, bool reset) {std::cout << "Called virtual function" << std::endl;};

            virtual bool flag(std::string key) {std::cout << "Called virtual function" << std::endl;};

            virtual int mode(std::string key) {std::cout << "Called virtual function" << std::endl;};

            virtual double parm(std::string key) {std::cout << "Called virtual function" << std::endl;};

            virtual std::string word(std::string key) {std::cout << "Called virtual function" << std::endl;};
        private:

            virtual Pythia8::Abstract_Pythia* operator_assignment_gambit(const Pythia8::Abstract_Pythia& arg_1) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Pythia* operator=(const Pythia8::Abstract_Pythia& arg_1)
            {
                return operator_assignment_gambit(arg_1);
            }

            virtual void banner(std::ostream& os) {std::cout << "Called virtual function" << std::endl;};

            virtual int readSubrun(std::string line, bool warn, std::ostream& os) {std::cout << "Called virtual function" << std::endl;};

            virtual int readCommented(std::string line) {std::cout << "Called virtual function" << std::endl;};

            virtual void checkSettings() {std::cout << "Called virtual function" << std::endl;};

            virtual bool checkBeams() {std::cout << "Called virtual function" << std::endl;};

            virtual bool initKinematics() {std::cout << "Called virtual function" << std::endl;};

            virtual bool initPDFs() {std::cout << "Called virtual function" << std::endl;};

            virtual void nextKinematics() {std::cout << "Called virtual function" << std::endl;};

            virtual void boostAndVertex(bool toLab, bool setVertex) {std::cout << "Called virtual function" << std::endl;};

            virtual bool doRHadronDecays() {std::cout << "Called virtual function" << std::endl;};

            virtual bool check(std::ostream& os) {std::cout << "Called virtual function" << std::endl;};

            virtual bool initSLHA() {std::cout << "Called virtual function" << std::endl;};

        public:
            virtual Abstract_Pythia* pointerCopy_gambit() {std::cout << "Called virtual function" << std::endl;};
            virtual ~Abstract_Pythia() {};
    };
}
#pragma GCC diagnostic pop


#endif /* __ABSTRACT_PYTHIA_H__ */
