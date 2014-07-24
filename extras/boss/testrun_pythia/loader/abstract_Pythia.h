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
            // IGNORED: Field  -- Name: settings  -- XML id: _24819
            // IGNORED: Field  -- Name: particleData  -- XML id: _24820
            // IGNORED: Field  -- Name: rndm  -- XML id: _24821
            // IGNORED: Field  -- Name: couplings  -- XML id: _24822
            // IGNORED: Field  -- Name: couplingsPtr  -- XML id: _24823
            // IGNORED: Field  -- Name: slhaInterface  -- XML id: _24824
            // IGNORED: Field  -- Name: partonSystems  -- XML id: _24825
            // IGNORED: Field  -- Name: merging  -- XML id: _24826
            // IGNORED: Field  -- Name: mergingHooksPtr  -- XML id: _24827
        private:
            // IGNORED: Variable  -- Name: VERSIONNUMBERCODE  -- XML id: _24828
            // IGNORED: Variable  -- Name: NTRY  -- XML id: _24829
            // IGNORED: Variable  -- Name: SUBRUNDEFAULT  -- XML id: _24830
            // IGNORED: Field  -- Name: xmlPath  -- XML id: _24831
            // IGNORED: Field  -- Name: doProcessLevel  -- XML id: _24832
            // IGNORED: Field  -- Name: doPartonLevel  -- XML id: _24833
            // IGNORED: Field  -- Name: doHadronLevel  -- XML id: _24834
            // IGNORED: Field  -- Name: doDiffraction  -- XML id: _24835
            // IGNORED: Field  -- Name: doResDec  -- XML id: _24836
            // IGNORED: Field  -- Name: doFSRinRes  -- XML id: _24837
            // IGNORED: Field  -- Name: decayRHadrons  -- XML id: _24838
            // IGNORED: Field  -- Name: abortIfVeto  -- XML id: _24839
            // IGNORED: Field  -- Name: checkEvent  -- XML id: _24840
            // IGNORED: Field  -- Name: checkHistory  -- XML id: _24841
            // IGNORED: Field  -- Name: nErrList  -- XML id: _24842
            // IGNORED: Field  -- Name: epTolErr  -- XML id: _24843
            // IGNORED: Field  -- Name: epTolWarn  -- XML id: _24844
            // IGNORED: Field  -- Name: isConstructed  -- XML id: _24845
            // IGNORED: Field  -- Name: isInit  -- XML id: _24846
            // IGNORED: Field  -- Name: isUnresolvedA  -- XML id: _24847
            // IGNORED: Field  -- Name: isUnresolvedB  -- XML id: _24848
            // IGNORED: Field  -- Name: showSaV  -- XML id: _24849
            // IGNORED: Field  -- Name: showMaD  -- XML id: _24850
            // IGNORED: Field  -- Name: idA  -- XML id: _24851
            // IGNORED: Field  -- Name: idB  -- XML id: _24852
            // IGNORED: Field  -- Name: frameType  -- XML id: _24853
            // IGNORED: Field  -- Name: boostType  -- XML id: _24854
            // IGNORED: Field  -- Name: nCount  -- XML id: _24855
            // IGNORED: Field  -- Name: nShowLHA  -- XML id: _24856
            // IGNORED: Field  -- Name: nShowInfo  -- XML id: _24857
            // IGNORED: Field  -- Name: nShowProc  -- XML id: _24858
            // IGNORED: Field  -- Name: nShowEvt  -- XML id: _24859
            // IGNORED: Field  -- Name: mA  -- XML id: _24860
            // IGNORED: Field  -- Name: mB  -- XML id: _24861
            // IGNORED: Field  -- Name: pxA  -- XML id: _24862
            // IGNORED: Field  -- Name: pxB  -- XML id: _24863
            // IGNORED: Field  -- Name: pyA  -- XML id: _24864
            // IGNORED: Field  -- Name: pyB  -- XML id: _24865
            // IGNORED: Field  -- Name: pzA  -- XML id: _24866
            // IGNORED: Field  -- Name: pzB  -- XML id: _24867
            // IGNORED: Field  -- Name: eA  -- XML id: _24868
            // IGNORED: Field  -- Name: eB  -- XML id: _24869
            // IGNORED: Field  -- Name: pzAcm  -- XML id: _24870
            // IGNORED: Field  -- Name: pzBcm  -- XML id: _24871
            // IGNORED: Field  -- Name: eCM  -- XML id: _24872
            // IGNORED: Field  -- Name: betaZ  -- XML id: _24873
            // IGNORED: Field  -- Name: gammaZ  -- XML id: _24874
            // IGNORED: Field  -- Name: pAinit  -- XML id: _24875
            // IGNORED: Field  -- Name: pBinit  -- XML id: _24876
            // IGNORED: Field  -- Name: pAnow  -- XML id: _24877
            // IGNORED: Field  -- Name: pBnow  -- XML id: _24878
            // IGNORED: Field  -- Name: MfromCM  -- XML id: _24879
            // IGNORED: Field  -- Name: MtoCM  -- XML id: _24880
            // IGNORED: Field  -- Name: nErrEvent  -- XML id: _24881
            // IGNORED: Field  -- Name: iErrId  -- XML id: _24882
            // IGNORED: Field  -- Name: iErrCol  -- XML id: _24883
            // IGNORED: Field  -- Name: iErrEpm  -- XML id: _24884
            // IGNORED: Field  -- Name: iErrNan  -- XML id: _24885
            // IGNORED: Field  -- Name: iErrNanVtx  -- XML id: _24886
            // IGNORED: Field  -- Name: pdfAPtr  -- XML id: _24887
            // IGNORED: Field  -- Name: pdfBPtr  -- XML id: _24888
            // IGNORED: Field  -- Name: pdfHardAPtr  -- XML id: _24889
            // IGNORED: Field  -- Name: pdfHardBPtr  -- XML id: _24890
            // IGNORED: Field  -- Name: pdfPomAPtr  -- XML id: _24891
            // IGNORED: Field  -- Name: pdfPomBPtr  -- XML id: _24892
            // IGNORED: Field  -- Name: useNewPdfA  -- XML id: _24893
            // IGNORED: Field  -- Name: useNewPdfB  -- XML id: _24894
            // IGNORED: Field  -- Name: useNewPdfHard  -- XML id: _24895
            // IGNORED: Field  -- Name: useNewPdfPomA  -- XML id: _24896
            // IGNORED: Field  -- Name: useNewPdfPomB  -- XML id: _24897
            // IGNORED: Field  -- Name: beamA  -- XML id: _24898
            // IGNORED: Field  -- Name: beamB  -- XML id: _24899
            // IGNORED: Field  -- Name: beamPomA  -- XML id: _24900
            // IGNORED: Field  -- Name: beamPomB  -- XML id: _24901
            // IGNORED: Field  -- Name: doLHA  -- XML id: _24902
            // IGNORED: Field  -- Name: useNewLHA  -- XML id: _24903
            // IGNORED: Field  -- Name: lhaUpPtr  -- XML id: _24904
            // IGNORED: Field  -- Name: decayHandlePtr  -- XML id: _24905
            // IGNORED: Field  -- Name: handledParticles  -- XML id: _24906
            // IGNORED: Field  -- Name: userHooksPtr  -- XML id: _24907
            // IGNORED: Field  -- Name: hasUserHooks  -- XML id: _24908
            // IGNORED: Field  -- Name: doVetoProcess  -- XML id: _24909
            // IGNORED: Field  -- Name: doVetoPartons  -- XML id: _24910
            // IGNORED: Field  -- Name: retryPartonLevel  -- XML id: _24911
            // IGNORED: Field  -- Name: beamShapePtr  -- XML id: _24912
            // IGNORED: Field  -- Name: useNewBeamShape  -- XML id: _24913
            // IGNORED: Field  -- Name: doMomentumSpread  -- XML id: _24914
            // IGNORED: Field  -- Name: doVertexSpread  -- XML id: _24915
            // IGNORED: Field  -- Name: sigmaPtrs  -- XML id: _24916
            // IGNORED: Field  -- Name: resonancePtrs  -- XML id: _24917
            // IGNORED: Field  -- Name: timesDecPtr  -- XML id: _24918
            // IGNORED: Field  -- Name: timesPtr  -- XML id: _24919
            // IGNORED: Field  -- Name: spacePtr  -- XML id: _24920
            // IGNORED: Field  -- Name: useNewTimes  -- XML id: _24921
            // IGNORED: Field  -- Name: useNewSpace  -- XML id: _24922
            // IGNORED: Field  -- Name: processLevel  -- XML id: _24923
            // IGNORED: Field  -- Name: partonLevel  -- XML id: _24924
            // IGNORED: Field  -- Name: trialPartonLevel  -- XML id: _24925
            // IGNORED: Field  -- Name: hasMergingHooks  -- XML id: _24926
            // IGNORED: Field  -- Name: hasOwnMergingHooks  -- XML id: _24927
            // IGNORED: Field  -- Name: doMerging  -- XML id: _24928
            // IGNORED: Field  -- Name: hadronLevel  -- XML id: _24929
            // IGNORED: Field  -- Name: sigmaTot  -- XML id: _24930
            // IGNORED: Field  -- Name: rHadrons  -- XML id: _24931
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

        public:
            virtual Abstract_Pythia* pointerCopy_gambit() {std::cout << "Called virtual function" << std::endl;};
            virtual ~Abstract_Pythia() {};
    };
}
#pragma GCC diagnostic pop


#endif /* __ABSTRACT_PYTHIA_H__ */
