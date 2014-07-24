#ifndef __ABSTRACT_INFO_H__
#define __ABSTRACT_INFO_H__

#include <iostream>  // FOR DEBUG: Allow virtual member functions to print a warning if executed.

#include "forward_decls_abstract_classes.h"
#include <vector>
#include <ostream>
#include <string>
#include <map>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wreturn-type"
namespace Pythia8
{
    class Abstract_Info
    {
        private:
            // IGNORED: Variable  -- Name: TIMESTOPRINT  -- XML id: _21792
            // IGNORED: Variable  -- Name: CONVERTMB2PB  -- XML id: _21793
            // IGNORED: Field  -- Name: idASave  -- XML id: _21794
            // IGNORED: Field  -- Name: idBSave  -- XML id: _21795
            // IGNORED: Field  -- Name: pzASave  -- XML id: _21796
            // IGNORED: Field  -- Name: eASave  -- XML id: _21797
            // IGNORED: Field  -- Name: mASave  -- XML id: _21798
            // IGNORED: Field  -- Name: pzBSave  -- XML id: _21799
            // IGNORED: Field  -- Name: eBSave  -- XML id: _21800
            // IGNORED: Field  -- Name: mBSave  -- XML id: _21801
            // IGNORED: Field  -- Name: eCMSave  -- XML id: _21802
            // IGNORED: Field  -- Name: sSave  -- XML id: _21803
            // IGNORED: Field  -- Name: lowPTmin  -- XML id: _21804
            // IGNORED: Field  -- Name: nTry  -- XML id: _21805
            // IGNORED: Field  -- Name: nSel  -- XML id: _21806
            // IGNORED: Field  -- Name: nAcc  -- XML id: _21807
            // IGNORED: Field  -- Name: sigGen  -- XML id: _21808
            // IGNORED: Field  -- Name: sigErr  -- XML id: _21809
            // IGNORED: Field  -- Name: wtAccSum  -- XML id: _21810
            // IGNORED: Field  -- Name: procNameM  -- XML id: _21811
            // IGNORED: Field  -- Name: nTryM  -- XML id: _21812
            // IGNORED: Field  -- Name: nSelM  -- XML id: _21813
            // IGNORED: Field  -- Name: nAccM  -- XML id: _21814
            // IGNORED: Field  -- Name: sigGenM  -- XML id: _21815
            // IGNORED: Field  -- Name: sigErrM  -- XML id: _21816
            // IGNORED: Field  -- Name: lhaStrategySave  -- XML id: _21817
            // IGNORED: Field  -- Name: a0MPISave  -- XML id: _21818
            // IGNORED: Field  -- Name: isRes  -- XML id: _21819
            // IGNORED: Field  -- Name: isDiffA  -- XML id: _21820
            // IGNORED: Field  -- Name: isDiffB  -- XML id: _21821
            // IGNORED: Field  -- Name: isDiffC  -- XML id: _21822
            // IGNORED: Field  -- Name: isND  -- XML id: _21823
            // IGNORED: Field  -- Name: isLH  -- XML id: _21824
            // IGNORED: Field  -- Name: hasSubSave  -- XML id: _21825
            // IGNORED: Field  -- Name: bIsSet  -- XML id: _21826
            // IGNORED: Field  -- Name: evolIsSet  -- XML id: _21827
            // IGNORED: Field  -- Name: atEOF  -- XML id: _21828
            // IGNORED: Field  -- Name: isVal1  -- XML id: _21829
            // IGNORED: Field  -- Name: isVal2  -- XML id: _21830
            // IGNORED: Field  -- Name: hasHistorySave  -- XML id: _21831
            // IGNORED: Field  -- Name: codeSave  -- XML id: _21832
            // IGNORED: Field  -- Name: codeSubSave  -- XML id: _21833
            // IGNORED: Field  -- Name: nFinalSave  -- XML id: _21834
            // IGNORED: Field  -- Name: nFinalSubSave  -- XML id: _21835
            // IGNORED: Field  -- Name: nTotal  -- XML id: _21836
            // IGNORED: Field  -- Name: id1Save  -- XML id: _21837
            // IGNORED: Field  -- Name: id2Save  -- XML id: _21838
            // IGNORED: Field  -- Name: id1pdfSave  -- XML id: _21839
            // IGNORED: Field  -- Name: id2pdfSave  -- XML id: _21840
            // IGNORED: Field  -- Name: nMPISave  -- XML id: _21841
            // IGNORED: Field  -- Name: nISRSave  -- XML id: _21842
            // IGNORED: Field  -- Name: nFSRinProcSave  -- XML id: _21843
            // IGNORED: Field  -- Name: nFSRinResSave  -- XML id: _21844
            // IGNORED: Field  -- Name: x1Save  -- XML id: _21845
            // IGNORED: Field  -- Name: x2Save  -- XML id: _21846
            // IGNORED: Field  -- Name: x1pdfSave  -- XML id: _21847
            // IGNORED: Field  -- Name: x2pdfSave  -- XML id: _21848
            // IGNORED: Field  -- Name: pdf1Save  -- XML id: _21849
            // IGNORED: Field  -- Name: pdf2Save  -- XML id: _21850
            // IGNORED: Field  -- Name: Q2FacSave  -- XML id: _21851
            // IGNORED: Field  -- Name: alphaEMSave  -- XML id: _21852
            // IGNORED: Field  -- Name: alphaSSave  -- XML id: _21853
            // IGNORED: Field  -- Name: Q2RenSave  -- XML id: _21854
            // IGNORED: Field  -- Name: scalupSave  -- XML id: _21855
            // IGNORED: Field  -- Name: sH  -- XML id: _21856
            // IGNORED: Field  -- Name: tH  -- XML id: _21857
            // IGNORED: Field  -- Name: uH  -- XML id: _21858
            // IGNORED: Field  -- Name: pTH  -- XML id: _21859
            // IGNORED: Field  -- Name: m3H  -- XML id: _21860
            // IGNORED: Field  -- Name: m4H  -- XML id: _21861
            // IGNORED: Field  -- Name: thetaH  -- XML id: _21862
            // IGNORED: Field  -- Name: phiH  -- XML id: _21863
            // IGNORED: Field  -- Name: weightSave  -- XML id: _21864
            // IGNORED: Field  -- Name: bMPISave  -- XML id: _21865
            // IGNORED: Field  -- Name: enhanceMPISave  -- XML id: _21866
            // IGNORED: Field  -- Name: pTmaxMPISave  -- XML id: _21867
            // IGNORED: Field  -- Name: pTmaxISRSave  -- XML id: _21868
            // IGNORED: Field  -- Name: pTmaxFSRSave  -- XML id: _21869
            // IGNORED: Field  -- Name: pTnowSave  -- XML id: _21870
            // IGNORED: Field  -- Name: zNowISRSave  -- XML id: _21871
            // IGNORED: Field  -- Name: pT2NowISRSave  -- XML id: _21872
            // IGNORED: Field  -- Name: nameSave  -- XML id: _21873
            // IGNORED: Field  -- Name: nameSubSave  -- XML id: _21874
            // IGNORED: Field  -- Name: codeMPISave  -- XML id: _21875
            // IGNORED: Field  -- Name: iAMPISave  -- XML id: _21876
            // IGNORED: Field  -- Name: iBMPISave  -- XML id: _21877
            // IGNORED: Field  -- Name: pTMPISave  -- XML id: _21878
            // IGNORED: Field  -- Name: eMPISave  -- XML id: _21879
            // IGNORED: Field  -- Name: counters  -- XML id: _21880
            // IGNORED: Field  -- Name: messages  -- XML id: _21881
            // IGNORED: Field  -- Name: headers  -- XML id: _21882
            // IGNORED: Field  -- Name: weightCKKWLSave  -- XML id: _21883
            // IGNORED: Field  -- Name: weightFIRSTSave  -- XML id: _21884
        public:

            virtual void list(std::ostream& os) const {std::cout << "Called virtual function" << std::endl;};

            virtual int idA() const {std::cout << "Called virtual function" << std::endl;};

            virtual int idB() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pzA() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pzB() const {std::cout << "Called virtual function" << std::endl;};

            virtual double eA() const {std::cout << "Called virtual function" << std::endl;};

            virtual double eB() const {std::cout << "Called virtual function" << std::endl;};

            virtual double mA() const {std::cout << "Called virtual function" << std::endl;};

            virtual double mB() const {std::cout << "Called virtual function" << std::endl;};

            virtual double eCM() const {std::cout << "Called virtual function" << std::endl;};

            virtual double s() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool tooLowPTmin() const {std::cout << "Called virtual function" << std::endl;};

            virtual std::string name() const {std::cout << "Called virtual function" << std::endl;};

            virtual int code() const {std::cout << "Called virtual function" << std::endl;};

            virtual int nFinal() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isResolved() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isDiffractiveA() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isDiffractiveB() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isDiffractiveC() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isNonDiffractive() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isMinBias() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isLHA() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool atEndOfFile() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool hasSub(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual std::string nameSub(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int codeSub(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int nFinalSub(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int id1(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int id2(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double x1(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double x2(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double y(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double tau(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int id1pdf(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int id2pdf(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double x1pdf(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double x2pdf(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double pdf1(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double pdf2(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double QFac(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double Q2Fac(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isValence1() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isValence2() const {std::cout << "Called virtual function" << std::endl;};

            virtual double alphaS(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double alphaEM(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double QRen(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double Q2Ren(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double scalup(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double mHat(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double sHat(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double tHat(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double uHat(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double pTHat(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double pT2Hat(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double m3Hat(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double m4Hat(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double thetaHat(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double phiHat(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double weight() const {std::cout << "Called virtual function" << std::endl;};

            virtual double weightSum() const {std::cout << "Called virtual function" << std::endl;};

            virtual double lhaStrategy() const {std::cout << "Called virtual function" << std::endl;};

            virtual int nISR() const {std::cout << "Called virtual function" << std::endl;};

            virtual int nFSRinProc() const {std::cout << "Called virtual function" << std::endl;};

            virtual int nFSRinRes() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pTmaxMPI() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pTmaxISR() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pTmaxFSR() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pTnow() const {std::cout << "Called virtual function" << std::endl;};

            virtual double a0MPI() const {std::cout << "Called virtual function" << std::endl;};

            virtual double bMPI() const {std::cout << "Called virtual function" << std::endl;};

            virtual double enhanceMPI() const {std::cout << "Called virtual function" << std::endl;};

            virtual double eMPI(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int nMPI() const {std::cout << "Called virtual function" << std::endl;};

            virtual int codeMPI(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual double pTMPI(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int iAMPI(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int iBMPI(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual std::vector<int, std::allocator<int> > codesHard() {std::cout << "Called virtual function" << std::endl;};

            virtual std::string nameProc(int i) {std::cout << "Called virtual function" << std::endl;};

            virtual long int nTried(int i) {std::cout << "Called virtual function" << std::endl;};

            virtual long int nSelected(int i) {std::cout << "Called virtual function" << std::endl;};

            virtual long int nAccepted(int i) {std::cout << "Called virtual function" << std::endl;};

            virtual double sigmaGen(int i) {std::cout << "Called virtual function" << std::endl;};

            virtual double sigmaErr(int i) {std::cout << "Called virtual function" << std::endl;};

            virtual int getCounter(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual void setCounter(int i, int value) {std::cout << "Called virtual function" << std::endl;};

            virtual void addCounter(int i, int value) {std::cout << "Called virtual function" << std::endl;};

            virtual void errorReset() {std::cout << "Called virtual function" << std::endl;};

            virtual void errorMsg(std::string messageIn, std::string extraIn, bool showAlways, std::ostream& os) {std::cout << "Called virtual function" << std::endl;};

            virtual int errorTotalNumber() {std::cout << "Called virtual function" << std::endl;};

            virtual void errorStatistics(std::ostream& os) {std::cout << "Called virtual function" << std::endl;};

            virtual void setTooLowPTmin(bool lowPTminIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void setValence(bool isVal1In, bool isVal2In) {std::cout << "Called virtual function" << std::endl;};

            virtual void hasHistory(bool hasHistoryIn) {std::cout << "Called virtual function" << std::endl;};

            virtual bool hasHistory() {std::cout << "Called virtual function" << std::endl;};

            virtual void zNowISR(double zNowIn) {std::cout << "Called virtual function" << std::endl;};

            virtual double zNowISR() {std::cout << "Called virtual function" << std::endl;};

            virtual void pT2NowISR(double pT2NowIn) {std::cout << "Called virtual function" << std::endl;};

            virtual double pT2NowISR() {std::cout << "Called virtual function" << std::endl;};

            virtual double getWeightCKKWL() const {std::cout << "Called virtual function" << std::endl;};

            virtual void setWeightCKKWL(double weightIn) {std::cout << "Called virtual function" << std::endl;};

            virtual double mergingWeight() const {std::cout << "Called virtual function" << std::endl;};

            virtual double mergingWeightNLO() const {std::cout << "Called virtual function" << std::endl;};

            virtual double getWeightFIRST() const {std::cout << "Called virtual function" << std::endl;};

            virtual void setWeightFIRST(double weightIn) {std::cout << "Called virtual function" << std::endl;};

            virtual std::string header(const std::string& key) {std::cout << "Called virtual function" << std::endl;};

            virtual std::vector<std::string, std::allocator<std::string> > headerKeys() {std::cout << "Called virtual function" << std::endl;};
        private:

            virtual void setBeamA(int idAin, double pzAin, double eAin, double mAin) {std::cout << "Called virtual function" << std::endl;};

            virtual void setBeamB(int idBin, double pzBin, double eBin, double mBin) {std::cout << "Called virtual function" << std::endl;};

            virtual void setECM(double eCMin) {std::cout << "Called virtual function" << std::endl;};

            virtual void clear() {std::cout << "Called virtual function" << std::endl;};

            virtual int sizeMPIarrays() const {std::cout << "Called virtual function" << std::endl;};

            virtual void resizeMPIarrays(int newSize) {std::cout << "Called virtual function" << std::endl;};

            virtual void setType(std::string nameIn, int codeIn, int nFinalIn, bool isNonDiffIn, bool isResolvedIn, bool isDiffractiveAin, bool isDiffractiveBin, bool isDiffractiveCin, bool isLHAin) {std::cout << "Called virtual function" << std::endl;};

            virtual void setSubType(int iDS, std::string nameSubIn, int codeSubIn, int nFinalSubIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void setPDFalpha(int iDS, int id1pdfIn, int id2pdfIn, double x1pdfIn, double x2pdfIn, double pdf1In, double pdf2In, double Q2FacIn, double alphaEMIn, double alphaSIn, double Q2RenIn, double scalupIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void setScalup(int iDS, double scalupIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void setKin(int iDS, int id1In, int id2In, double x1In, double x2In, double sHatIn, double tHatIn, double uHatIn, double pTHatIn, double m3HatIn, double m4HatIn, double thetaHatIn, double phiHatIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void setTypeMPI(int codeMPIIn, double pTMPIIn, int iAMPIIn, int iBMPIIn, double eMPIIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void sigmaReset() {std::cout << "Called virtual function" << std::endl;};

            virtual void setSigma(int i, std::string procNameIn, long int nTryIn, long int nSelIn, long int nAccIn, double sigGenIn, double sigErrIn, double wtAccSumIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void addSigma(int i, long int nTryIn, long int nSelIn, long int nAccIn, double sigGenIn, double sigErrIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void setImpact(double bMPIIn, double enhanceMPIIn, bool bIsSetIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void setPartEvolved(int nMPIIn, int nISRIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void setEvolution(double pTmaxMPIIn, double pTmaxISRIn, double pTmaxFSRIn, int nMPIIn, int nISRIn, int nFSRinProcIn, int nFSRinResIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void setPTnow(double pTnowIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void seta0MPI(double a0MPIin) {std::cout << "Called virtual function" << std::endl;};

            virtual void setEndOfFile(bool atEOFin) {std::cout << "Called virtual function" << std::endl;};

            virtual void setWeight(double weightIn, int lhaStrategyIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void setHeader(const std::string& key, const std::string& val) {std::cout << "Called virtual function" << std::endl;};

        public:
            virtual void pointerAssign_gambit(Abstract_Info* in) {std::cout << "Called virtual function" << std::endl;};
            virtual Abstract_Info* pointerCopy_gambit() {std::cout << "Called virtual function" << std::endl;};
            virtual ~Abstract_Info() {};
    };
}
#pragma GCC diagnostic pop


#endif /* __ABSTRACT_INFO_H__ */
