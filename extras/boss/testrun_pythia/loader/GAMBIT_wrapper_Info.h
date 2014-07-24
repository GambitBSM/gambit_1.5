#ifndef __GAMBIT_WRAPPER_INFO_H__
#define __GAMBIT_WRAPPER_INFO_H__

#include "GAMBIT_wrapper_WrapperBase.h"
#include <map>
#include <string>
#include <ostream>
#include "abstract_Info.h"
#include <vector>


// Factory function pointers to be filled by dynamic loading
Pythia8::Abstract_Info* (*Factory_Info_0)() = NULL;

class Info_gambit : public WrapperBase<Pythia8::Abstract_Info>
{
    public:
        // Member variables: 

        // Member functions: 
        const void list(std::ostream& os) const
        {
            BEptr->list(os);
        }

        const int idA() const
        {
            return BEptr->idA();
        }

        const int idB() const
        {
            return BEptr->idB();
        }

        const double pzA() const
        {
            return BEptr->pzA();
        }

        const double pzB() const
        {
            return BEptr->pzB();
        }

        const double eA() const
        {
            return BEptr->eA();
        }

        const double eB() const
        {
            return BEptr->eB();
        }

        const double mA() const
        {
            return BEptr->mA();
        }

        const double mB() const
        {
            return BEptr->mB();
        }

        const double eCM() const
        {
            return BEptr->eCM();
        }

        const double s() const
        {
            return BEptr->s();
        }

        const bool tooLowPTmin() const
        {
            return BEptr->tooLowPTmin();
        }

        const std::string name() const
        {
            return BEptr->name();
        }

        const int code() const
        {
            return BEptr->code();
        }

        const int nFinal() const
        {
            return BEptr->nFinal();
        }

        const bool isResolved() const
        {
            return BEptr->isResolved();
        }

        const bool isDiffractiveA() const
        {
            return BEptr->isDiffractiveA();
        }

        const bool isDiffractiveB() const
        {
            return BEptr->isDiffractiveB();
        }

        const bool isDiffractiveC() const
        {
            return BEptr->isDiffractiveC();
        }

        const bool isNonDiffractive() const
        {
            return BEptr->isNonDiffractive();
        }

        const bool isMinBias() const
        {
            return BEptr->isMinBias();
        }

        const bool isLHA() const
        {
            return BEptr->isLHA();
        }

        const bool atEndOfFile() const
        {
            return BEptr->atEndOfFile();
        }

        const bool hasSub(int i) const
        {
            return BEptr->hasSub(i);
        }

        const std::string nameSub(int i) const
        {
            return BEptr->nameSub(i);
        }

        const int codeSub(int i) const
        {
            return BEptr->codeSub(i);
        }

        const int nFinalSub(int i) const
        {
            return BEptr->nFinalSub(i);
        }

        const int id1(int i) const
        {
            return BEptr->id1(i);
        }

        const int id2(int i) const
        {
            return BEptr->id2(i);
        }

        const double x1(int i) const
        {
            return BEptr->x1(i);
        }

        const double x2(int i) const
        {
            return BEptr->x2(i);
        }

        const double y(int i) const
        {
            return BEptr->y(i);
        }

        const double tau(int i) const
        {
            return BEptr->tau(i);
        }

        const int id1pdf(int i) const
        {
            return BEptr->id1pdf(i);
        }

        const int id2pdf(int i) const
        {
            return BEptr->id2pdf(i);
        }

        const double x1pdf(int i) const
        {
            return BEptr->x1pdf(i);
        }

        const double x2pdf(int i) const
        {
            return BEptr->x2pdf(i);
        }

        const double pdf1(int i) const
        {
            return BEptr->pdf1(i);
        }

        const double pdf2(int i) const
        {
            return BEptr->pdf2(i);
        }

        const double QFac(int i) const
        {
            return BEptr->QFac(i);
        }

        const double Q2Fac(int i) const
        {
            return BEptr->Q2Fac(i);
        }

        const bool isValence1() const
        {
            return BEptr->isValence1();
        }

        const bool isValence2() const
        {
            return BEptr->isValence2();
        }

        const double alphaS(int i) const
        {
            return BEptr->alphaS(i);
        }

        const double alphaEM(int i) const
        {
            return BEptr->alphaEM(i);
        }

        const double QRen(int i) const
        {
            return BEptr->QRen(i);
        }

        const double Q2Ren(int i) const
        {
            return BEptr->Q2Ren(i);
        }

        const double scalup(int i) const
        {
            return BEptr->scalup(i);
        }

        const double mHat(int i) const
        {
            return BEptr->mHat(i);
        }

        const double sHat(int i) const
        {
            return BEptr->sHat(i);
        }

        const double tHat(int i) const
        {
            return BEptr->tHat(i);
        }

        const double uHat(int i) const
        {
            return BEptr->uHat(i);
        }

        const double pTHat(int i) const
        {
            return BEptr->pTHat(i);
        }

        const double pT2Hat(int i) const
        {
            return BEptr->pT2Hat(i);
        }

        const double m3Hat(int i) const
        {
            return BEptr->m3Hat(i);
        }

        const double m4Hat(int i) const
        {
            return BEptr->m4Hat(i);
        }

        const double thetaHat(int i) const
        {
            return BEptr->thetaHat(i);
        }

        const double phiHat(int i) const
        {
            return BEptr->phiHat(i);
        }

        const double weight() const
        {
            return BEptr->weight();
        }

        const double weightSum() const
        {
            return BEptr->weightSum();
        }

        const double lhaStrategy() const
        {
            return BEptr->lhaStrategy();
        }

        const int nISR() const
        {
            return BEptr->nISR();
        }

        const int nFSRinProc() const
        {
            return BEptr->nFSRinProc();
        }

        const int nFSRinRes() const
        {
            return BEptr->nFSRinRes();
        }

        const double pTmaxMPI() const
        {
            return BEptr->pTmaxMPI();
        }

        const double pTmaxISR() const
        {
            return BEptr->pTmaxISR();
        }

        const double pTmaxFSR() const
        {
            return BEptr->pTmaxFSR();
        }

        const double pTnow() const
        {
            return BEptr->pTnow();
        }

        const double a0MPI() const
        {
            return BEptr->a0MPI();
        }

        const double bMPI() const
        {
            return BEptr->bMPI();
        }

        const double enhanceMPI() const
        {
            return BEptr->enhanceMPI();
        }

        const double eMPI(int i) const
        {
            return BEptr->eMPI(i);
        }

        const int nMPI() const
        {
            return BEptr->nMPI();
        }

        const int codeMPI(int i) const
        {
            return BEptr->codeMPI(i);
        }

        const double pTMPI(int i) const
        {
            return BEptr->pTMPI(i);
        }

        const int iAMPI(int i) const
        {
            return BEptr->iAMPI(i);
        }

        const int iBMPI(int i) const
        {
            return BEptr->iBMPI(i);
        }

        std::vector<int, std::allocator<int> > codesHard()        {
            return BEptr->codesHard();
        }

        std::string nameProc(int i)        {
            return BEptr->nameProc(i);
        }

        long int nTried(int i)        {
            return BEptr->nTried(i);
        }

        long int nSelected(int i)        {
            return BEptr->nSelected(i);
        }

        long int nAccepted(int i)        {
            return BEptr->nAccepted(i);
        }

        double sigmaGen(int i)        {
            return BEptr->sigmaGen(i);
        }

        double sigmaErr(int i)        {
            return BEptr->sigmaErr(i);
        }

        const int getCounter(int i) const
        {
            return BEptr->getCounter(i);
        }

        void setCounter(int i, int value)        {
            BEptr->setCounter(i, value);
        }

        void addCounter(int i, int value)        {
            BEptr->addCounter(i, value);
        }

        void errorReset()        {
            BEptr->errorReset();
        }

        void errorMsg(std::string messageIn, std::string extraIn, bool showAlways, std::ostream& os)        {
            BEptr->errorMsg(messageIn, extraIn, showAlways, os);
        }

        int errorTotalNumber()        {
            return BEptr->errorTotalNumber();
        }

        void errorStatistics(std::ostream& os)        {
            BEptr->errorStatistics(os);
        }

        void setTooLowPTmin(bool lowPTminIn)        {
            BEptr->setTooLowPTmin(lowPTminIn);
        }

        void setValence(bool isVal1In, bool isVal2In)        {
            BEptr->setValence(isVal1In, isVal2In);
        }

        void hasHistory(bool hasHistoryIn)        {
            BEptr->hasHistory(hasHistoryIn);
        }

        bool hasHistory()        {
            return BEptr->hasHistory();
        }

        void zNowISR(double zNowIn)        {
            BEptr->zNowISR(zNowIn);
        }

        double zNowISR()        {
            return BEptr->zNowISR();
        }

        void pT2NowISR(double pT2NowIn)        {
            BEptr->pT2NowISR(pT2NowIn);
        }

        double pT2NowISR()        {
            return BEptr->pT2NowISR();
        }

        const double getWeightCKKWL() const
        {
            return BEptr->getWeightCKKWL();
        }

        void setWeightCKKWL(double weightIn)        {
            BEptr->setWeightCKKWL(weightIn);
        }

        const double mergingWeight() const
        {
            return BEptr->mergingWeight();
        }

        const double mergingWeightNLO() const
        {
            return BEptr->mergingWeightNLO();
        }

        const double getWeightFIRST() const
        {
            return BEptr->getWeightFIRST();
        }

        void setWeightFIRST(double weightIn)        {
            BEptr->setWeightFIRST(weightIn);
        }

        std::string header(const std::string& key)        {
            return BEptr->header(key);
        }

        std::vector<std::string, std::allocator<std::string> > headerKeys()        {
            return BEptr->headerKeys();
        }


        // Wrappers for original constructors: 
        Info_gambit() :
            WrapperBase<Pythia8::Abstract_Info>( Factory_Info_0(), false )
        {
        }

        // Special pointer-based constructor: 
        Info_gambit(Pythia8::Abstract_Info* in) :
            WrapperBase<Pythia8::Abstract_Info>( in, false )
        {
        }

        // Copy constructor: 
        Info_gambit(const Info_gambit& in) :
            WrapperBase<Pythia8::Abstract_Info>( in.BEptr->pointerCopy_gambit(), false )
        {
        }

        // Assignment operator: 
        Info_gambit& operator=(const Info_gambit& in)
        {
            if (this != &in) { BEptr->pointerAssign_gambit(in.BEptr); }
        }
};

#endif /* __GAMBIT_WRAPPER_INFO_H__ */
