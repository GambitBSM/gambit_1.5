#ifndef __GAMBIT_WRAPPER_PARTICLE_H__
#define __GAMBIT_WRAPPER_PARTICLE_H__

#include "abstract_Particle.h"
#include "GAMBIT_wrapper_WrapperBase.h"
#include "GAMBIT_wrapper_Vec4.h"
#include <vector>
#include <string>


// Factory function pointers to be filled by dynamic loading
Pythia8::Abstract_Particle* (*Factory_Particle_0)() = NULL;
Pythia8::Abstract_Particle* (*Factory_Particle_1)(int, int, int, int, int, int, int, int, double, double, double, double, double, double, double) = NULL;
Pythia8::Abstract_Particle* (*Factory_Particle_2)(int, int, int, int, int, int, int, int, Pythia8::Abstract_Vec4&, double, double, double) = NULL;

class Particle_gambit : public WrapperBase<Pythia8::Abstract_Particle>
{
    public:
        // Member variables: 

        // Member functions: 
        Particle_gambit operator=(const WrapperBase< Pythia8::Abstract_Particle > pt)        {
            return BEptr->operator=(*pt.BEptr);
        }

        void setEvtPtr(WrapperBase< Pythia8::Abstract_Event >* evtPtrIn)        {
            BEptr->setEvtPtr((*evtPtrIn).BEptr);
        }

        void id(int idIn)        {
            BEptr->id(idIn);
        }

        void status(int statusIn)        {
            BEptr->status(statusIn);
        }

        void statusPos()        {
            BEptr->statusPos();
        }

        void statusNeg()        {
            BEptr->statusNeg();
        }

        void statusCode(int statusIn)        {
            BEptr->statusCode(statusIn);
        }

        void mother1(int mother1In)        {
            BEptr->mother1(mother1In);
        }

        void mother2(int mother2In)        {
            BEptr->mother2(mother2In);
        }

        void mothers(int mother1In, int mother2In)        {
            BEptr->mothers(mother1In, mother2In);
        }

        void daughter1(int daughter1In)        {
            BEptr->daughter1(daughter1In);
        }

        void daughter2(int daughter2In)        {
            BEptr->daughter2(daughter2In);
        }

        void daughters(int daughter1In, int daughter2In)        {
            BEptr->daughters(daughter1In, daughter2In);
        }

        void col(int colIn)        {
            BEptr->col(colIn);
        }

        void acol(int acolIn)        {
            BEptr->acol(acolIn);
        }

        void cols(int colIn, int acolIn)        {
            BEptr->cols(colIn, acolIn);
        }

        void p(WrapperBase< Pythia8::Abstract_Vec4 > pIn)        {
            BEptr->p(*pIn.BEptr);
        }

        void p(double pxIn, double pyIn, double pzIn, double eIn)        {
            BEptr->p(pxIn, pyIn, pzIn, eIn);
        }

        void px(double pxIn)        {
            BEptr->px(pxIn);
        }

        void py(double pyIn)        {
            BEptr->py(pyIn);
        }

        void pz(double pzIn)        {
            BEptr->pz(pzIn);
        }

        void e(double eIn)        {
            BEptr->e(eIn);
        }

        void m(double mIn)        {
            BEptr->m(mIn);
        }

        void scale(double scaleIn)        {
            BEptr->scale(scaleIn);
        }

        void pol(double polIn)        {
            BEptr->pol(polIn);
        }

        void vProd(WrapperBase< Pythia8::Abstract_Vec4 > vProdIn)        {
            BEptr->vProd(*vProdIn.BEptr);
        }

        void vProd(double xProdIn, double yProdIn, double zProdIn, double tProdIn)        {
            BEptr->vProd(xProdIn, yProdIn, zProdIn, tProdIn);
        }

        void xProd(double xProdIn)        {
            BEptr->xProd(xProdIn);
        }

        void yProd(double yProdIn)        {
            BEptr->yProd(yProdIn);
        }

        void zProd(double zProdIn)        {
            BEptr->zProd(zProdIn);
        }

        void tProd(double tProdIn)        {
            BEptr->tProd(tProdIn);
        }

        void tau(double tauIn)        {
            BEptr->tau(tauIn);
        }

        const int id() const
        {
            return BEptr->id();
        }

        const int status() const
        {
            return BEptr->status();
        }

        const int mother1() const
        {
            return BEptr->mother1();
        }

        const int mother2() const
        {
            return BEptr->mother2();
        }

        const int daughter1() const
        {
            return BEptr->daughter1();
        }

        const int daughter2() const
        {
            return BEptr->daughter2();
        }

        const int col() const
        {
            return BEptr->col();
        }

        const int acol() const
        {
            return BEptr->acol();
        }

        const Vec4_gambit p() const
        {
            return BEptr->p();
        }

        const double px() const
        {
            return BEptr->px();
        }

        const double py() const
        {
            return BEptr->py();
        }

        const double pz() const
        {
            return BEptr->pz();
        }

        const double e() const
        {
            return BEptr->e();
        }

        const double m() const
        {
            return BEptr->m();
        }

        const double scale() const
        {
            return BEptr->scale();
        }

        const double pol() const
        {
            return BEptr->pol();
        }

        const bool hasVertex() const
        {
            return BEptr->hasVertex();
        }

        const Vec4_gambit vProd() const
        {
            return BEptr->vProd();
        }

        const double xProd() const
        {
            return BEptr->xProd();
        }

        const double yProd() const
        {
            return BEptr->yProd();
        }

        const double zProd() const
        {
            return BEptr->zProd();
        }

        const double tProd() const
        {
            return BEptr->tProd();
        }

        const double tau() const
        {
            return BEptr->tau();
        }

        const int idAbs() const
        {
            return BEptr->idAbs();
        }

        const int statusAbs() const
        {
            return BEptr->statusAbs();
        }

        const bool isFinal() const
        {
            return BEptr->isFinal();
        }

        const bool isRescatteredIncoming() const
        {
            return BEptr->isRescatteredIncoming();
        }

        const double m2() const
        {
            return BEptr->m2();
        }

        const double mCalc() const
        {
            return BEptr->mCalc();
        }

        const double m2Calc() const
        {
            return BEptr->m2Calc();
        }

        const double eCalc() const
        {
            return BEptr->eCalc();
        }

        const double pT() const
        {
            return BEptr->pT();
        }

        const double pT2() const
        {
            return BEptr->pT2();
        }

        const double mT() const
        {
            return BEptr->mT();
        }

        const double mT2() const
        {
            return BEptr->mT2();
        }

        const double pAbs() const
        {
            return BEptr->pAbs();
        }

        const double pAbs2() const
        {
            return BEptr->pAbs2();
        }

        const double eT() const
        {
            return BEptr->eT();
        }

        const double eT2() const
        {
            return BEptr->eT2();
        }

        const double theta() const
        {
            return BEptr->theta();
        }

        const double phi() const
        {
            return BEptr->phi();
        }

        const double thetaXZ() const
        {
            return BEptr->thetaXZ();
        }

        const double pPos() const
        {
            return BEptr->pPos();
        }

        const double pNeg() const
        {
            return BEptr->pNeg();
        }

        const double y() const
        {
            return BEptr->y();
        }

        const double eta() const
        {
            return BEptr->eta();
        }

        const Vec4_gambit vDec() const
        {
            return BEptr->vDec();
        }

        const double xDec() const
        {
            return BEptr->xDec();
        }

        const double yDec() const
        {
            return BEptr->yDec();
        }

        const double zDec() const
        {
            return BEptr->zDec();
        }

        const double tDec() const
        {
            return BEptr->tDec();
        }

        const int index() const
        {
            return BEptr->index();
        }

        const int statusHepMC() const
        {
            return BEptr->statusHepMC();
        }

        const int iTopCopy() const
        {
            return BEptr->iTopCopy();
        }

        const int iBotCopy() const
        {
            return BEptr->iBotCopy();
        }

        const int iTopCopyId() const
        {
            return BEptr->iTopCopyId();
        }

        const int iBotCopyId() const
        {
            return BEptr->iBotCopyId();
        }

        const std::vector<int, std::allocator<int> > motherList() const
        {
            return BEptr->motherList();
        }

        const std::vector<int, std::allocator<int> > daughterList() const
        {
            return BEptr->daughterList();
        }

        const std::vector<int, std::allocator<int> > sisterList(bool traceTopBot) const
        {
            return BEptr->sisterList(traceTopBot);
        }

        const bool isAncestor(int iAncestor) const
        {
            return BEptr->isAncestor(iAncestor);
        }

        bool undoDecay()        {
            return BEptr->undoDecay();
        }

        const std::string name() const
        {
            return BEptr->name();
        }

        const std::string nameWithStatus(int maxLen) const
        {
            return BEptr->nameWithStatus(maxLen);
        }

        const int spinType() const
        {
            return BEptr->spinType();
        }

        const int chargeType() const
        {
            return BEptr->chargeType();
        }

        const double charge() const
        {
            return BEptr->charge();
        }

        const bool isCharged() const
        {
            return BEptr->isCharged();
        }

        const bool isNeutral() const
        {
            return BEptr->isNeutral();
        }

        const int colType() const
        {
            return BEptr->colType();
        }

        const double m0() const
        {
            return BEptr->m0();
        }

        const double mWidth() const
        {
            return BEptr->mWidth();
        }

        const double mMin() const
        {
            return BEptr->mMin();
        }

        const double mMax() const
        {
            return BEptr->mMax();
        }

        const double mSel() const
        {
            return BEptr->mSel();
        }

        const double constituentMass() const
        {
            return BEptr->constituentMass();
        }

        const double tau0() const
        {
            return BEptr->tau0();
        }

        const bool mayDecay() const
        {
            return BEptr->mayDecay();
        }

        const bool canDecay() const
        {
            return BEptr->canDecay();
        }

        const bool doExternalDecay() const
        {
            return BEptr->doExternalDecay();
        }

        const bool isResonance() const
        {
            return BEptr->isResonance();
        }

        const bool isVisible() const
        {
            return BEptr->isVisible();
        }

        const bool isLepton() const
        {
            return BEptr->isLepton();
        }

        const bool isQuark() const
        {
            return BEptr->isQuark();
        }

        const bool isGluon() const
        {
            return BEptr->isGluon();
        }

        const bool isDiquark() const
        {
            return BEptr->isDiquark();
        }

        const bool isParton() const
        {
            return BEptr->isParton();
        }

        const bool isHadron() const
        {
            return BEptr->isHadron();
        }

        void rescale3(double fac)        {
            BEptr->rescale3(fac);
        }

        void rescale4(double fac)        {
            BEptr->rescale4(fac);
        }

        void rescale5(double fac)        {
            BEptr->rescale5(fac);
        }

        void rot(double thetaIn, double phiIn)        {
            BEptr->rot(thetaIn, phiIn);
        }

        void bst(double betaX, double betaY, double betaZ)        {
            BEptr->bst(betaX, betaY, betaZ);
        }

        void bst(double betaX, double betaY, double betaZ, double gamma)        {
            BEptr->bst(betaX, betaY, betaZ, gamma);
        }

        void bst(const WrapperBase< Pythia8::Abstract_Vec4 > pBst)        {
            BEptr->bst(*pBst.BEptr);
        }

        void bst(const WrapperBase< Pythia8::Abstract_Vec4 > pBst, double mBst)        {
            BEptr->bst(*pBst.BEptr, mBst);
        }

        void bstback(const WrapperBase< Pythia8::Abstract_Vec4 > pBst)        {
            BEptr->bstback(*pBst.BEptr);
        }

        void bstback(const WrapperBase< Pythia8::Abstract_Vec4 > pBst, double mBst)        {
            BEptr->bstback(*pBst.BEptr, mBst);
        }

        void offsetHistory(int minMother, int addMother, int minDaughter, int addDaughter)        {
            BEptr->offsetHistory(minMother, addMother, minDaughter, addDaughter);
        }

        void offsetCol(int addCol)        {
            BEptr->offsetCol(addCol);
        }


        // Wrappers for original constructors: 
        Particle_gambit() :
            WrapperBase<Pythia8::Abstract_Particle>( Factory_Particle_0(), false )
        {
        }

        Particle_gambit(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, double pxIn, double pyIn, double pzIn, double eIn, double mIn, double scaleIn, double polIn) :
            WrapperBase<Pythia8::Abstract_Particle>( Factory_Particle_1(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, pxIn, pyIn, pzIn, eIn, mIn, scaleIn, polIn), false )
        {
        }

        Particle_gambit(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, Vec4_gambit pIn, double mIn, double scaleIn, double polIn) :
            WrapperBase<Pythia8::Abstract_Particle>( Factory_Particle_2(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, *pIn.BEptr, mIn, scaleIn, polIn), false )
        {
        }

        // Special pointer-based constructor: 
        Particle_gambit(Pythia8::Abstract_Particle* in) :
            WrapperBase<Pythia8::Abstract_Particle>( in, false )
        {
        }

        // Copy constructor: 
        Particle_gambit(const Particle_gambit& in) :
            WrapperBase<Pythia8::Abstract_Particle>( in.BEptr->pointerCopy_gambit(), false )
        {
        }
};

#endif /* __GAMBIT_WRAPPER_PARTICLE_H__ */
