#ifndef __GAMBIT_WRAPPER_VEC4_H__
#define __GAMBIT_WRAPPER_VEC4_H__

#include "GAMBIT_wrapper_WrapperBase.h"
#include "abstract_Vec4.h"


// Factory function pointers to be filled by dynamic loading
Pythia8::Abstract_Vec4* (*Factory_Vec4_0)(double, double, double, double) = NULL;

class Vec4_gambit : public WrapperBase<Pythia8::Abstract_Vec4>
{
    public:
        // Member variables: 

        // Member functions: 
        Vec4_gambit operator=(const WrapperBase< Pythia8::Abstract_Vec4 > v)        {
            return BEptr->operator=(*v.BEptr);
        }

        Vec4_gambit operator=(double value)        {
            return BEptr->operator=(value);
        }

        void reset()        {
            BEptr->reset();
        }

        void p(double xIn, double yIn, double zIn, double tIn)        {
            BEptr->p(xIn, yIn, zIn, tIn);
        }

        void p(WrapperBase< Pythia8::Abstract_Vec4 > pIn)        {
            BEptr->p(*pIn.BEptr);
        }

        void px(double xIn)        {
            BEptr->px(xIn);
        }

        void py(double yIn)        {
            BEptr->py(yIn);
        }

        void pz(double zIn)        {
            BEptr->pz(zIn);
        }

        void e(double tIn)        {
            BEptr->e(tIn);
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

        const double mCalc() const
        {
            return BEptr->mCalc();
        }

        const double m2Calc() const
        {
            return BEptr->m2Calc();
        }

        const double pT() const
        {
            return BEptr->pT();
        }

        const double pT2() const
        {
            return BEptr->pT2();
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

        const double rap() const
        {
            return BEptr->rap();
        }

        const double eta() const
        {
            return BEptr->eta();
        }

        void rescale3(double fac)        {
            BEptr->rescale3(fac);
        }

        void rescale4(double fac)        {
            BEptr->rescale4(fac);
        }

        void flip3()        {
            BEptr->flip3();
        }

        void flip4()        {
            BEptr->flip4();
        }

        void rot(double thetaIn, double phiIn)        {
            BEptr->rot(thetaIn, phiIn);
        }

        void rotaxis(double phiIn, double nx, double ny, double nz)        {
            BEptr->rotaxis(phiIn, nx, ny, nz);
        }

        void rotaxis(double phiIn, const WrapperBase< Pythia8::Abstract_Vec4 > n)        {
            BEptr->rotaxis(phiIn, *n.BEptr);
        }

        void bst(double betaX, double betaY, double betaZ)        {
            BEptr->bst(betaX, betaY, betaZ);
        }

        void bst(double betaX, double betaY, double betaZ, double gamma)        {
            BEptr->bst(betaX, betaY, betaZ, gamma);
        }

        void bst(const WrapperBase< Pythia8::Abstract_Vec4 > pIn)        {
            BEptr->bst(*pIn.BEptr);
        }

        void bst(const WrapperBase< Pythia8::Abstract_Vec4 > pIn, double mIn)        {
            BEptr->bst(*pIn.BEptr, mIn);
        }

        void bstback(const WrapperBase< Pythia8::Abstract_Vec4 > pIn)        {
            BEptr->bstback(*pIn.BEptr);
        }

        void bstback(const WrapperBase< Pythia8::Abstract_Vec4 > pIn, double mIn)        {
            BEptr->bstback(*pIn.BEptr, mIn);
        }

        Vec4_gambit operator-()        {
            return BEptr->operator-();
        }

        Vec4_gambit operator+=(const WrapperBase< Pythia8::Abstract_Vec4 > v)        {
            return BEptr->operator+=(*v.BEptr);
        }

        Vec4_gambit operator-=(const WrapperBase< Pythia8::Abstract_Vec4 > v)        {
            return BEptr->operator-=(*v.BEptr);
        }

        Vec4_gambit operator*=(double f)        {
            return BEptr->operator*=(f);
        }

        Vec4_gambit operator/=(double f)        {
            return BEptr->operator/=(f);
        }


        // Wrappers for original constructors: 
        Vec4_gambit(double xIn, double yIn, double zIn, double tIn) :
            WrapperBase<Pythia8::Abstract_Vec4>( Factory_Vec4_0(xIn, yIn, zIn, tIn), false )
        {
        }

        // Special pointer-based constructor: 
        Vec4_gambit(Pythia8::Abstract_Vec4* in) :
            WrapperBase<Pythia8::Abstract_Vec4>( in, false )
        {
        }

        // Copy constructor: 
        Vec4_gambit(const Vec4_gambit& in) :
            WrapperBase<Pythia8::Abstract_Vec4>( in.BEptr->pointerCopy_gambit(), false )
        {
        }
};

#endif /* __GAMBIT_WRAPPER_VEC4_H__ */
