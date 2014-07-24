#ifndef __GAMBIT_WRAPPER_HIST_H__
#define __GAMBIT_WRAPPER_HIST_H__

#include "GAMBIT_wrapper_WrapperBase.h"
#include "abstract_Hist.h"
#include <vector>
#include <ostream>
#include <string>


// Factory function pointers to be filled by dynamic loading
Pythia8::Abstract_Hist* (*Factory_Hist_0)() = NULL;
Pythia8::Abstract_Hist* (*Factory_Hist_1)(std::string, int, double, double) = NULL;
Pythia8::Abstract_Hist* (*Factory_Hist_2)(std::string, const Pythia8::Abstract_Hist&) = NULL;

class Hist_gambit : public WrapperBase<Pythia8::Abstract_Hist>
{
    public:
        // Member variables: 

        // Member functions: 
        Hist_gambit operator=(const WrapperBase< Pythia8::Abstract_Hist > h)        {
            return BEptr->operator=(*h.BEptr);
        }

        void book(std::string titleIn, int nBinIn, double xMinIn, double xMaxIn)        {
            BEptr->book(titleIn, nBinIn, xMinIn, xMaxIn);
        }

        void name(std::string titleIn)        {
            BEptr->name(titleIn);
        }

        void null()        {
            BEptr->null();
        }

        void fill(double x, double w)        {
            BEptr->fill(x, w);
        }

        const void table(std::ostream& os, bool printOverUnder, bool xMidBin) const
        {
            BEptr->table(os, printOverUnder, xMidBin);
        }

        const void table(std::string fileName, bool printOverUnder, bool xMidBin) const
        {
            BEptr->table(fileName, printOverUnder, xMidBin);
        }

        const double getBinContent(int iBin) const
        {
            return BEptr->getBinContent(iBin);
        }

        const int getEntries() const
        {
            return BEptr->getEntries();
        }

        const bool sameSize(const WrapperBase< Pythia8::Abstract_Hist > h) const
        {
            return BEptr->sameSize(*h.BEptr);
        }

        void takeLog(bool tenLog)        {
            BEptr->takeLog(tenLog);
        }

        void takeSqrt()        {
            BEptr->takeSqrt();
        }

        Hist_gambit operator+=(const WrapperBase< Pythia8::Abstract_Hist > h)        {
            return BEptr->operator+=(*h.BEptr);
        }

        Hist_gambit operator-=(const WrapperBase< Pythia8::Abstract_Hist > h)        {
            return BEptr->operator-=(*h.BEptr);
        }

        Hist_gambit operator*=(const WrapperBase< Pythia8::Abstract_Hist > h)        {
            return BEptr->operator*=(*h.BEptr);
        }

        Hist_gambit operator/=(const WrapperBase< Pythia8::Abstract_Hist > h)        {
            return BEptr->operator/=(*h.BEptr);
        }

        Hist_gambit operator+=(double f)        {
            return BEptr->operator+=(f);
        }

        Hist_gambit operator-=(double f)        {
            return BEptr->operator-=(f);
        }

        Hist_gambit operator*=(double f)        {
            return BEptr->operator*=(f);
        }

        Hist_gambit operator/=(double f)        {
            return BEptr->operator/=(f);
        }


        // Wrappers for original constructors: 
        Hist_gambit() :
            WrapperBase<Pythia8::Abstract_Hist>( Factory_Hist_0(), false )
        {
        }

        Hist_gambit(std::string titleIn, int nBinIn, double xMinIn, double xMaxIn) :
            WrapperBase<Pythia8::Abstract_Hist>( Factory_Hist_1(titleIn, nBinIn, xMinIn, xMaxIn), false )
        {
        }

        Hist_gambit(std::string titleIn, const Hist_gambit h) :
            WrapperBase<Pythia8::Abstract_Hist>( Factory_Hist_2(titleIn, *h.BEptr), false )
        {
        }

        // Special pointer-based constructor: 
        Hist_gambit(Pythia8::Abstract_Hist* in) :
            WrapperBase<Pythia8::Abstract_Hist>( in, false )
        {
        }

        // Copy constructor: 
        Hist_gambit(const Hist_gambit& in) :
            WrapperBase<Pythia8::Abstract_Hist>( in.BEptr->pointerCopy_gambit(), false )
        {
        }
};

#endif /* __GAMBIT_WRAPPER_HIST_H__ */
