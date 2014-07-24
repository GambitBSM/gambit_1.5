#ifndef __GAMBIT_WRAPPER_PYTHIA_H__
#define __GAMBIT_WRAPPER_PYTHIA_H__

#include "GAMBIT_wrapper_WrapperBase.h"
#include "abstract_Pythia.h"
#include <string>
#include <istream>
#include "GAMBIT_wrapper_Event.h"
#include "GAMBIT_wrapper_Vec4.h"
#include "GAMBIT_wrapper_Info.h"
#include <ostream>
#include <vector>


// Factory function pointers to be filled by dynamic loading
Pythia8::Abstract_Pythia* (*Factory_Pythia_0)(std::string, bool) = NULL;

class Pythia_gambit : public WrapperBase<Pythia8::Abstract_Pythia>
{
    public:
        // Member variables: 
        Event_gambit process;
        Event_gambit event;
        Info_gambit info;

        // Member functions: 
        bool readString(std::string line, bool warn)        {
            return BEptr->readString(line, warn);
        }

        bool readFile(std::string fileName, bool warn, int subrun)        {
            return BEptr->readFile(fileName, warn, subrun);
        }

        bool readFile(std::string fileName, int subrun)        {
            return BEptr->readFile(fileName, subrun);
        }

        bool readFile(std::istream& is, bool warn, int subrun)        {
            return BEptr->readFile(is, warn, subrun);
        }

        bool readFile(std::istream& is, int subrun)        {
            return BEptr->readFile(is, subrun);
        }

        bool init()        {
            return BEptr->init();
        }

        bool init(int idAin, int idBin, double eCMin)        {
            return BEptr->init(idAin, idBin, eCMin);
        }

        bool init(int idAin, int idBin, double eAin, double eBin)        {
            return BEptr->init(idAin, idBin, eAin, eBin);
        }

        bool init(int idAin, int idBin, double pxAin, double pyAin, double pzAin, double pxBin, double pyBin, double pzBin)        {
            return BEptr->init(idAin, idBin, pxAin, pyAin, pzAin, pxBin, pyBin, pzBin);
        }

        bool init(std::string LesHouchesEventFile, bool skipInit)        {
            return BEptr->init(LesHouchesEventFile, skipInit);
        }

        bool next()        {
            return BEptr->next();
        }

        int forceTimeShower(int iBeg, int iEnd, double pTmax, int nBranchMax)        {
            return BEptr->forceTimeShower(iBeg, iEnd, pTmax, nBranchMax);
        }

        bool forceHadronLevel(bool findJunctions)        {
            return BEptr->forceHadronLevel(findJunctions);
        }

        bool moreDecays()        {
            return BEptr->moreDecays();
        }

        bool forceRHadronDecays()        {
            return BEptr->forceRHadronDecays();
        }

        void LHAeventList(std::ostream& os)        {
            BEptr->LHAeventList(os);
        }

        bool LHAeventSkip(int nSkip)        {
            return BEptr->LHAeventSkip(nSkip);
        }

        void stat()        {
            BEptr->stat();
        }

        void statistics(bool all, bool reset)        {
            BEptr->statistics(all, reset);
        }

        bool flag(std::string key)        {
            return BEptr->flag(key);
        }

        int mode(std::string key)        {
            return BEptr->mode(key);
        }

        double parm(std::string key)        {
            return BEptr->parm(key);
        }

        std::string word(std::string key)        {
            return BEptr->word(key);
        }


        // Wrappers for original constructors: 
        Pythia_gambit(std::string xmlDir, bool printBanner) :
            WrapperBase<Pythia8::Abstract_Pythia>( Factory_Pythia_0(xmlDir, printBanner), false ),
            process(&(BEptr->process_ref_gambit())),
            event(&(BEptr->event_ref_gambit())),
            info(&(BEptr->info_ref_gambit()))
        {
            process._set_member_variable(true);
            event._set_member_variable(true);
            info._set_member_variable(true);
        }

        // Special pointer-based constructor: 
        Pythia_gambit(Pythia8::Abstract_Pythia* in) :
            WrapperBase<Pythia8::Abstract_Pythia>( in, false ),
            process(&(BEptr->process_ref_gambit())),
            event(&(BEptr->event_ref_gambit())),
            info(&(BEptr->info_ref_gambit()))
        {
            process._set_member_variable(true);
            event._set_member_variable(true);
            info._set_member_variable(true);
        }

        // Copy constructor: 
        Pythia_gambit(const Pythia_gambit& in) :
            WrapperBase<Pythia8::Abstract_Pythia>( in.BEptr->pointerCopy_gambit(), false ),
            process(&(BEptr->process_ref_gambit())),
            event(&(BEptr->event_ref_gambit())),
            info(&(BEptr->info_ref_gambit()))
        {
            process._set_member_variable(true);
            event._set_member_variable(true);
            info._set_member_variable(true);
        }
};

#endif /* __GAMBIT_WRAPPER_PYTHIA_H__ */
