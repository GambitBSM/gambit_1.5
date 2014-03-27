
        public:
            bool* readString_GAMBIT(std::string line, bool warn)
            {
                return new bool(readString(line,warn));
            }
            bool* readFile_GAMBIT(std::string fileName, bool warn, int subrun)
            {
                return new bool(readFile(fileName,warn,subrun));
            }
            bool* readFile_GAMBIT(std::string fileName, int subrun)
            {
                return new bool(readFile(fileName,subrun));
            }
            bool* readFile_GAMBIT(std::istream& is, bool warn, int subrun)
            {
                return new bool(readFile(is,warn,subrun));
            }
            bool* readFile_GAMBIT(std::istream& is, int subrun)
            {
                return new bool(readFile(is,subrun));
            }
            bool* init_GAMBIT()
            {
                return new bool(init());
            }
            bool* init_GAMBIT(int idAin, int idBin, double eCMin)
            {
                return new bool(init(idAin,idBin,eCMin));
            }
            bool* init_GAMBIT(int idAin, int idBin, double eAin, double eBin)
            {
                return new bool(init(idAin,idBin,eAin,eBin));
            }
            bool* init_GAMBIT(int idAin, int idBin, double pxAin, double pyAin, double pzAin, double pxBin, double pyBin, double pzBin)
            {
                return new bool(init(idAin,idBin,pxAin,pyAin,pzAin,pxBin,pyBin,pzBin));
          } 
#include "Pythia8/abstract_Pythia.h"
namespace Pythia8 { 
  }
         : public virtual Abstract__Pythia    bool* init_GAMBIT(std::string LesHouchesEventFile, bool skipInit)
            {
                return new bool(init(LesHouchesEventFile,skipInit));
            }
            bool* next_GAMBIT()
            {
                return new bool(next());
            }
            int* forceTimeShower_GAMBIT(int iBeg, int iEnd, double pTmax, int nBranchMax)
            {
                return new int(forceTimeShower(iBeg,iEnd,pTmax,nBranchMax));
            }
            bool* forceHadronLevel_GAMBIT(bool findJunctions)
            {
                return new bool(forceHadronLevel(findJunctions));
            }
            bool* moreDecays_GAMBIT()
            {
                return new bool(moreDecays());
            }
            bool* forceRHadronDecays_GAMBIT()
            {
                return new bool(forceRHadronDecays());
            }
            void LHAeventList_GAMBIT(std::ostream& os)
            {
                LHAeventList(os);
            }
            bool* LHAeventSkip_GAMBIT(int nSkip)
            {
                return new bool(LHAeventSkip(nSkip));
            }
            void stat_GAMBIT()
            {
                stat();
            }
            void statistics_GAMBIT(bool all, bool reset)
            {
                statistics(all,reset);
            }
            bool* flag_GAMBIT(std::string key)
            {
                return new bool(flag(key));
            }
            int* mode_GAMBIT(std::string key)
            {
                return new int(mode(key));
            }
            double* parm_GAMBIT(std::string key)
            {
                return new double(parm(key));
            }
            std::string* word_GAMBIT(std::string key)
            {
                return new std::string(word(key));
            }
        private:
            void banner_GAMBIT(std::ostream& os)
            {
                banner(os);
            }
            int* readSubrun_GAMBIT(std::string line, bool warn, std::ostream& os)
            {
                return new int(readSubrun(line,warn,os));
            }
            int* readCommented_GAMBIT(std::string line)
            {
                return new int(readCommented(line));
            }
            void checkSettings_GAMBIT()
            {
                checkSettings();
            }
            bool* checkBeams_GAMBIT()
            {
                return new bool(checkBeams());
            }
            bool* initKinematics_GAMBIT()
            {
                return new bool(initKinematics());
            }
            bool* initPDFs_GAMBIT()
            {
                return new bool(initPDFs());
            }
            void nextKinematics_GAMBIT()
            {
                nextKinematics();
            }
            void boostAndVertex_GAMBIT(bool toLab, bool setVertex)
            {
                boostAndVertex(toLab,setVertex);
            }
            bool* doRHadronDecays_GAMBIT()
            {
                return new bool(doRHadronDecays());
            }
            bool* check_GAMBIT(std::ostream& os)
            {
                return new bool(check(os));
            }
            bool* initSLHA_GAMBIT()
            {
                return new bool(initSLHA());
            }
