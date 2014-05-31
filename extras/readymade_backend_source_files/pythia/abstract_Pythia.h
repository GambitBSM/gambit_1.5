// Forward declarations:
namespace Pythia8
{
    class Abstract__Particle;
    class Particle;
    class Abstract__ParticleData;
    class ParticleData;
    class Abstract__Vec4;
    class Vec4;
    class Abstract__RotBstMatrix;
    class RotBstMatrix;
    class Abstract__Event;
    class Event;
    class Abstract__Pythia;
    class Pythia;
}


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wreturn-type"
namespace Pythia8
{
    class Abstract__Pythia
    {
        public:
        private:
        public:
            virtual ~Abstract__Pythia() {};

            virtual bool readString(std::string line, bool warn) {};

            virtual bool readFile(std::string fileName, bool warn, int subrun) {};

            virtual bool readFile(std::string fileName, int subrun) {};

            virtual bool readFile(std::istream& is, bool warn, int subrun) {};

            virtual bool readFile(std::istream& is, int subrun) {};

            virtual bool init() {};

            virtual bool init(int idAin, int idBin, double eCMin) {};

            virtual bool init(int idAin, int idBin, double eAin, double eBin) {};

            virtual bool init(int idAin, int idBin, double pxAin, double pyAin, double pzAin, double pxBin, double pyBin, double pzBin) {};

            virtual bool init(std::string LesHouchesEventFile, bool skipInit) {};

            virtual bool next() {};

            virtual int forceTimeShower(int iBeg, int iEnd, double pTmax, int nBranchMax) {};

            virtual bool forceHadronLevel(bool findJunctions) {};

            virtual bool moreDecays() {};

            virtual bool forceRHadronDecays() {};

            virtual void LHAeventList(std::ostream& os) {};

            virtual bool LHAeventSkip(int nSkip) {};

            virtual void stat() {};

            virtual void statistics(bool all, bool reset) {};

            virtual bool flag(std::string key) {};

            virtual int mode(std::string key) {};

            virtual double parm(std::string key) {};

            virtual std::string word(std::string key) {};
        private:
            // UNKNOWN: OperatorMethod

            virtual void banner(std::ostream& os) {};

            virtual int readSubrun(std::string line, bool warn, std::ostream& os) {};

            virtual int readCommented(std::string line) {};

            virtual void checkSettings() {};

            virtual bool checkBeams() {};

            virtual bool initKinematics() {};

            virtual bool initPDFs() {};

            virtual void nextKinematics() {};

            virtual void boostAndVertex(bool toLab, bool setVertex) {};

            virtual bool doRHadronDecays() {};

            virtual bool check(std::ostream& os) {};

            virtual bool initSLHA() {};

        public:
            Pythia* downcast()
            {
                return reinterpret_cast<Pythia*>(this);
            }
    };
}
#pragma GCC diagnostic pop

