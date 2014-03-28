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
    class Abstract__ParticleData
    {
        private:
        public:

            virtual bool init(std::string startFile) {};

            virtual bool reInit(std::string startFile, bool xmlFormat) {};

            virtual bool readXML(std::string inFile, bool reset) {};

            virtual void listXML(std::string outFile) {};

            virtual bool readFF(std::string inFile, bool reset) {};

            virtual void listFF(std::string outFile) {};

            virtual bool readString(std::string lineIn, bool warn, std::ostream& os) {};

            virtual bool readingFailed() {};

            virtual void listAll(std::ostream& os) {};

            virtual void listChanged(std::ostream& os) {};

            virtual void listChanged(bool changedRes, std::ostream& os) {};

            virtual void list(bool changedOnly, bool changedRes, std::ostream& os) {};

            virtual void list(int idList, std::ostream& os) {};

            virtual void list(std::vector<int, std::allocator<int> > idList, std::ostream& os) {};

            virtual void checkTable(std::ostream& os) {};

            virtual void checkTable(int verbosity, std::ostream& os) {};

            virtual void addParticle(int idIn, std::string nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In) {};

            virtual void addParticle(int idIn, std::string nameIn, std::string antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In) {};

            virtual void setAll(int idIn, std::string nameIn, std::string antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In) {};

            virtual bool isParticle(int idIn) {};

            virtual int nextId(int idIn) {};

            virtual void name(int idIn, std::string nameIn) {};

            virtual void antiName(int idIn, std::string antiNameIn) {};

            virtual void names(int idIn, std::string nameIn, std::string antiNameIn) {};

            virtual void spinType(int idIn, int spinTypeIn) {};

            virtual void chargeType(int idIn, int chargeTypeIn) {};

            virtual void colType(int idIn, int colTypeIn) {};

            virtual void m0(int idIn, double m0In) {};

            virtual void mWidth(int idIn, double mWidthIn) {};

            virtual void mMin(int idIn, double mMinIn) {};

            virtual void mMax(int idIn, double mMaxIn) {};

            virtual void tau0(int idIn, double tau0In) {};

            virtual void isResonance(int idIn, bool isResonanceIn) {};

            virtual void mayDecay(int idIn, bool mayDecayIn) {};

            virtual void doExternalDecay(int idIn, bool doExternalDecayIn) {};

            virtual void isVisible(int idIn, bool isVisibleIn) {};

            virtual void doForceWidth(int idIn, bool doForceWidthIn) {};

            virtual void hasChanged(int idIn, bool hasChangedIn) {};

            virtual bool hasAnti(int idIn) {};

            virtual std::string name(int idIn) {};

            virtual int spinType(int idIn) {};

            virtual int chargeType(int idIn) {};

            virtual double charge(int idIn) {};

            virtual int colType(int idIn) {};

            virtual double m0(int idIn) {};

            virtual double mWidth(int idIn) {};

            virtual double mMin(int idIn) {};

            virtual double m0Min(int idIn) {};

            virtual double mMax(int idIn) {};

            virtual double m0Max(int idIn) {};

            virtual double tau0(int idIn) {};

            virtual bool isResonance(int idIn) {};

            virtual bool mayDecay(int idIn) {};

            virtual bool doExternalDecay(int idIn) {};

            virtual bool isVisible(int idIn) {};

            virtual bool doForceWidth(int idIn) {};

            virtual bool hasChanged(int idIn) {};

            virtual bool useBreitWigner(int idIn) {};

            virtual double constituentMass(int idIn) {};

            virtual double mSel(int idIn) {};

            virtual double mRun(int idIn, double mH) {};

            virtual bool canDecay(int idIn) {};

            virtual bool isLepton(int idIn) {};

            virtual bool isQuark(int idIn) {};

            virtual bool isGluon(int idIn) {};

            virtual bool isDiquark(int idIn) {};

            virtual bool isParton(int idIn) {};

            virtual bool isHadron(int idIn) {};

            virtual bool isMeson(int idIn) {};

            virtual bool isBaryon(int idIn) {};

            virtual bool isOctetHadron(int idIn) {};

            virtual int heaviestQuark(int idIn) {};

            virtual int baryonNumberType(int idIn) {};

            virtual void rescaleBR(int idIn, double newSumBR) {};

            virtual void resInit(int idIn) {};

            virtual double resWidth(int idIn, double mHat, int idInFlav, bool openOnly, bool setBR) {};

            virtual double resWidthOpen(int idIn, double mHat, int idInFlav) {};

            virtual double resWidthStore(int idIn, double mHat, int idInFlav) {};

            virtual double resOpenFrac(int id1In, int id2In, int id3In) {};

            virtual double resWidthRescaleFactor(int idIn) {};

            virtual double resWidthChan(int idIn, double mHat, int idAbs1, int idAbs2) {};
        private:

            virtual void initCommon() {};

            virtual std::string toLower(const std::string& nameConv) {};

            virtual bool boolString(std::string tag) {};

            virtual std::string attributeValue(std::string line, std::string attribute) {};

            virtual bool boolAttributeValue(std::string line, std::string attribute) {};

            virtual int intAttributeValue(std::string line, std::string attribute) {};

            virtual double doubleAttributeValue(std::string line, std::string attribute) {};

        public:
            ParticleData* downcast()
            {
                return reinterpret_cast<ParticleData*>(this);
            }
    };
}
#pragma GCC diagnostic pop

