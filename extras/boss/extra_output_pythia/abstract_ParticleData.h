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

            virtual bool* init_GAMBIT(std::string startFile) {};
            bool* init(std::string startFile)
            {
                return init_GAMBIT( startFile);
            }

            virtual bool* reInit_GAMBIT(std::string startFile, bool xmlFormat) {};
            bool* reInit(std::string startFile, bool xmlFormat)
            {
                return reInit_GAMBIT( startFile,  xmlFormat);
            }

            virtual bool* readXML_GAMBIT(std::string inFile, bool reset) {};
            bool* readXML(std::string inFile, bool reset)
            {
                return readXML_GAMBIT( inFile,  reset);
            }

            virtual void listXML_GAMBIT(std::string outFile) {};
            void listXML(std::string outFile)
            {
                listXML_GAMBIT( outFile);
            }

            virtual bool* readFF_GAMBIT(std::string inFile, bool reset) {};
            bool* readFF(std::string inFile, bool reset)
            {
                return readFF_GAMBIT( inFile,  reset);
            }

            virtual void listFF_GAMBIT(std::string outFile) {};
            void listFF(std::string outFile)
            {
                listFF_GAMBIT( outFile);
            }

            virtual bool* readString_GAMBIT(std::string lineIn, bool warn, std::ostream& os) {};
            bool* readString(std::string lineIn, bool warn, std::ostream& os)
            {
                return readString_GAMBIT( lineIn,  warn,  os);
            }

            virtual bool* readingFailed_GAMBIT() {};
            bool* readingFailed()
            {
                return readingFailed_GAMBIT();
            }

            virtual void listAll_GAMBIT(std::ostream& os) {};
            void listAll(std::ostream& os)
            {
                listAll_GAMBIT( os);
            }

            virtual void listChanged_GAMBIT(std::ostream& os) {};
            void listChanged(std::ostream& os)
            {
                listChanged_GAMBIT( os);
            }

            virtual void listChanged_GAMBIT(bool changedRes, std::ostream& os) {};
            void listChanged(bool changedRes, std::ostream& os)
            {
                listChanged_GAMBIT( changedRes,  os);
            }

            virtual void list_GAMBIT(bool changedOnly, bool changedRes, std::ostream& os) {};
            void list(bool changedOnly, bool changedRes, std::ostream& os)
            {
                list_GAMBIT( changedOnly,  changedRes,  os);
            }

            virtual void list_GAMBIT(int idList, std::ostream& os) {};
            void list(int idList, std::ostream& os)
            {
                list_GAMBIT( idList,  os);
            }

            virtual void list_GAMBIT(std::vector<int, std::allocator<int> > idList, std::ostream& os) {};
            void list(std::vector<int, std::allocator<int> > idList, std::ostream& os)
            {
                list_GAMBIT( idList,  os);
            }

            virtual void checkTable_GAMBIT(std::ostream& os) {};
            void checkTable(std::ostream& os)
            {
                checkTable_GAMBIT( os);
            }

            virtual void checkTable_GAMBIT(int verbosity, std::ostream& os) {};
            void checkTable(int verbosity, std::ostream& os)
            {
                checkTable_GAMBIT( verbosity,  os);
            }

            virtual void addParticle_GAMBIT(int idIn, std::string nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In) {};
            void addParticle(int idIn, std::string nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In)
            {
                addParticle_GAMBIT( idIn,  nameIn,  spinTypeIn,  chargeTypeIn,  colTypeIn,  m0In,  mWidthIn,  mMinIn,  mMaxIn,  tau0In);
            }

            virtual void addParticle_GAMBIT(int idIn, std::string nameIn, std::string antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In) {};
            void addParticle(int idIn, std::string nameIn, std::string antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In)
            {
                addParticle_GAMBIT( idIn,  nameIn,  antiNameIn,  spinTypeIn,  chargeTypeIn,  colTypeIn,  m0In,  mWidthIn,  mMinIn,  mMaxIn,  tau0In);
            }

            virtual void setAll_GAMBIT(int idIn, std::string nameIn, std::string antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In) {};
            void setAll(int idIn, std::string nameIn, std::string antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In)
            {
                setAll_GAMBIT( idIn,  nameIn,  antiNameIn,  spinTypeIn,  chargeTypeIn,  colTypeIn,  m0In,  mWidthIn,  mMinIn,  mMaxIn,  tau0In);
            }

            virtual bool* isParticle_GAMBIT(int idIn) {};
            bool* isParticle(int idIn)
            {
                return isParticle_GAMBIT( idIn);
            }

            virtual int* nextId_GAMBIT(int idIn) {};
            int* nextId(int idIn)
            {
                return nextId_GAMBIT( idIn);
            }

            virtual void name_GAMBIT(int idIn, std::string nameIn) {};
            void name(int idIn, std::string nameIn)
            {
                name_GAMBIT( idIn,  nameIn);
            }

            virtual void antiName_GAMBIT(int idIn, std::string antiNameIn) {};
            void antiName(int idIn, std::string antiNameIn)
            {
                antiName_GAMBIT( idIn,  antiNameIn);
            }

            virtual void names_GAMBIT(int idIn, std::string nameIn, std::string antiNameIn) {};
            void names(int idIn, std::string nameIn, std::string antiNameIn)
            {
                names_GAMBIT( idIn,  nameIn,  antiNameIn);
            }

            virtual void spinType_GAMBIT(int idIn, int spinTypeIn) {};
            void spinType(int idIn, int spinTypeIn)
            {
                spinType_GAMBIT( idIn,  spinTypeIn);
            }

            virtual void chargeType_GAMBIT(int idIn, int chargeTypeIn) {};
            void chargeType(int idIn, int chargeTypeIn)
            {
                chargeType_GAMBIT( idIn,  chargeTypeIn);
            }

            virtual void colType_GAMBIT(int idIn, int colTypeIn) {};
            void colType(int idIn, int colTypeIn)
            {
                colType_GAMBIT( idIn,  colTypeIn);
            }

            virtual void m0_GAMBIT(int idIn, double m0In) {};
            void m0(int idIn, double m0In)
            {
                m0_GAMBIT( idIn,  m0In);
            }

            virtual void mWidth_GAMBIT(int idIn, double mWidthIn) {};
            void mWidth(int idIn, double mWidthIn)
            {
                mWidth_GAMBIT( idIn,  mWidthIn);
            }

            virtual void mMin_GAMBIT(int idIn, double mMinIn) {};
            void mMin(int idIn, double mMinIn)
            {
                mMin_GAMBIT( idIn,  mMinIn);
            }

            virtual void mMax_GAMBIT(int idIn, double mMaxIn) {};
            void mMax(int idIn, double mMaxIn)
            {
                mMax_GAMBIT( idIn,  mMaxIn);
            }

            virtual void tau0_GAMBIT(int idIn, double tau0In) {};
            void tau0(int idIn, double tau0In)
            {
                tau0_GAMBIT( idIn,  tau0In);
            }

            virtual void isResonance_GAMBIT(int idIn, bool isResonanceIn) {};
            void isResonance(int idIn, bool isResonanceIn)
            {
                isResonance_GAMBIT( idIn,  isResonanceIn);
            }

            virtual void mayDecay_GAMBIT(int idIn, bool mayDecayIn) {};
            void mayDecay(int idIn, bool mayDecayIn)
            {
                mayDecay_GAMBIT( idIn,  mayDecayIn);
            }

            virtual void doExternalDecay_GAMBIT(int idIn, bool doExternalDecayIn) {};
            void doExternalDecay(int idIn, bool doExternalDecayIn)
            {
                doExternalDecay_GAMBIT( idIn,  doExternalDecayIn);
            }

            virtual void isVisible_GAMBIT(int idIn, bool isVisibleIn) {};
            void isVisible(int idIn, bool isVisibleIn)
            {
                isVisible_GAMBIT( idIn,  isVisibleIn);
            }

            virtual void doForceWidth_GAMBIT(int idIn, bool doForceWidthIn) {};
            void doForceWidth(int idIn, bool doForceWidthIn)
            {
                doForceWidth_GAMBIT( idIn,  doForceWidthIn);
            }

            virtual void hasChanged_GAMBIT(int idIn, bool hasChangedIn) {};
            void hasChanged(int idIn, bool hasChangedIn)
            {
                hasChanged_GAMBIT( idIn,  hasChangedIn);
            }

            virtual bool* hasAnti_GAMBIT(int idIn) {};
            bool* hasAnti(int idIn)
            {
                return hasAnti_GAMBIT( idIn);
            }

            virtual std::string* name_GAMBIT(int idIn) {};
            std::string* name(int idIn)
            {
                return name_GAMBIT( idIn);
            }

            virtual int* spinType_GAMBIT(int idIn) {};
            int* spinType(int idIn)
            {
                return spinType_GAMBIT( idIn);
            }

            virtual int* chargeType_GAMBIT(int idIn) {};
            int* chargeType(int idIn)
            {
                return chargeType_GAMBIT( idIn);
            }

            virtual double* charge_GAMBIT(int idIn) {};
            double* charge(int idIn)
            {
                return charge_GAMBIT( idIn);
            }

            virtual int* colType_GAMBIT(int idIn) {};
            int* colType(int idIn)
            {
                return colType_GAMBIT( idIn);
            }

            virtual double* m0_GAMBIT(int idIn) {};
            double* m0(int idIn)
            {
                return m0_GAMBIT( idIn);
            }

            virtual double* mWidth_GAMBIT(int idIn) {};
            double* mWidth(int idIn)
            {
                return mWidth_GAMBIT( idIn);
            }

            virtual double* mMin_GAMBIT(int idIn) {};
            double* mMin(int idIn)
            {
                return mMin_GAMBIT( idIn);
            }

            virtual double* m0Min_GAMBIT(int idIn) {};
            double* m0Min(int idIn)
            {
                return m0Min_GAMBIT( idIn);
            }

            virtual double* mMax_GAMBIT(int idIn) {};
            double* mMax(int idIn)
            {
                return mMax_GAMBIT( idIn);
            }

            virtual double* m0Max_GAMBIT(int idIn) {};
            double* m0Max(int idIn)
            {
                return m0Max_GAMBIT( idIn);
            }

            virtual double* tau0_GAMBIT(int idIn) {};
            double* tau0(int idIn)
            {
                return tau0_GAMBIT( idIn);
            }

            virtual bool* isResonance_GAMBIT(int idIn) {};
            bool* isResonance(int idIn)
            {
                return isResonance_GAMBIT( idIn);
            }

            virtual bool* mayDecay_GAMBIT(int idIn) {};
            bool* mayDecay(int idIn)
            {
                return mayDecay_GAMBIT( idIn);
            }

            virtual bool* doExternalDecay_GAMBIT(int idIn) {};
            bool* doExternalDecay(int idIn)
            {
                return doExternalDecay_GAMBIT( idIn);
            }

            virtual bool* isVisible_GAMBIT(int idIn) {};
            bool* isVisible(int idIn)
            {
                return isVisible_GAMBIT( idIn);
            }

            virtual bool* doForceWidth_GAMBIT(int idIn) {};
            bool* doForceWidth(int idIn)
            {
                return doForceWidth_GAMBIT( idIn);
            }

            virtual bool* hasChanged_GAMBIT(int idIn) {};
            bool* hasChanged(int idIn)
            {
                return hasChanged_GAMBIT( idIn);
            }

            virtual bool* useBreitWigner_GAMBIT(int idIn) {};
            bool* useBreitWigner(int idIn)
            {
                return useBreitWigner_GAMBIT( idIn);
            }

            virtual double* constituentMass_GAMBIT(int idIn) {};
            double* constituentMass(int idIn)
            {
                return constituentMass_GAMBIT( idIn);
            }

            virtual double* mSel_GAMBIT(int idIn) {};
            double* mSel(int idIn)
            {
                return mSel_GAMBIT( idIn);
            }

            virtual double* mRun_GAMBIT(int idIn, double mH) {};
            double* mRun(int idIn, double mH)
            {
                return mRun_GAMBIT( idIn,  mH);
            }

            virtual bool* canDecay_GAMBIT(int idIn) {};
            bool* canDecay(int idIn)
            {
                return canDecay_GAMBIT( idIn);
            }

            virtual bool* isLepton_GAMBIT(int idIn) {};
            bool* isLepton(int idIn)
            {
                return isLepton_GAMBIT( idIn);
            }

            virtual bool* isQuark_GAMBIT(int idIn) {};
            bool* isQuark(int idIn)
            {
                return isQuark_GAMBIT( idIn);
            }

            virtual bool* isGluon_GAMBIT(int idIn) {};
            bool* isGluon(int idIn)
            {
                return isGluon_GAMBIT( idIn);
            }

            virtual bool* isDiquark_GAMBIT(int idIn) {};
            bool* isDiquark(int idIn)
            {
                return isDiquark_GAMBIT( idIn);
            }

            virtual bool* isParton_GAMBIT(int idIn) {};
            bool* isParton(int idIn)
            {
                return isParton_GAMBIT( idIn);
            }

            virtual bool* isHadron_GAMBIT(int idIn) {};
            bool* isHadron(int idIn)
            {
                return isHadron_GAMBIT( idIn);
            }

            virtual bool* isMeson_GAMBIT(int idIn) {};
            bool* isMeson(int idIn)
            {
                return isMeson_GAMBIT( idIn);
            }

            virtual bool* isBaryon_GAMBIT(int idIn) {};
            bool* isBaryon(int idIn)
            {
                return isBaryon_GAMBIT( idIn);
            }

            virtual bool* isOctetHadron_GAMBIT(int idIn) {};
            bool* isOctetHadron(int idIn)
            {
                return isOctetHadron_GAMBIT( idIn);
            }

            virtual int* heaviestQuark_GAMBIT(int idIn) {};
            int* heaviestQuark(int idIn)
            {
                return heaviestQuark_GAMBIT( idIn);
            }

            virtual int* baryonNumberType_GAMBIT(int idIn) {};
            int* baryonNumberType(int idIn)
            {
                return baryonNumberType_GAMBIT( idIn);
            }

            virtual void rescaleBR_GAMBIT(int idIn, double newSumBR) {};
            void rescaleBR(int idIn, double newSumBR)
            {
                rescaleBR_GAMBIT( idIn,  newSumBR);
            }

            virtual void resInit_GAMBIT(int idIn) {};
            void resInit(int idIn)
            {
                resInit_GAMBIT( idIn);
            }

            virtual double* resWidth_GAMBIT(int idIn, double mHat, int idInFlav, bool openOnly, bool setBR) {};
            double* resWidth(int idIn, double mHat, int idInFlav, bool openOnly, bool setBR)
            {
                return resWidth_GAMBIT( idIn,  mHat,  idInFlav,  openOnly,  setBR);
            }

            virtual double* resWidthOpen_GAMBIT(int idIn, double mHat, int idInFlav) {};
            double* resWidthOpen(int idIn, double mHat, int idInFlav)
            {
                return resWidthOpen_GAMBIT( idIn,  mHat,  idInFlav);
            }

            virtual double* resWidthStore_GAMBIT(int idIn, double mHat, int idInFlav) {};
            double* resWidthStore(int idIn, double mHat, int idInFlav)
            {
                return resWidthStore_GAMBIT( idIn,  mHat,  idInFlav);
            }

            virtual double* resOpenFrac_GAMBIT(int id1In, int id2In, int id3In) {};
            double* resOpenFrac(int id1In, int id2In, int id3In)
            {
                return resOpenFrac_GAMBIT( id1In,  id2In,  id3In);
            }

            virtual double* resWidthRescaleFactor_GAMBIT(int idIn) {};
            double* resWidthRescaleFactor(int idIn)
            {
                return resWidthRescaleFactor_GAMBIT( idIn);
            }

            virtual double* resWidthChan_GAMBIT(int idIn, double mHat, int idAbs1, int idAbs2) {};
            double* resWidthChan(int idIn, double mHat, int idAbs1, int idAbs2)
            {
                return resWidthChan_GAMBIT( idIn,  mHat,  idAbs1,  idAbs2);
            }
        private:

            virtual void initCommon_GAMBIT() {};
            void initCommon()
            {
                initCommon_GAMBIT();
            }

            virtual std::string* toLower_GAMBIT(const std::string& nameConv) {};
            std::string* toLower(const std::string& nameConv)
            {
                return toLower_GAMBIT( nameConv);
            }

            virtual bool* boolString_GAMBIT(std::string tag) {};
            bool* boolString(std::string tag)
            {
                return boolString_GAMBIT( tag);
            }

            virtual std::string* attributeValue_GAMBIT(std::string line, std::string attribute) {};
            std::string* attributeValue(std::string line, std::string attribute)
            {
                return attributeValue_GAMBIT( line,  attribute);
            }

            virtual bool* boolAttributeValue_GAMBIT(std::string line, std::string attribute) {};
            bool* boolAttributeValue(std::string line, std::string attribute)
            {
                return boolAttributeValue_GAMBIT( line,  attribute);
            }

            virtual int* intAttributeValue_GAMBIT(std::string line, std::string attribute) {};
            int* intAttributeValue(std::string line, std::string attribute)
            {
                return intAttributeValue_GAMBIT( line,  attribute);
            }

            virtual double* doubleAttributeValue_GAMBIT(std::string line, std::string attribute) {};
            double* doubleAttributeValue(std::string line, std::string attribute)
            {
                return doubleAttributeValue_GAMBIT( line,  attribute);
            }

        public:
            ParticleData* downcast()
            {
                return reinterpret_cast<ParticleData*>(this);
            }
    };
}
#pragma GCC diagnostic pop

