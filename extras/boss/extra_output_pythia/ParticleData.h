
        public:
            bool* init_GAMBIT(std::string startFile)
            {
                return new bool(init(startFile));
            }
            bool* reInit_GAMBIT(std::string startFile, bool xmlFormat)
            {
                return new bool(reInit(startFile,xmlFormat));
            }
            bool* readXML_GAMBIT(std::string inFile, bool reset)
            {
                return new bool(readXML(inFile,reset));
            }
            void listXML_GAMBIT(std::string outFile)
            {
                listXML(outFile);
            }
            bool* readFF_GAMBIT(std::string inFile, bool reset)
            {
                return new bool(readFF(inFile,reset));
            }
            void listFF_GAMBIT(std::string outFile)
            {
                listFF(outFile);
            }
            bool* readString_GAMBIT(std::string lineIn, bool warn, std::ostream& os)
            {
                return new bool(readString(lineIn,warn,os));
            }
            bool* readingFailed_GAMBIT()
            {
                return new bool(readingFailed());
            }
            void listAll_GAMBIT(std::ostream& os)
            {
                listAll(os);
            }
            void listChanged_GAMBIT(std::ostream& os)
            {
                listChanged(os);
            }
            void listChanged_GAMBIT(bool changedRes, std::ostream& os)
            {
                listChanged(changedRes,os);
            }
            void list_GAMBIT(bool changedOnly, bool changedRes, std::ostream& os)
            {
                list(changedOnly,changedRes,os);
            }
            void list_GAMBIT(int idList, std::ostream& os)
            {
                list(idList,os);
            }
            void list_GAMBIT(std::vector<int, std::allocator<int> > idList, std::ostream& os)
            {
                list(idList,os);
            }
            void checkTable_GAMBIT(std::ostream& os)
            {
                checkTable(os);
            }
            void checkTable_GAMBIT(int verbosity, std::ostream& os)
            {
                checkTable(verbosity,os);
            }
            void addParticle_GAMBIT(int idIn, std::string nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In)
            {
                addParticle(idIn,nameIn,spinTypeIn,chargeTypeIn,colTypeIn,m0In,mWidthIn,mMinIn,mMaxIn,tau0In);
            }
            void addParticle_GAMBIT(int idIn, std::string nameIn, std::string antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In)
            {
                addParticle(idIn,nameIn,antiNameIn,spinTypeIn,chargeTypeIn,colTypeIn,m0In,mWidthIn,mMinIn,mMaxIn,tau0In);
            }
            void setAll_GAMBIT(int idIn, std::string nameIn, std::string antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In)
            {
                setAll(idIn,nameIn,antiNameIn,spinTypeIn,chargeTypeIn,colTypeIn,m0In,mWidthIn,mMinIn,mMaxIn,tau0In);
            }
            bool* isParticle_GAMBIT(int idIn)
            {
                return new bool(isParticle(idIn));
            }
            int* nextId_GAMBIT(int idIn)
            {
                return new int(nextId(idIn));
            }
            void name_GAMBIT(int idIn, std::string nameIn)
            {
                name(idIn,nameIn);
            }
            void antiName_GAMBIT(int idIn, std::string antiNameIn)
            {
                antiName(idIn,antiNameIn);
            }
            void names_GAMBIT(int idIn, std::string nameIn, std::string antiNameIn)
            {
                names(idIn,nameIn,antiNameIn);
            }
            void spinType_GAMBIT(int idIn, int spinTypeIn)
            {
                spinType(idIn,spinTypeIn);
            }
            void chargeType_GAMBIT(int idIn, int chargeTypeIn)
            {
                chargeType(idIn,chargeTypeIn);
            }
            void colType_GAMBIT(int idIn, int colTypeIn)
            {
                colType(idIn,colTypeIn);
            }
            void m0_GAMBIT(int idIn, double m0In)
            {
                m0(idIn,m0In);
            }
            void mWidth_GAMBIT(int idIn, double mWidthIn)
            {
                mWidth(idIn,mWidthIn);
            }
            void mMin_GAMBIT(int idIn, double mMinIn)
            {
                mMin(idIn,mMinIn);
            }
            void mMax_GAMBIT(int idIn, double mMaxIn)
            {
                mMax(idIn,mMaxIn);
            }
            void tau0_GAMBIT(int idIn, double tau0In)
            {
                tau0(idIn,tau0In);
            }
            void isResonance_GAMBIT(int idIn, bool isResonanceIn)
            {
                isResonance(idIn,isResonanceIn);
            }
            void mayDecay_GAMBIT(int idIn, bool mayDecayIn)
            {
                mayDecay(idIn,mayDecayIn);
            }
            void doExternalDecay_GAMBIT(int idIn, bool doExternalDecayIn)
            {
                doExternalDecay(idIn,doExternalDecayIn);
            }
            void isVisible_GAMBIT(int idIn, bool isVisibleIn)
            {
                isVisible(idIn,isVisibleIn);
            }
            void doForceWidth_GAMBIT(int idIn, bool doForceWidthIn)
            {
                doForceWidth(idIn,doForceWidthIn);
            }
            void hasChanged_GAMBIT(int idIn, bool hasChangedIn)
            {
                hasChanged(idIn,hasChangedIn);
            }
            bool* hasAnti_GAMBIT(int idIn)
            {
                return new bool(hasAnti(idIn));
            }
            std::string* name_GAMBIT(int idIn)
            {
                return new std::string(name(idIn));
            }
            int* spinType_GAMBIT(int idIn)
            {
                return new int(spinType(idIn));
            }
            int* chargeType_GAMBIT(int idIn)
            {
                return new int(chargeType(idIn));
            }
            double* charge_GAMBIT(int idIn)
            {
                return new double(charge(idIn));
            }
            int* colType_GAMBIT(int idIn)
            {
                return new int(colType(idIn));
            }
            double* m0_GAMBIT(int idIn)
            {
                return new double(m0(idIn));
            }
            double* mWidth_GAMBIT(int idIn)
            {
                return new double(mWidth(idIn));
            }
            double* mMin_GAMBIT(int idIn)
            {
                return new double(mMin(idIn));
            }
            double* m0Min_GAMBIT(int idIn)
            {
                return new double(m0Min(idIn));
            }
            double* mMax_GAMBIT(int idIn)
            {
                return new double(mMax(idIn));
            }
            double* m0Max_GAMBIT(int idIn)
            {
                return new double(m0Max(idIn));
            }
            double* tau0_GAMBIT(int idIn)
            {
                return new double(tau0(idIn));
            }
            bool* isResonance_GAMBIT(int idIn)
            {
                return new bool(isResonance(idIn));
            }
            bool* mayDecay_GAMBIT(int idIn)
            {
                return new bool(mayDecay(idIn));
            }
            bool* doExternalDecay_GAMBIT(int idIn)
            {
                return new bool(doExternalDecay(idIn));
            }
            bool* isVisible_GAMBIT(int idIn)
            {
                return new bool(isVisible(idIn));
            }
            bool* doForceWidth_GAMBIT(int idIn)
            {
                return new bool(doForceWidth(idIn));
            }
            bool* hasChanged_GAMBIT(int idIn)
            {
                return new bool(hasChanged(idIn));
            }
            bool* useBreitWigner_GAMBIT(int idIn)
            {
                return new bool(useBreitWigner(idIn));
            }
            double* constituentMass_GAMBIT(int idIn)
            {
                return new double(constituentMass(idIn));
            }
            double* mSel_GAMBIT(int idIn)
            {
                return new double(mSel(idIn));
            }
            double* mRun_GAMBIT(int idIn, double mH)
            {
                return new double(mRun(idIn,mH));
            }
            bool* canDecay_GAMBIT(int idIn)
            {
                return new bool(canDecay(idIn));
            }
            bool* isLepton_GAMBIT(int idIn)
            {
                return new bool(isLepton(idIn));
            }
            bool* isQuark_GAMBIT(int idIn)
            {
                return new bool(isQuark(idIn));
            }
            bool* isGluon_GAMBIT(int idIn)
            {
                return new bool(isGluon(idIn));
            }
            bool* isDiquark_GAMBIT(int idIn)
            {
                return new bool(isDiquark(idIn));
            }
            bool* isParton_GAMBIT(int idIn)
            {
                return new bool(isParton(idIn));
            }
            bool* isHadron_GAMBIT(int idIn)
            {
                return new bool(isHadron(idIn));
            }
            bool* isMeson_GAMBIT(int idIn)
            {
                return new bool(isMeson(idIn));
            }
            bool* isBaryon_GAMBIT(int idIn)
            {
                return new bool(isBaryon(idIn));
            }
            bool* isOctetHadron_GAMBIT(int idIn)
            {
                return new bool(isOctetHadron(idIn));
            }
            int* heaviestQuark_GAMBIT(int idIn)
            {
                return new int(heaviestQuark(idIn));
            }
            int* baryonNumberType_GAMBIT(int idIn)
            {
                return new int(baryonNumberType(idIn));
            }
            void rescaleBR_GAMBIT(int idIn, double newSumBR)
            {
                rescaleBR(idIn,newSumBR);
            }
            void resInit_GAMBIT(int idIn)
            {
                resInit(idIn);
            }
            double* resWidth_GAMBIT(int idIn, double mHat, int idInFlav, bool openOnly, bool setBR)
            {
                return new double(resWidth(idIn,mHat,idInFlav,openOnly,setBR));
            }
            double* resWidthOpen_GAMBIT(int idIn, double mHat, int idInFlav)
            {
                return new double(resWidthOpen(idIn,mHat,idInFlav));
            }
            double* resWidthStore_GAMBIT(int idIn, double mHat, int idInFlav)
            {
                return new double(resWidthStore(idIn,mHat,idInFlav));
            }
            double* resOpenFrac_GAMBIT(int id1In, int id2In, int id3In)
            {
                return new double(resOpenFrac(id1In,id2In,id3In));
            }
            double* resWidthRescaleFactor_GAMBIT(int idIn)
            {
                return new double(resWidthRescaleFactor(idIn));
            }
            double* resWidthChan_GAMBIT(int idIn, double mHat, int idAbs1, int idAbs2)
            {
                return new double(resWidthChan(idIn,mHat,idAbs1,idAbs2));
            }
        private:
            void initCommon_GAMBIT()
            {
                initCommon();
            }
            std::string* toLower_GAMBIT(const std::string& nameConv)
            {
                return new std::string(toLower(nameConv));
            }
            bool* boolString_GAMBIT(std::string tag)
            {
                return new bool(boolString(tag));
            }
            std::string* attributeValue_GAMBIT(std::string line, std::string attribute)
            {
                return new std::string(attributeValue(line,attribute));
            }
            bool* boolAttributeValue_GAMBIT(std::string line, std::string attribute)
            {
                return new bool(boolAttributeValue(line,attribute));
            }
            int* intAttributeValue_GAMBIT(std::string line, std::string attribute)
            {
                return new int(intAttributeValue(line,attribute));
            }
            double* doubleAttributeValue_GAMBIT(std::string line, std::string attribute)
            {
                return new double(doubleAttributeValue(line,attribute));
            }
 : public virtual Abstract__ParticleData} 
#include "Pythia8/abstract_ParticleData.h"
namespace Pythia8 { 
