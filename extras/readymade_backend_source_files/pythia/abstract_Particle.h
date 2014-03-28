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
    class Abstract__Particle
    {
        private:
        public:
            // UNKNOWN: OperatorMethod

            virtual void setEvtPtr_GAMBIT(Pythia8::Abstract__Event* evtPtrIn) {};
            void setEvtPtr(Pythia8::Abstract__Event* evtPtrIn)
            {
                setEvtPtr_GAMBIT( evtPtrIn);
            }

            virtual void id(int idIn) {};

            virtual void status(int statusIn) {};

            virtual void statusPos() {};

            virtual void statusNeg() {};

            virtual void statusCode(int statusIn) {};

            virtual void mother1(int mother1In) {};

            virtual void mother2(int mother2In) {};

            virtual void mothers(int mother1In, int mother2In) {};

            virtual void daughter1(int daughter1In) {};

            virtual void daughter2(int daughter2In) {};

            virtual void daughters(int daughter1In, int daughter2In) {};

            virtual void col(int colIn) {};

            virtual void acol(int acolIn) {};

            virtual void cols(int colIn, int acolIn) {};

            virtual void p_GAMBIT(Pythia8::Abstract__Vec4 pIn) {};
            void p(Pythia8::Abstract__Vec4 pIn)
            {
                p_GAMBIT( pIn);
            }

            virtual void p(double pxIn, double pyIn, double pzIn, double eIn) {};

            virtual void px(double pxIn) {};

            virtual void py(double pyIn) {};

            virtual void pz(double pzIn) {};

            virtual void e(double eIn) {};

            virtual void m(double mIn) {};

            virtual void scale(double scaleIn) {};

            virtual void pol(double polIn) {};

            virtual void vProd_GAMBIT(Pythia8::Abstract__Vec4 vProdIn) {};
            void vProd(Pythia8::Abstract__Vec4 vProdIn)
            {
                vProd_GAMBIT( vProdIn);
            }

            virtual void vProd(double xProdIn, double yProdIn, double zProdIn, double tProdIn) {};

            virtual void xProd(double xProdIn) {};

            virtual void yProd(double yProdIn) {};

            virtual void zProd(double zProdIn) {};

            virtual void tProd(double tProdIn) {};

            virtual void tau(double tauIn) {};

            virtual int id() {};

            virtual int status() {};

            virtual int mother1() {};

            virtual int mother2() {};

            virtual int daughter1() {};

            virtual int daughter2() {};

            virtual int col() {};

            virtual int acol() {};

            virtual Pythia8::Vec4* p_GAMBIT() {};
            Pythia8::Abstract__Vec4* p()
            {
                return reinterpret_cast<Pythia8::Abstract__Vec4*>(p_GAMBIT());
            }

            virtual double px() {};

            virtual double py() {};

            virtual double pz() {};

            virtual double e() {};

            virtual double m() {};

            virtual double scale() {};

            virtual double pol() {};

            virtual bool hasVertex() {};

            virtual Pythia8::Vec4* vProd_GAMBIT() {};
            Pythia8::Abstract__Vec4* vProd()
            {
                return reinterpret_cast<Pythia8::Abstract__Vec4*>(vProd_GAMBIT());
            }

            virtual double xProd() {};

            virtual double yProd() {};

            virtual double zProd() {};

            virtual double tProd() {};

            virtual double tau() {};

            virtual int idAbs() {};

            virtual int statusAbs() {};

            virtual bool isFinal() {};

            virtual bool isRescatteredIncoming() {};

            virtual double m2() {};

            virtual double mCalc() {};

            virtual double m2Calc() {};

            virtual double eCalc() {};

            virtual double pT() {};

            virtual double pT2() {};

            virtual double mT() {};

            virtual double mT2() {};

            virtual double pAbs() {};

            virtual double pAbs2() {};

            virtual double eT() {};

            virtual double eT2() {};

            virtual double theta() {};

            virtual double phi() {};

            virtual double thetaXZ() {};

            virtual double pPos() {};

            virtual double pNeg() {};

            virtual double y() {};

            virtual double eta() {};

            virtual Pythia8::Vec4* vDec_GAMBIT() {};
            Pythia8::Abstract__Vec4* vDec()
            {
                return reinterpret_cast<Pythia8::Abstract__Vec4*>(vDec_GAMBIT());
            }

            virtual double xDec() {};

            virtual double yDec() {};

            virtual double zDec() {};

            virtual double tDec() {};

            virtual int index() {};

            virtual int statusHepMC() {};

            virtual int iTopCopy() {};

            virtual int iBotCopy() {};

            virtual int iTopCopyId() {};

            virtual int iBotCopyId() {};

            virtual bool isAncestor(int iAncestor) {};

            virtual bool undoDecay() {};

            virtual std::string name() {};

            virtual std::string nameWithStatus(int maxLen) {};

            virtual int spinType() {};

            virtual int chargeType() {};

            virtual double charge() {};

            virtual bool isCharged() {};

            virtual bool isNeutral() {};

            virtual int colType() {};

            virtual double m0() {};

            virtual double mWidth() {};

            virtual double mMin() {};

            virtual double mMax() {};

            virtual double mSel() {};

            virtual double constituentMass() {};

            virtual double tau0() {};

            virtual bool mayDecay() {};

            virtual bool canDecay() {};

            virtual bool doExternalDecay() {};

            virtual bool isResonance() {};

            virtual bool isVisible() {};

            virtual bool isLepton() {};

            virtual bool isQuark() {};

            virtual bool isGluon() {};

            virtual bool isDiquark() {};

            virtual bool isParton() {};

            virtual bool isHadron() {};

            virtual void rescale3(double fac) {};

            virtual void rescale4(double fac) {};

            virtual void rescale5(double fac) {};

            virtual void rot(double thetaIn, double phiIn) {};

            virtual void bst(double betaX, double betaY, double betaZ) {};

            virtual void bst(double betaX, double betaY, double betaZ, double gamma) {};

            virtual void bst_GAMBIT(const Pythia8::Abstract__Vec4& pBst) {};
            void bst(const Pythia8::Abstract__Vec4& pBst)
            {
                bst_GAMBIT( pBst);
            }

            virtual void bst_GAMBIT(const Pythia8::Abstract__Vec4& pBst, double mBst) {};
            void bst(const Pythia8::Abstract__Vec4& pBst, double mBst)
            {
                bst_GAMBIT( pBst,  mBst);
            }

            virtual void bstback_GAMBIT(const Pythia8::Abstract__Vec4& pBst) {};
            void bstback(const Pythia8::Abstract__Vec4& pBst)
            {
                bstback_GAMBIT( pBst);
            }

            virtual void bstback_GAMBIT(const Pythia8::Abstract__Vec4& pBst, double mBst) {};
            void bstback(const Pythia8::Abstract__Vec4& pBst, double mBst)
            {
                bstback_GAMBIT( pBst,  mBst);
            }

            virtual void rotbst_GAMBIT(const Pythia8::Abstract__RotBstMatrix& M) {};
            void rotbst(const Pythia8::Abstract__RotBstMatrix& M)
            {
                rotbst_GAMBIT( M);
            }

            virtual void offsetHistory(int minMother, int addMother, int minDaughter, int addDaughter) {};

            virtual void offsetCol(int addCol) {};

        public:
            Particle* downcast()
            {
                return reinterpret_cast<Particle*>(this);
            }
    };
}
#pragma GCC diagnostic pop

