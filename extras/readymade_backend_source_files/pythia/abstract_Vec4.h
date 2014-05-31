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
    class Abstract__Vec4
    {
        private:
        public:
            // UNKNOWN: OperatorMethod
            // UNKNOWN: OperatorMethod

            virtual void reset() {};

            virtual void p(double xIn, double yIn, double zIn, double tIn) {};

            virtual void p_GAMBIT(Pythia8::Abstract__Vec4 pIn) {};
            void p(Pythia8::Abstract__Vec4 pIn)
            {
                p_GAMBIT( pIn);
            }

            virtual void px(double xIn) {};

            virtual void py(double yIn) {};

            virtual void pz(double zIn) {};

            virtual void e(double tIn) {};

            virtual double px() {};

            virtual double py() {};

            virtual double pz() {};

            virtual double e() {};

            virtual double mCalc() {};

            virtual double m2Calc() {};

            virtual double pT() {};

            virtual double pT2() {};

            virtual double pAbs() {};

            virtual double pAbs2() {};

            virtual double eT() {};

            virtual double eT2() {};

            virtual double theta() {};

            virtual double phi() {};

            virtual double thetaXZ() {};

            virtual double pPos() {};

            virtual double pNeg() {};

            virtual double rap() {};

            virtual double eta() {};

            virtual void rescale3(double fac) {};

            virtual void rescale4(double fac) {};

            virtual void flip3() {};

            virtual void flip4() {};

            virtual void rot(double thetaIn, double phiIn) {};

            virtual void rotaxis(double phiIn, double nx, double ny, double nz) {};

            virtual void rotaxis_GAMBIT(double phiIn, const Pythia8::Abstract__Vec4& n) {};
            void rotaxis(double phiIn, const Pythia8::Abstract__Vec4& n)
            {
                rotaxis_GAMBIT( phiIn,  n);
            }

            virtual void bst(double betaX, double betaY, double betaZ) {};

            virtual void bst(double betaX, double betaY, double betaZ, double gamma) {};

            virtual void bst_GAMBIT(const Pythia8::Abstract__Vec4& pIn) {};
            void bst(const Pythia8::Abstract__Vec4& pIn)
            {
                bst_GAMBIT( pIn);
            }

            virtual void bst_GAMBIT(const Pythia8::Abstract__Vec4& pIn, double mIn) {};
            void bst(const Pythia8::Abstract__Vec4& pIn, double mIn)
            {
                bst_GAMBIT( pIn,  mIn);
            }

            virtual void bstback_GAMBIT(const Pythia8::Abstract__Vec4& pIn) {};
            void bstback(const Pythia8::Abstract__Vec4& pIn)
            {
                bstback_GAMBIT( pIn);
            }

            virtual void bstback_GAMBIT(const Pythia8::Abstract__Vec4& pIn, double mIn) {};
            void bstback(const Pythia8::Abstract__Vec4& pIn, double mIn)
            {
                bstback_GAMBIT( pIn,  mIn);
            }

            virtual void rotbst_GAMBIT(const Pythia8::Abstract__RotBstMatrix& M) {};
            void rotbst(const Pythia8::Abstract__RotBstMatrix& M)
            {
                rotbst_GAMBIT( M);
            }
            // UNKNOWN: OperatorMethod
            // UNKNOWN: OperatorMethod
            // UNKNOWN: OperatorMethod
            // UNKNOWN: OperatorMethod
            // UNKNOWN: OperatorMethod

        public:
            Vec4* downcast()
            {
                return reinterpret_cast<Vec4*>(this);
            }
    };
}
#pragma GCC diagnostic pop

