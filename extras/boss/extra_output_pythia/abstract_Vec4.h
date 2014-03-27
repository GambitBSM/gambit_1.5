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

namespace Pythia8
{
    class Abstract__Vec4
    {
        private:
        public:
            // UNKNOWN: OperatorMethod
            // UNKNOWN: OperatorMethod

            virtual void reset_GAMBIT() {};
            void reset()
            {
                reset_GAMBIT();
            }

            virtual void p_GAMBIT(double xIn, double yIn, double zIn, double tIn) {};
            void p(double xIn, double yIn, double zIn, double tIn)
            {
                p_GAMBIT( xIn,  yIn,  zIn,  tIn);
            }

            virtual void p_GAMBIT(Pythia8::Abstract__Vec4 pIn) {};
            void p(Pythia8::Abstract__Vec4 pIn)
            {
                p_GAMBIT( pIn);
            }

            virtual void px_GAMBIT(double xIn) {};
            void px(double xIn)
            {
                px_GAMBIT( xIn);
            }

            virtual void py_GAMBIT(double yIn) {};
            void py(double yIn)
            {
                py_GAMBIT( yIn);
            }

            virtual void pz_GAMBIT(double zIn) {};
            void pz(double zIn)
            {
                pz_GAMBIT( zIn);
            }

            virtual void e_GAMBIT(double tIn) {};
            void e(double tIn)
            {
                e_GAMBIT( tIn);
            }

            virtual double* px_GAMBIT() {};
            double* px()
            {
                return px_GAMBIT();
            }

            virtual double* py_GAMBIT() {};
            double* py()
            {
                return py_GAMBIT();
            }

            virtual double* pz_GAMBIT() {};
            double* pz()
            {
                return pz_GAMBIT();
            }

            virtual double* e_GAMBIT() {};
            double* e()
            {
                return e_GAMBIT();
            }

            virtual double* mCalc_GAMBIT() {};
            double* mCalc()
            {
                return mCalc_GAMBIT();
            }

            virtual double* m2Calc_GAMBIT() {};
            double* m2Calc()
            {
                return m2Calc_GAMBIT();
            }

            virtual double* pT_GAMBIT() {};
            double* pT()
            {
                return pT_GAMBIT();
            }

            virtual double* pT2_GAMBIT() {};
            double* pT2()
            {
                return pT2_GAMBIT();
            }

            virtual double* pAbs_GAMBIT() {};
            double* pAbs()
            {
                return pAbs_GAMBIT();
            }

            virtual double* pAbs2_GAMBIT() {};
            double* pAbs2()
            {
                return pAbs2_GAMBIT();
            }

            virtual double* eT_GAMBIT() {};
            double* eT()
            {
                return eT_GAMBIT();
            }

            virtual double* eT2_GAMBIT() {};
            double* eT2()
            {
                return eT2_GAMBIT();
            }

            virtual double* theta_GAMBIT() {};
            double* theta()
            {
                return theta_GAMBIT();
            }

            virtual double* phi_GAMBIT() {};
            double* phi()
            {
                return phi_GAMBIT();
            }

            virtual double* thetaXZ_GAMBIT() {};
            double* thetaXZ()
            {
                return thetaXZ_GAMBIT();
            }

            virtual double* pPos_GAMBIT() {};
            double* pPos()
            {
                return pPos_GAMBIT();
            }

            virtual double* pNeg_GAMBIT() {};
            double* pNeg()
            {
                return pNeg_GAMBIT();
            }

            virtual double* rap_GAMBIT() {};
            double* rap()
            {
                return rap_GAMBIT();
            }

            virtual double* eta_GAMBIT() {};
            double* eta()
            {
                return eta_GAMBIT();
            }

            virtual void rescale3_GAMBIT(double fac) {};
            void rescale3(double fac)
            {
                rescale3_GAMBIT( fac);
            }

            virtual void rescale4_GAMBIT(double fac) {};
            void rescale4(double fac)
            {
                rescale4_GAMBIT( fac);
            }

            virtual void flip3_GAMBIT() {};
            void flip3()
            {
                flip3_GAMBIT();
            }

            virtual void flip4_GAMBIT() {};
            void flip4()
            {
                flip4_GAMBIT();
            }

            virtual void rot_GAMBIT(double thetaIn, double phiIn) {};
            void rot(double thetaIn, double phiIn)
            {
                rot_GAMBIT( thetaIn,  phiIn);
            }

            virtual void rotaxis_GAMBIT(double phiIn, double nx, double ny, double nz) {};
            void rotaxis(double phiIn, double nx, double ny, double nz)
            {
                rotaxis_GAMBIT( phiIn,  nx,  ny,  nz);
            }

            virtual void rotaxis_GAMBIT(double phiIn, const Pythia8::Abstract__Vec4& n) {};
            void rotaxis(double phiIn, const Pythia8::Abstract__Vec4& n)
            {
                rotaxis_GAMBIT( phiIn,  n);
            }

            virtual void bst_GAMBIT(double betaX, double betaY, double betaZ) {};
            void bst(double betaX, double betaY, double betaZ)
            {
                bst_GAMBIT( betaX,  betaY,  betaZ);
            }

            virtual void bst_GAMBIT(double betaX, double betaY, double betaZ, double gamma) {};
            void bst(double betaX, double betaY, double betaZ, double gamma)
            {
                bst_GAMBIT( betaX,  betaY,  betaZ,  gamma);
            }

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

