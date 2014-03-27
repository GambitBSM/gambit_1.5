
        public:
            void rot_GAMBIT(double arg_1, double arg_2)
            {
                rot(arg_1,arg_2);
            }
            void rot_GAMBIT(const Pythia8::Abstract__Vec4& p)
            {
                const Pythia8::Vec4& _temp_p = *(reinterpret_cast<const Pythia8::Vec4*>(&p));
                rot(_temp_p);
            }
            void bst_GAMBIT(double arg_1, double arg_2, double arg_3)
            {
                bst(arg_1,arg_2,arg_3);
            }
            void bst_GAMBIT(const Pythia8::Abstract__Vec4& arg_1)
            {
                const Pythia8::Vec4& _temp_arg_1 = *(reinterpret_cast<const Pythia8::Vec4*>(&arg_1));
                bst(_temp_arg_1);
            }
            void bstback_GAMBIT(const Pythia8::Abstract__Vec4& arg_1)
            {
                const Pythia8::Vec4& _temp_arg_1 = *(reinterpret_cast<const Pythia8::Vec4*>(&arg_1));
                bstback(_temp_arg_1);
            }
            void bst_GAMBIT(const Pythia8::Abstract__Vec4& arg_1, const Pythia8::Abstract__Vec4& arg_2)
            {
                const Pythia8::Vec4& _temp_arg_1 = *(reinterpret_cast<const Pythia8::Vec4*>(&arg_1));
                const Pythia8::Vec4& _temp_arg_2 = *(reinterpret_cast<const Pythia8::Vec4*>(&arg_2));
                bst(_temp_arg_1,_temp_arg_2);
            }
            void toCMframe_GAMBIT(const Pythia8::Abstract__Vec4& arg_1, const Pythia8::Abstract__Vec4& arg_2)
            {
                const Pythia8::Vec4& _temp_arg_1 = *(reinterpret_cast<const Pythia8::Vec4*>(&arg_1));
                const Pythia8::Vec4& _temp_arg_2 = *(reinterpret_cast<const Pythia8::Vec4*>(&arg_2));
                toCMframe(_temp_arg_1,_temp_arg_2);
            }
            void fromCMframe_GAMBIT(const Pythia8::Abstract__Vec4& arg_1, const Pythia8::Abstract__Vec4& arg_2)
            {
                const Pythia8::Vec4& _temp_arg_1 = *(reinterpret_cast<const Pythia8::Vec4*>(&arg_1));
                const Pythia8::Vec4& _temp_arg_2 = *(reinterpret_cast<const Pythia8::Vec4*>(&arg_2));
                fromCMframe(_temp_arg_1,_temp_arg_2);
            }
            void rotbst_GAMBIT(const Pythia8::Abstract__RotBstMatrix& arg_1)
            {
                const Pythia8::RotBstMatrix& _temp_arg_1 = *(reinterpret_cast<const Pythia8::RotBstMatrix*>(&arg_1));
                rotbst(_temp_arg_1);
            }
            void invert_GAMBIT()
            {
                invert();
            }
            void reset_GAMBIT()
            {
                reset();
            }
            double* deviation_GAMBIT()
            {
                return new double(deviation());
            }
 : public virtual Abstract__RotBstMatrix} 
#include "Pythia8/abstract_RotBstMatrix.h"
namespace Pythia8 { 

        public:
            void reset_GAMBIT()
            {
                reset();
            }
            void p_GAMBIT(double xIn, double yIn, double zIn, double tIn)
            {
                p(xIn,yIn,zIn,tIn);
            }
            void p_GAMBIT(Pythia8::Abstract__Vec4 pIn)
            {
                Pythia8::Vec4 _temp_pIn = *(reinterpret_cast<Pythia8::Vec4*>(&pIn));
                p(_temp_pIn);
            }
            void px_GAMBIT(double xIn)
            {
                px(xIn);
            }
            void py_GAMBIT(double yIn)
            {
            } 
#include "Pythia8/abstract_Vec4.h"
namespace Pythia8 { 
    py(yIn : public virtual Abstract__Vec4);
            }
            void pz_GAMBIT(double zIn)
            {
                pz(zIn);
            }
            void e_GAMBIT(double tIn)
            {
                e(tIn);
            }
            double* px_GAMBIT()
            {
                return new double(px());
            }
            double* py_GAMBIT()
            {
                return new double(py());
            }
            double* pz_GAMBIT()
            {
                return new double(pz());
            }
            double* e_GAMBIT()
            {
                return new double(e());
            }
            double* mCalc_GAMBIT()
            {
                return new double(mCalc());
            }
            double* m2Calc_GAMBIT()
            {
                return new double(m2Calc());
            }
            double* pT_GAMBIT()
            {
                return new double(pT());
            }
            double* pT2_GAMBIT()
            {
                return new double(pT2());
            }
            double* pAbs_GAMBIT()
            {
                return new double(pAbs());
            }
            double* pAbs2_GAMBIT()
            {
                return new double(pAbs2());
            }
            double* eT_GAMBIT()
            {
                return new double(eT());
            }
            double* eT2_GAMBIT()
            {
                return new double(eT2());
            }
            double* theta_GAMBIT()
            {
                return new double(theta());
            }
            double* phi_GAMBIT()
            {
                return new double(phi());
            }
            double* thetaXZ_GAMBIT()
            {
                return new double(thetaXZ());
            }
            double* pPos_GAMBIT()
            {
                return new double(pPos());
            }
            double* pNeg_GAMBIT()
            {
                return new double(pNeg());
            }
            double* rap_GAMBIT()
            {
                return new double(rap());
            }
            double* eta_GAMBIT()
            {
                return new double(eta());
            }
            void rescale3_GAMBIT(double fac)
            {
                rescale3(fac);
            }
            void rescale4_GAMBIT(double fac)
            {
                rescale4(fac);
            }
            void flip3_GAMBIT()
            {
                flip3();
            }
            void flip4_GAMBIT()
            {
                flip4();
            }
            void rot_GAMBIT(double thetaIn, double phiIn)
            {
                rot(thetaIn,phiIn);
            }
            void rotaxis_GAMBIT(double phiIn, double nx, double ny, double nz)
            {
                rotaxis(phiIn,nx,ny,nz);
            }
            void rotaxis_GAMBIT(double phiIn, const Pythia8::Abstract__Vec4& n)
            {
                const Pythia8::Vec4& _temp_n = *(reinterpret_cast<const Pythia8::Vec4*>(&n));
                rotaxis(phiIn,_temp_n);
            }
            void bst_GAMBIT(double betaX, double betaY, double betaZ)
            {
                bst(betaX,betaY,betaZ);
            }
            void bst_GAMBIT(double betaX, double betaY, double betaZ, double gamma)
            {
                bst(betaX,betaY,betaZ,gamma);
            }
            void bst_GAMBIT(const Pythia8::Abstract__Vec4& pIn)
            {
                const Pythia8::Vec4& _temp_pIn = *(reinterpret_cast<const Pythia8::Vec4*>(&pIn));
                bst(_temp_pIn);
            }
            void bst_GAMBIT(const Pythia8::Abstract__Vec4& pIn, double mIn)
            {
                const Pythia8::Vec4& _temp_pIn = *(reinterpret_cast<const Pythia8::Vec4*>(&pIn));
                bst(_temp_pIn,mIn);
            }
            void bstback_GAMBIT(const Pythia8::Abstract__Vec4& pIn)
            {
                const Pythia8::Vec4& _temp_pIn = *(reinterpret_cast<const Pythia8::Vec4*>(&pIn));
                bstback(_temp_pIn);
            }
            void bstback_GAMBIT(const Pythia8::Abstract__Vec4& pIn, double mIn)
            {
                const Pythia8::Vec4& _temp_pIn = *(reinterpret_cast<const Pythia8::Vec4*>(&pIn));
                bstback(_temp_pIn,mIn);
            }
            void rotbst_GAMBIT(const Pythia8::Abstract__RotBstMatrix& M)
            {
                const Pythia8::RotBstMatrix& _temp_M = *(reinterpret_cast<const Pythia8::RotBstMatrix*>(&M));
                rotbst(_temp_M);
            }
