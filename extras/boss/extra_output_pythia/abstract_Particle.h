// Forward declarations:
class Pythia8::Abstract__Particle;
class Pythia8::Particle;
class Pythia8::Abstract__ParticleData;
class Pythia8::ParticleData;
class Pythia8::Abstract__Vec4;
class Pythia8::Vec4;
class Pythia8::Abstract__RotBstMatrix;
class Pythia8::RotBstMatrix;
class Pythia8::Abstract__Event;
class Pythia8::Event;
class Pythia8::Abstract__Pythia;
class Pythia8::Pythia;


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

            virtual void id_GAMBIT(int idIn) {};
            void id(int idIn)
            {
                id_GAMBIT( idIn);
            }

            virtual void status_GAMBIT(int statusIn) {};
            void status(int statusIn)
            {
                status_GAMBIT( statusIn);
            }

            virtual void statusPos_GAMBIT() {};
            void statusPos()
            {
                statusPos_GAMBIT();
            }

            virtual void statusNeg_GAMBIT() {};
            void statusNeg()
            {
                statusNeg_GAMBIT();
            }

            virtual void statusCode_GAMBIT(int statusIn) {};
            void statusCode(int statusIn)
            {
                statusCode_GAMBIT( statusIn);
            }

            virtual void mother1_GAMBIT(int mother1In) {};
            void mother1(int mother1In)
            {
                mother1_GAMBIT( mother1In);
            }

            virtual void mother2_GAMBIT(int mother2In) {};
            void mother2(int mother2In)
            {
                mother2_GAMBIT( mother2In);
            }

            virtual void mothers_GAMBIT(int mother1In, int mother2In) {};
            void mothers(int mother1In, int mother2In)
            {
                mothers_GAMBIT( mother1In,  mother2In);
            }

            virtual void daughter1_GAMBIT(int daughter1In) {};
            void daughter1(int daughter1In)
            {
                daughter1_GAMBIT( daughter1In);
            }

            virtual void daughter2_GAMBIT(int daughter2In) {};
            void daughter2(int daughter2In)
            {
                daughter2_GAMBIT( daughter2In);
            }

            virtual void daughters_GAMBIT(int daughter1In, int daughter2In) {};
            void daughters(int daughter1In, int daughter2In)
            {
                daughters_GAMBIT( daughter1In,  daughter2In);
            }

            virtual void col_GAMBIT(int colIn) {};
            void col(int colIn)
            {
                col_GAMBIT( colIn);
            }

            virtual void acol_GAMBIT(int acolIn) {};
            void acol(int acolIn)
            {
                acol_GAMBIT( acolIn);
            }

            virtual void cols_GAMBIT(int colIn, int acolIn) {};
            void cols(int colIn, int acolIn)
            {
                cols_GAMBIT( colIn,  acolIn);
            }

            virtual void p_GAMBIT(Pythia8::Abstract__Vec4 pIn) {};
            void p(Pythia8::Abstract__Vec4 pIn)
            {
                p_GAMBIT( pIn);
            }

            virtual void p_GAMBIT(double pxIn, double pyIn, double pzIn, double eIn) {};
            void p(double pxIn, double pyIn, double pzIn, double eIn)
            {
                p_GAMBIT( pxIn,  pyIn,  pzIn,  eIn);
            }

            virtual void px_GAMBIT(double pxIn) {};
            void px(double pxIn)
            {
                px_GAMBIT( pxIn);
            }

            virtual void py_GAMBIT(double pyIn) {};
            void py(double pyIn)
            {
                py_GAMBIT( pyIn);
            }

            virtual void pz_GAMBIT(double pzIn) {};
            void pz(double pzIn)
            {
                pz_GAMBIT( pzIn);
            }

            virtual void e_GAMBIT(double eIn) {};
            void e(double eIn)
            {
                e_GAMBIT( eIn);
            }

            virtual void m_GAMBIT(double mIn) {};
            void m(double mIn)
            {
                m_GAMBIT( mIn);
            }

            virtual void scale_GAMBIT(double scaleIn) {};
            void scale(double scaleIn)
            {
                scale_GAMBIT( scaleIn);
            }

            virtual void pol_GAMBIT(double polIn) {};
            void pol(double polIn)
            {
                pol_GAMBIT( polIn);
            }

            virtual void vProd_GAMBIT(Pythia8::Abstract__Vec4 vProdIn) {};
            void vProd(Pythia8::Abstract__Vec4 vProdIn)
            {
                vProd_GAMBIT( vProdIn);
            }

            virtual void vProd_GAMBIT(double xProdIn, double yProdIn, double zProdIn, double tProdIn) {};
            void vProd(double xProdIn, double yProdIn, double zProdIn, double tProdIn)
            {
                vProd_GAMBIT( xProdIn,  yProdIn,  zProdIn,  tProdIn);
            }

            virtual void xProd_GAMBIT(double xProdIn) {};
            void xProd(double xProdIn)
            {
                xProd_GAMBIT( xProdIn);
            }

            virtual void yProd_GAMBIT(double yProdIn) {};
            void yProd(double yProdIn)
            {
                yProd_GAMBIT( yProdIn);
            }

            virtual void zProd_GAMBIT(double zProdIn) {};
            void zProd(double zProdIn)
            {
                zProd_GAMBIT( zProdIn);
            }

            virtual void tProd_GAMBIT(double tProdIn) {};
            void tProd(double tProdIn)
            {
                tProd_GAMBIT( tProdIn);
            }

            virtual void tau_GAMBIT(double tauIn) {};
            void tau(double tauIn)
            {
                tau_GAMBIT( tauIn);
            }

            virtual int* id_GAMBIT() {};
            int* id()
            {
                return id_GAMBIT();
            }

            virtual int* status_GAMBIT() {};
            int* status()
            {
                return status_GAMBIT();
            }

            virtual int* mother1_GAMBIT() {};
            int* mother1()
            {
                return mother1_GAMBIT();
            }

            virtual int* mother2_GAMBIT() {};
            int* mother2()
            {
                return mother2_GAMBIT();
            }

            virtual int* daughter1_GAMBIT() {};
            int* daughter1()
            {
                return daughter1_GAMBIT();
            }

            virtual int* daughter2_GAMBIT() {};
            int* daughter2()
            {
                return daughter2_GAMBIT();
            }

            virtual int* col_GAMBIT() {};
            int* col()
            {
                return col_GAMBIT();
            }

            virtual int* acol_GAMBIT() {};
            int* acol()
            {
                return acol_GAMBIT();
            }

            virtual Pythia8::Vec4* p_GAMBIT() {};
            Pythia8::Abstract__Vec4* p()
            {
                return reinterpret_cast<Pythia8::Abstract__Vec4*>(p_GAMBIT());
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

            virtual double* m_GAMBIT() {};
            double* m()
            {
                return m_GAMBIT();
            }

            virtual double* scale_GAMBIT() {};
            double* scale()
            {
                return scale_GAMBIT();
            }

            virtual double* pol_GAMBIT() {};
            double* pol()
            {
                return pol_GAMBIT();
            }

            virtual bool* hasVertex_GAMBIT() {};
            bool* hasVertex()
            {
                return hasVertex_GAMBIT();
            }

            virtual Pythia8::Vec4* vProd_GAMBIT() {};
            Pythia8::Abstract__Vec4* vProd()
            {
                return reinterpret_cast<Pythia8::Abstract__Vec4*>(vProd_GAMBIT());
            }

            virtual double* xProd_GAMBIT() {};
            double* xProd()
            {
                return xProd_GAMBIT();
            }

            virtual double* yProd_GAMBIT() {};
            double* yProd()
            {
                return yProd_GAMBIT();
            }

            virtual double* zProd_GAMBIT() {};
            double* zProd()
            {
                return zProd_GAMBIT();
            }

            virtual double* tProd_GAMBIT() {};
            double* tProd()
            {
                return tProd_GAMBIT();
            }

            virtual double* tau_GAMBIT() {};
            double* tau()
            {
                return tau_GAMBIT();
            }

            virtual int* idAbs_GAMBIT() {};
            int* idAbs()
            {
                return idAbs_GAMBIT();
            }

            virtual int* statusAbs_GAMBIT() {};
            int* statusAbs()
            {
                return statusAbs_GAMBIT();
            }

            virtual bool* isFinal_GAMBIT() {};
            bool* isFinal()
            {
                return isFinal_GAMBIT();
            }

            virtual bool* isRescatteredIncoming_GAMBIT() {};
            bool* isRescatteredIncoming()
            {
                return isRescatteredIncoming_GAMBIT();
            }

            virtual double* m2_GAMBIT() {};
            double* m2()
            {
                return m2_GAMBIT();
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

            virtual double* eCalc_GAMBIT() {};
            double* eCalc()
            {
                return eCalc_GAMBIT();
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

            virtual double* mT_GAMBIT() {};
            double* mT()
            {
                return mT_GAMBIT();
            }

            virtual double* mT2_GAMBIT() {};
            double* mT2()
            {
                return mT2_GAMBIT();
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

            virtual double* y_GAMBIT() {};
            double* y()
            {
                return y_GAMBIT();
            }

            virtual double* eta_GAMBIT() {};
            double* eta()
            {
                return eta_GAMBIT();
            }

            virtual Pythia8::Vec4* vDec_GAMBIT() {};
            Pythia8::Abstract__Vec4* vDec()
            {
                return reinterpret_cast<Pythia8::Abstract__Vec4*>(vDec_GAMBIT());
            }

            virtual double* xDec_GAMBIT() {};
            double* xDec()
            {
                return xDec_GAMBIT();
            }

            virtual double* yDec_GAMBIT() {};
            double* yDec()
            {
                return yDec_GAMBIT();
            }

            virtual double* zDec_GAMBIT() {};
            double* zDec()
            {
                return zDec_GAMBIT();
            }

            virtual double* tDec_GAMBIT() {};
            double* tDec()
            {
                return tDec_GAMBIT();
            }

            virtual int* index_GAMBIT() {};
            int* index()
            {
                return index_GAMBIT();
            }

            virtual int* statusHepMC_GAMBIT() {};
            int* statusHepMC()
            {
                return statusHepMC_GAMBIT();
            }

            virtual int* iTopCopy_GAMBIT() {};
            int* iTopCopy()
            {
                return iTopCopy_GAMBIT();
            }

            virtual int* iBotCopy_GAMBIT() {};
            int* iBotCopy()
            {
                return iBotCopy_GAMBIT();
            }

            virtual int* iTopCopyId_GAMBIT() {};
            int* iTopCopyId()
            {
                return iTopCopyId_GAMBIT();
            }

            virtual int* iBotCopyId_GAMBIT() {};
            int* iBotCopyId()
            {
                return iBotCopyId_GAMBIT();
            }

            virtual bool* isAncestor_GAMBIT(int iAncestor) {};
            bool* isAncestor(int iAncestor)
            {
                return isAncestor_GAMBIT( iAncestor);
            }

            virtual bool* undoDecay_GAMBIT() {};
            bool* undoDecay()
            {
                return undoDecay_GAMBIT();
            }

            virtual std::string* name_GAMBIT() {};
            std::string* name()
            {
                return name_GAMBIT();
            }

            virtual std::string* nameWithStatus_GAMBIT(int maxLen) {};
            std::string* nameWithStatus(int maxLen)
            {
                return nameWithStatus_GAMBIT( maxLen);
            }

            virtual int* spinType_GAMBIT() {};
            int* spinType()
            {
                return spinType_GAMBIT();
            }

            virtual int* chargeType_GAMBIT() {};
            int* chargeType()
            {
                return chargeType_GAMBIT();
            }

            virtual double* charge_GAMBIT() {};
            double* charge()
            {
                return charge_GAMBIT();
            }

            virtual bool* isCharged_GAMBIT() {};
            bool* isCharged()
            {
                return isCharged_GAMBIT();
            }

            virtual bool* isNeutral_GAMBIT() {};
            bool* isNeutral()
            {
                return isNeutral_GAMBIT();
            }

            virtual int* colType_GAMBIT() {};
            int* colType()
            {
                return colType_GAMBIT();
            }

            virtual double* m0_GAMBIT() {};
            double* m0()
            {
                return m0_GAMBIT();
            }

            virtual double* mWidth_GAMBIT() {};
            double* mWidth()
            {
                return mWidth_GAMBIT();
            }

            virtual double* mMin_GAMBIT() {};
            double* mMin()
            {
                return mMin_GAMBIT();
            }

            virtual double* mMax_GAMBIT() {};
            double* mMax()
            {
                return mMax_GAMBIT();
            }

            virtual double* mSel_GAMBIT() {};
            double* mSel()
            {
                return mSel_GAMBIT();
            }

            virtual double* constituentMass_GAMBIT() {};
            double* constituentMass()
            {
                return constituentMass_GAMBIT();
            }

            virtual double* tau0_GAMBIT() {};
            double* tau0()
            {
                return tau0_GAMBIT();
            }

            virtual bool* mayDecay_GAMBIT() {};
            bool* mayDecay()
            {
                return mayDecay_GAMBIT();
            }

            virtual bool* canDecay_GAMBIT() {};
            bool* canDecay()
            {
                return canDecay_GAMBIT();
            }

            virtual bool* doExternalDecay_GAMBIT() {};
            bool* doExternalDecay()
            {
                return doExternalDecay_GAMBIT();
            }

            virtual bool* isResonance_GAMBIT() {};
            bool* isResonance()
            {
                return isResonance_GAMBIT();
            }

            virtual bool* isVisible_GAMBIT() {};
            bool* isVisible()
            {
                return isVisible_GAMBIT();
            }

            virtual bool* isLepton_GAMBIT() {};
            bool* isLepton()
            {
                return isLepton_GAMBIT();
            }

            virtual bool* isQuark_GAMBIT() {};
            bool* isQuark()
            {
                return isQuark_GAMBIT();
            }

            virtual bool* isGluon_GAMBIT() {};
            bool* isGluon()
            {
                return isGluon_GAMBIT();
            }

            virtual bool* isDiquark_GAMBIT() {};
            bool* isDiquark()
            {
                return isDiquark_GAMBIT();
            }

            virtual bool* isParton_GAMBIT() {};
            bool* isParton()
            {
                return isParton_GAMBIT();
            }

            virtual bool* isHadron_GAMBIT() {};
            bool* isHadron()
            {
                return isHadron_GAMBIT();
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

            virtual void rescale5_GAMBIT(double fac) {};
            void rescale5(double fac)
            {
                rescale5_GAMBIT( fac);
            }

            virtual void rot_GAMBIT(double thetaIn, double phiIn) {};
            void rot(double thetaIn, double phiIn)
            {
                rot_GAMBIT( thetaIn,  phiIn);
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

            virtual void offsetHistory_GAMBIT(int minMother, int addMother, int minDaughter, int addDaughter) {};
            void offsetHistory(int minMother, int addMother, int minDaughter, int addDaughter)
            {
                offsetHistory_GAMBIT( minMother,  addMother,  minDaughter,  addDaughter);
            }

            virtual void offsetCol_GAMBIT(int addCol) {};
            void offsetCol(int addCol)
            {
                offsetCol_GAMBIT( addCol);
            }

        public:
            Particle* downcast()
            {
                return reinterpret_cast<Particle*>(this);
            }
    };
}

