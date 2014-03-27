
        public:
            void init_GAMBIT(std::string headerIn, Pythia8::Abstract__ParticleData* particleDataPtrIn, int startColTagIn)
            {
                Pythia8::ParticleData* _temp_particleDataPtrIn = reinterpret_cast<Pythia8::ParticleData*>(particleDataPtrIn);
                init(headerIn,_temp_particleDataPtrIn,startColTagIn);
            }
            void clear_GAMBIT()
            {
                clear();
            }
            void reset_GAMBIT()
            {
                reset();
            }
            Pythia8::Particle*& front_GAMBIT()
            {
                Pythia8::Particle* temp_p = new Pythia8::Particle(front());
                Pythia8::Particle*& temp_pr = temp_p;
                return temp_pr;
            }
            Pythia8::Particle*& at_GAMBIT(int i)
            {
                Pythia8::Particle* temp_p = new Pythia8::Particle(at(i));
                Pythia8::Particle*& temp_pr = temp_p;
                return te} 
#include "Pythia8/abstract_Particle.h"
namespace Pythia8 { 
mp_pr;
        : public virtual Abstract__Particle     }
            Pythia8::Particle*& back_GAMBIT()
            {
                Pythia8::Particle* temp_p = new Pythia8::Particle(back());
                Pythia8::Particle*& temp_pr = temp_p;
                return temp_pr;
            }
            int* size_GAMBIT()
            {
                return new int(size());
            }
            int* append_GAMBIT(Pythia8::Abstract__Particle entryIn)
            {
                Pythia8::Particle _temp_entryIn = *(reinterpret_cast<Pythia8::Particle*>(&entryIn));
                return new int(append(_temp_entryIn));
            }
            int* append_GAMBIT(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn, double polIn)
            {
                return new int(append(id,status,mother1,mother2,daughter1,daughter2,col,acol,px,py,pz,e,m,scaleIn,polIn));
            }
            int* append_GAMBIT(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract__Vec4 p, double m, double scaleIn, double polIn)
            {
                Pythia8::Vec4 _temp_p = *(reinterpret_cast<Pythia8::Vec4*>(&p));
                return new int(append(id,status,mother1,mother2,daughter1,daughter2,col,acol,_temp_p,m,scaleIn,polIn));
            }
            int* append_GAMBIT(int id, int status, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn, double polIn)
            {
                return new int(append(id,status,col,acol,px,py,pz,e,m,scaleIn,polIn));
            }
            int* append_GAMBIT(int id, int status, int col, int acol, Pythia8::Abstract__Vec4 p, double m, double scaleIn, double polIn)
            {
                Pythia8::Vec4 _temp_p = *(reinterpret_cast<Pythia8::Vec4*>(&p));
                return new int(append(id,status,col,acol,_temp_p,m,scaleIn,polIn));
            }
            void setEvtPtr_GAMBIT(int iSet)
            {
                setEvtPtr(iSet);
            }
            int* copy_GAMBIT(int iCopy, int newStatus)
            {
                return new int(copy(iCopy,newStatus));
            }
            void list_GAMBIT()
            {
                list();
            }
            void list_GAMBIT(std::ostream& os)
            {
                list(os);
            }
            void list_GAMBIT(bool showScaleAndVertex, bool showMothersAndDaughters)
            {
                list(showScaleAndVertex,showMothersAndDaughters);
            }
            void list_GAMBIT(bool showScaleAndVertex, bool showMothersAndDaughters, std::ostream& os)
            {
                list(showScaleAndVertex,showMothersAndDaughters,os);
            }
            void popBack_GAMBIT(int nRemove)
            {
                popBack(nRemove);
            }
            void remove_GAMBIT(int iFirst, int iLast)
            {
                remove(iFirst,iLast);
            }
            bool* undoDecay_GAMBIT(int i)
            {
                return new bool(undoDecay(i));
            }
            void restorePtrs_GAMBIT()
            {
                restorePtrs();
            }
            void saveSize_GAMBIT()
            {
                saveSize();
            }
            void restoreSize_GAMBIT()
            {
                restoreSize();
            }
            int* savedSizeValue_GAMBIT()
            {
                return new int(savedSizeValue());
            }
            void initColTag_GAMBIT(int colTag)
            {
                initColTag(colTag);
            }
            int* lastColTag_GAMBIT()
            {
                return new int(lastColTag());
            }
            int* nextColTag_GAMBIT()
            {
                return new int(nextColTag());
            }
            void scale_GAMBIT(double scaleIn)
            {
                scale(scaleIn);
            }
            double* scale_GAMBIT()
            {
                return new double(scale());
            }
            void scaleSecond_GAMBIT(double scaleSecondIn)
            {
                scaleSecond(scaleSecondIn);
            }
            double* scaleSecond_GAMBIT()
            {
                return new double(scaleSecond());
            }
            int* statusHepMC_GAMBIT(int i)
            {
                return new int(statusHepMC(i));
            }
            int* iTopCopy_GAMBIT(int i)
            {
                return new int(iTopCopy(i));
            }
            int* iBotCopy_GAMBIT(int i)
            {
                return new int(iBotCopy(i));
            }
            int* iTopCopyId_GAMBIT(int i)
            {
                return new int(iTopCopyId(i));
            }
            int* iBotCopyId_GAMBIT(int i)
            {
                return new int(iBotCopyId(i));
            }
            bool* isAncestor_GAMBIT(int i, int iAncestor)
            {
                return new bool(isAncestor(i,iAncestor));
            }
            void rot_GAMBIT(double theta, double phi)
            {
                rot(theta,phi);
            }
            void bst_GAMBIT(double betaX, double betaY, double betaZ)
            {
                bst(betaX,betaY,betaZ);
            }
            void bst_GAMBIT(double betaX, double betaY, double betaZ, double gamma)
            {
                bst(betaX,betaY,betaZ,gamma);
            }
            void bst_GAMBIT(const Pythia8::Abstract__Vec4& vec)
            {
                const Pythia8::Vec4& _temp_vec = *(reinterpret_cast<const Pythia8::Vec4*>(&vec));
                bst(_temp_vec);
            }
            void rotbst_GAMBIT(const Pythia8::Abstract__RotBstMatrix& M)
            {
                const Pythia8::RotBstMatrix& _temp_M = *(reinterpret_cast<const Pythia8::RotBstMatrix*>(&M));
                rotbst(_temp_M);
            }
            void clearJunctions_GAMBIT()
            {
                clearJunctions();
            }
            void appendJunction_GAMBIT(int kind, int col0, int col1, int col2)
            {
                appendJunction(kind,col0,col1,col2);
            }
            int* sizeJunction_GAMBIT()
            {
                return new int(sizeJunction());
            }
            bool* remainsJunction_GAMBIT(int i)
            {
                return new bool(remainsJunction(i));
            }
            void remainsJunction_GAMBIT(int i, bool remainsIn)
            {
                remainsJunction(i,remainsIn);
            }
            int* kindJunction_GAMBIT(int i)
            {
                return new int(kindJunction(i));
            }
            int* colJunction_GAMBIT(int i, int j)
            {
                return new int(colJunction(i,j));
            }
            void colJunction_GAMBIT(int i, int j, int colIn)
            {
                colJunction(i,j,colIn);
            }
            int* endColJunction_GAMBIT(int i, int j)
            {
                return new int(endColJunction(i,j));
            }
            void endColJunction_GAMBIT(int i, int j, int endColIn)
            {
                endColJunction(i,j,endColIn);
            }
            int* statusJunction_GAMBIT(int i, int j)
            {
                return new int(statusJunction(i,j));
            }
            void statusJunction_GAMBIT(int i, int j, int statusIn)
            {
                statusJunction(i,j,statusIn);
            }
            void eraseJunction_GAMBIT(int i)
            {
                eraseJunction(i);
            }
            void saveJunctionSize_GAMBIT()
            {
                saveJunctionSize();
            }
            void restoreJunctionSize_GAMBIT()
            {
                restoreJunctionSize();
            }
            void listJunctions_GAMBIT(std::ostream& os)
            {
                listJunctions(os);
            }
 : public virtual Abstract__Event} 
#include "Pythia8/abstract_Event.h"
namespace Pythia8 { 

        public:
            void setEvtPtr_GAMBIT(Pythia8::Abstract__Event* evtPtrIn)
            {
                Pythia8::Event* _temp_evtPtrIn = reinterpret_cast<Pythia8::Event*>(evtPtrIn);
                setEvtPtr(_temp_evtPtrIn);
            }
            void id_GAMBIT(int idIn)
            {
                id(idIn);
            }
            void status_GAMBIT(int statusIn)
            {
                status(statusIn);
            }
            void statusPos_GAMBIT()
            {
                statusPos();
            }
            void statusNeg_GAMBIT()
            {
                statusNeg();
            }
            void statusCode_GAMBIT(int statusIn)
            {
                statusCode(statusIn);
            }
            void mother1_GAMBIT(int mother1In)
            {
                mother1(mother1In);
            }
            void mother2_GAMBIT(int mother2In)
            {
                mother2(mother2In);
            }
            void mothers_GAMBIT(int mother1In, int mother2In)
            {
                mothers(mother1In,mother2In);
            }
            void daughter1_GAMBIT(int daughter1In)
            {
                daughter1(daughter1In);
            }
            void daughter2_GAMBIT(int daughter2In)
            {
                daughter2(daughter2In);
            }
            void daughters_GAMBIT(int daughter1In, int daughter2In)
            {
                daughters(daughter1In,daughter2In);
            }
            void col_GAMBIT(int colIn)
            {
                col(colIn);
            }
            void acol_GAMBIT(int acolIn)
            {
                acol(acolIn);
            }
            void cols_GAMBIT(int colIn, int acolIn)
            {
                cols(colIn,acolIn);
            }
            void p_GAMBIT(Pythia8::Abstract__Vec4 pIn)
            {
                Pythia8::Vec4 _temp_pIn = *(reinterpret_cast<Pythia8::Vec4*>(&pIn));
                p(_temp_pIn);
            }
            void p_GAMBIT(double pxIn, double pyIn, double pzIn, double eIn)
            {
                p(pxIn,pyIn,pzIn,eIn);
            }
            void px_GAMBIT(double pxIn)
            {
                px(pxIn);
            }
            void py_GAMBIT(double pyIn)
            {
                py(pyIn);
            }
            void pz_GAMBIT(double pzIn)
            {
                pz(pzIn);
            }
            void e_GAMBIT(double eIn)
            {
                e(eIn);
            }
            void m_GAMBIT(double mIn)
            {
                m(mIn);
            }
            void scale_GAMBIT(double scaleIn)
            {
                scale(scaleIn);
            }
            void pol_GAMBIT(double polIn)
            {
                pol(polIn);
            }
            void vProd_GAMBIT(Pythia8::Abstract__Vec4 vProdIn)
            {
                Pythia8::Vec4 _temp_vProdIn = *(reinterpret_cast<Pythia8::Vec4*>(&vProdIn));
                vProd(_temp_vProdIn);
            }
            void vProd_GAMBIT(double xProdIn, double yProdIn, double zProdIn, double tProdIn)
            {
                vProd(xProdIn,yProdIn,zProdIn,tProdIn);
            }
            void xProd_GAMBIT(double xProdIn)
            {
                xProd(xProdIn);
            }
            void yProd_GAMBIT(double yProdIn)
            {
                yProd(yProdIn);
            }
            void zProd_GAMBIT(double zProdIn)
            {
                zProd(zProdIn);
            }
            void tProd_GAMBIT(double tProdIn)
            {
                tProd(tProdIn);
            }
            void tau_GAMBIT(double tauIn)
            {
                tau(tauIn);
            }
            int* id_GAMBIT()
            {
                return new int(id());
            }
            int* status_GAMBIT()
            {
                return new int(status());
            }
            int* mother1_GAMBIT()
            {
                return new int(mother1());
            }
            int* mother2_GAMBIT()
            {
                return new int(mother2());
            }
            int* daughter1_GAMBIT()
            {
                return new int(daughter1());
            }
            int* daughter2_GAMBIT()
            {
                return new int(daughter2());
            }
            int* col_GAMBIT()
            {
                return new int(col());
            }
            int* acol_GAMBIT()
            {
                return new int(acol());
            }
            Pythia8::Vec4* p_GAMBIT()
            {
                return new Pythia8::Vec4(p());
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
            double* m_GAMBIT()
            {
                return new double(m());
            }
            double* scale_GAMBIT()
            {
                return new double(scale());
            }
            double* pol_GAMBIT()
            {
                return new double(pol());
            }
            bool* hasVertex_GAMBIT()
            {
                return new bool(hasVertex());
            }
            Pythia8::Vec4* vProd_GAMBIT()
            {
                return new Pythia8::Vec4(vProd());
            }
            double* xProd_GAMBIT()
            {
                return new double(xProd());
            }
            double* yProd_GAMBIT()
            {
                return new double(yProd());
            }
            double* zProd_GAMBIT()
            {
                return new double(zProd());
            }
            double* tProd_GAMBIT()
            {
                return new double(tProd());
            }
            double* tau_GAMBIT()
            {
                return new double(tau());
            }
            int* idAbs_GAMBIT()
            {
                return new int(idAbs());
            }
            int* statusAbs_GAMBIT()
            {
                return new int(statusAbs());
            }
            bool* isFinal_GAMBIT()
            {
                return new bool(isFinal());
            }
            bool* isRescatteredIncoming_GAMBIT()
            {
                return new bool(isRescatteredIncoming());
            }
            double* m2_GAMBIT()
            {
                return new double(m2());
            }
            double* mCalc_GAMBIT()
            {
                return new double(mCalc());
            }
            double* m2Calc_GAMBIT()
            {
                return new double(m2Calc());
            }
            double* eCalc_GAMBIT()
            {
                return new double(eCalc());
            }
            double* pT_GAMBIT()
            {
                return new double(pT());
            }
            double* pT2_GAMBIT()
            {
                return new double(pT2());
            }
            double* mT_GAMBIT()
            {
                return new double(mT());
            }
            double* mT2_GAMBIT()
            {
                return new double(mT2());
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
            double* y_GAMBIT()
            {
                return new double(y());
            }
            double* eta_GAMBIT()
            {
                return new double(eta());
            }
            Pythia8::Vec4* vDec_GAMBIT()
            {
                return new Pythia8::Vec4(vDec());
            }
            double* xDec_GAMBIT()
            {
                return new double(xDec());
            }
            double* yDec_GAMBIT()
            {
                return new double(yDec());
            }
            double* zDec_GAMBIT()
            {
                return new double(zDec());
            }
            double* tDec_GAMBIT()
            {
                return new double(tDec());
            }
            int* index_GAMBIT()
            {
                return new int(index());
            }
            int* statusHepMC_GAMBIT()
            {
                return new int(statusHepMC());
            }
            int* iTopCopy_GAMBIT()
            {
                return new int(iTopCopy());
            }
            int* iBotCopy_GAMBIT()
            {
                return new int(iBotCopy());
            }
            int* iTopCopyId_GAMBIT()
            {
                return new int(iTopCopyId());
            }
            int* iBotCopyId_GAMBIT()
            {
                return new int(iBotCopyId());
            }
            bool* isAncestor_GAMBIT(int iAncestor)
            {
                return new bool(isAncestor(iAncestor));
            }
            bool* undoDecay_GAMBIT()
            {
                return new bool(undoDecay());
            }
            std::string* name_GAMBIT()
            {
                return new std::string(name());
            }
            std::string* nameWithStatus_GAMBIT(int maxLen)
            {
                return new std::string(nameWithStatus(maxLen));
            }
            int* spinType_GAMBIT()
            {
                return new int(spinType());
            }
            int* chargeType_GAMBIT()
            {
                return new int(chargeType());
            }
            double* charge_GAMBIT()
            {
                return new double(charge());
            }
            bool* isCharged_GAMBIT()
            {
                return new bool(isCharged());
            }
            bool* isNeutral_GAMBIT()
            {
                return new bool(isNeutral());
            }
            int* colType_GAMBIT()
            {
                return new int(colType());
            }
            double* m0_GAMBIT()
            {
                return new double(m0());
            }
            double* mWidth_GAMBIT()
            {
                return new double(mWidth());
            }
            double* mMin_GAMBIT()
            {
                return new double(mMin());
            }
            double* mMax_GAMBIT()
            {
                return new double(mMax());
            }
            double* mSel_GAMBIT()
            {
                return new double(mSel());
            }
            double* constituentMass_GAMBIT()
            {
                return new double(constituentMass());
            }
            double* tau0_GAMBIT()
            {
                return new double(tau0());
            }
            bool* mayDecay_GAMBIT()
            {
                return new bool(mayDecay());
            }
            bool* canDecay_GAMBIT()
            {
                return new bool(canDecay());
            }
            bool* doExternalDecay_GAMBIT()
            {
                return new bool(doExternalDecay());
            }
            bool* isResonance_GAMBIT()
            {
                return new bool(isResonance());
            }
            bool* isVisible_GAMBIT()
            {
                return new bool(isVisible());
            }
            bool* isLepton_GAMBIT()
            {
                return new bool(isLepton());
            }
            bool* isQuark_GAMBIT()
            {
                return new bool(isQuark());
            }
            bool* isGluon_GAMBIT()
            {
                return new bool(isGluon());
            }
            bool* isDiquark_GAMBIT()
            {
                return new bool(isDiquark());
            }
            bool* isParton_GAMBIT()
            {
                return new bool(isParton());
            }
            bool* isHadron_GAMBIT()
            {
                return new bool(isHadron());
            }
            void rescale3_GAMBIT(double fac)
            {
                rescale3(fac);
            }
            void rescale4_GAMBIT(double fac)
            {
                rescale4(fac);
            }
            void rescale5_GAMBIT(double fac)
            {
                rescale5(fac);
            }
            void rot_GAMBIT(double thetaIn, double phiIn)
            {
                rot(thetaIn,phiIn);
            }
            void bst_GAMBIT(double betaX, double betaY, double betaZ)
            {
                bst(betaX,betaY,betaZ);
            }
            void bst_GAMBIT(double betaX, double betaY, double betaZ, double gamma)
            {
                bst(betaX,betaY,betaZ,gamma);
            }
            void bst_GAMBIT(const Pythia8::Abstract__Vec4& pBst)
            {
                const Pythia8::Vec4& _temp_pBst = *(reinterpret_cast<const Pythia8::Vec4*>(&pBst));
                bst(_temp_pBst);
            }
            void bst_GAMBIT(const Pythia8::Abstract__Vec4& pBst, double mBst)
            {
                const Pythia8::Vec4& _temp_pBst = *(reinterpret_cast<const Pythia8::Vec4*>(&pBst));
                bst(_temp_pBst,mBst);
            }
            void bstback_GAMBIT(const Pythia8::Abstract__Vec4& pBst)
            {
                const Pythia8::Vec4& _temp_pBst = *(reinterpret_cast<const Pythia8::Vec4*>(&pBst));
                bstback(_temp_pBst);
            }
            void bstback_GAMBIT(const Pythia8::Abstract__Vec4& pBst, double mBst)
            {
                const Pythia8::Vec4& _temp_pBst = *(reinterpret_cast<const Pythia8::Vec4*>(&pBst));
                bstback(_temp_pBst,mBst);
            }
            void rotbst_GAMBIT(const Pythia8::Abstract__RotBstMatrix& M)
            {
                const Pythia8::RotBstMatrix& _temp_M = *(reinterpret_cast<const Pythia8::RotBstMatrix*>(&M));
                rotbst(_temp_M);
            }
            void offsetHistory_GAMBIT(int minMother, int addMother, int minDaughter, int addDaughter)
            {
                offsetHistory(minMother,addMother,minDaughter,addDaughter);
            }
            void offsetCol_GAMBIT(int addCol)
            {
                offsetCol(addCol);
            }
