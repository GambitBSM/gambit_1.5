#ifndef __ABSTRACT_PARTICLE_H__
#define __ABSTRACT_PARTICLE_H__

#include <iostream>  // FOR DEBUG: Allow virtual member functions to print a warning if executed.

#include "forward_decls_abstract_classes.h"
#include <vector>
#include "abstract_Vec4.h"
#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wreturn-type"
namespace Pythia8
{
    class Abstract_Particle
    {
        private:
            // IGNORED: Variable  -- Name: TINY  -- XML id: _19566
            // IGNORED: Field  -- Name: idSave  -- XML id: _19567
            // IGNORED: Field  -- Name: statusSave  -- XML id: _19568
            // IGNORED: Field  -- Name: mother1Save  -- XML id: _19569
            // IGNORED: Field  -- Name: mother2Save  -- XML id: _19570
            // IGNORED: Field  -- Name: daughter1Save  -- XML id: _19571
            // IGNORED: Field  -- Name: daughter2Save  -- XML id: _19572
            // IGNORED: Field  -- Name: colSave  -- XML id: _19573
            // IGNORED: Field  -- Name: acolSave  -- XML id: _19574
            // IGNORED: Field  -- Name: pSave  -- XML id: _19575
            // IGNORED: Field  -- Name: mSave  -- XML id: _19576
            // IGNORED: Field  -- Name: scaleSave  -- XML id: _19577
            // IGNORED: Field  -- Name: polSave  -- XML id: _19578
            // IGNORED: Field  -- Name: hasVertexSave  -- XML id: _19579
            // IGNORED: Field  -- Name: vProdSave  -- XML id: _19580
            // IGNORED: Field  -- Name: tauSave  -- XML id: _19581
            // IGNORED: Field  -- Name: pdePtr  -- XML id: _19582
            // IGNORED: Field  -- Name: evtPtr  -- XML id: _19583
        public:

            virtual Pythia8::Abstract_Particle* operator_assignment_gambit(const Pythia8::Abstract_Particle& pt) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Particle* operator=(const Pythia8::Abstract_Particle& pt)
            {
                return operator_assignment_gambit(pt);
            }

            virtual void setEvtPtr_gambit(Pythia8::Abstract_Event* evtPtrIn) {std::cout << "Called virtual function" << std::endl;};
            void setEvtPtr(Pythia8::Abstract_Event* evtPtrIn)
            {
                setEvtPtr_gambit(evtPtrIn);
            }

            virtual void id(int idIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void status(int statusIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void statusPos() {std::cout << "Called virtual function" << std::endl;};

            virtual void statusNeg() {std::cout << "Called virtual function" << std::endl;};

            virtual void statusCode(int statusIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void mother1(int mother1In) {std::cout << "Called virtual function" << std::endl;};

            virtual void mother2(int mother2In) {std::cout << "Called virtual function" << std::endl;};

            virtual void mothers(int mother1In, int mother2In) {std::cout << "Called virtual function" << std::endl;};

            virtual void daughter1(int daughter1In) {std::cout << "Called virtual function" << std::endl;};

            virtual void daughter2(int daughter2In) {std::cout << "Called virtual function" << std::endl;};

            virtual void daughters(int daughter1In, int daughter2In) {std::cout << "Called virtual function" << std::endl;};

            virtual void col(int colIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void acol(int acolIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void cols(int colIn, int acolIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void p_gambit(Pythia8::Abstract_Vec4& pIn) {std::cout << "Called virtual function" << std::endl;};
            void p(Pythia8::Abstract_Vec4& pIn)
            {
                p_gambit(pIn);
            }

            virtual void p(double pxIn, double pyIn, double pzIn, double eIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void px(double pxIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void py(double pyIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void pz(double pzIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void e(double eIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void m(double mIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void scale(double scaleIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void pol(double polIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void vProd_gambit(Pythia8::Abstract_Vec4& vProdIn) {std::cout << "Called virtual function" << std::endl;};
            void vProd(Pythia8::Abstract_Vec4& vProdIn)
            {
                vProd_gambit(vProdIn);
            }

            virtual void vProd(double xProdIn, double yProdIn, double zProdIn, double tProdIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void xProd(double xProdIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void yProd(double yProdIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void zProd(double zProdIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void tProd(double tProdIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void tau(double tauIn) {std::cout << "Called virtual function" << std::endl;};

            virtual int id() const {std::cout << "Called virtual function" << std::endl;};

            virtual int status() const {std::cout << "Called virtual function" << std::endl;};

            virtual int mother1() const {std::cout << "Called virtual function" << std::endl;};

            virtual int mother2() const {std::cout << "Called virtual function" << std::endl;};

            virtual int daughter1() const {std::cout << "Called virtual function" << std::endl;};

            virtual int daughter2() const {std::cout << "Called virtual function" << std::endl;};

            virtual int col() const {std::cout << "Called virtual function" << std::endl;};

            virtual int acol() const {std::cout << "Called virtual function" << std::endl;};

            virtual Pythia8::Abstract_Vec4* p_gambit() const {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Vec4* p() const
            {
                return p_gambit();
            }

            virtual double px() const {std::cout << "Called virtual function" << std::endl;};

            virtual double py() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pz() const {std::cout << "Called virtual function" << std::endl;};

            virtual double e() const {std::cout << "Called virtual function" << std::endl;};

            virtual double m() const {std::cout << "Called virtual function" << std::endl;};

            virtual double scale() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pol() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool hasVertex() const {std::cout << "Called virtual function" << std::endl;};

            virtual Pythia8::Abstract_Vec4* vProd_gambit() const {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Vec4* vProd() const
            {
                return vProd_gambit();
            }

            virtual double xProd() const {std::cout << "Called virtual function" << std::endl;};

            virtual double yProd() const {std::cout << "Called virtual function" << std::endl;};

            virtual double zProd() const {std::cout << "Called virtual function" << std::endl;};

            virtual double tProd() const {std::cout << "Called virtual function" << std::endl;};

            virtual double tau() const {std::cout << "Called virtual function" << std::endl;};

            virtual int idAbs() const {std::cout << "Called virtual function" << std::endl;};

            virtual int statusAbs() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isFinal() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isRescatteredIncoming() const {std::cout << "Called virtual function" << std::endl;};

            virtual double m2() const {std::cout << "Called virtual function" << std::endl;};

            virtual double mCalc() const {std::cout << "Called virtual function" << std::endl;};

            virtual double m2Calc() const {std::cout << "Called virtual function" << std::endl;};

            virtual double eCalc() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pT() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pT2() const {std::cout << "Called virtual function" << std::endl;};

            virtual double mT() const {std::cout << "Called virtual function" << std::endl;};

            virtual double mT2() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pAbs() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pAbs2() const {std::cout << "Called virtual function" << std::endl;};

            virtual double eT() const {std::cout << "Called virtual function" << std::endl;};

            virtual double eT2() const {std::cout << "Called virtual function" << std::endl;};

            virtual double theta() const {std::cout << "Called virtual function" << std::endl;};

            virtual double phi() const {std::cout << "Called virtual function" << std::endl;};

            virtual double thetaXZ() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pPos() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pNeg() const {std::cout << "Called virtual function" << std::endl;};

            virtual double y() const {std::cout << "Called virtual function" << std::endl;};

            virtual double eta() const {std::cout << "Called virtual function" << std::endl;};

            virtual Pythia8::Abstract_Vec4* vDec_gambit() const {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Vec4* vDec() const
            {
                return vDec_gambit();
            }

            virtual double xDec() const {std::cout << "Called virtual function" << std::endl;};

            virtual double yDec() const {std::cout << "Called virtual function" << std::endl;};

            virtual double zDec() const {std::cout << "Called virtual function" << std::endl;};

            virtual double tDec() const {std::cout << "Called virtual function" << std::endl;};

            virtual int index() const {std::cout << "Called virtual function" << std::endl;};

            virtual int statusHepMC() const {std::cout << "Called virtual function" << std::endl;};

            virtual int iTopCopy() const {std::cout << "Called virtual function" << std::endl;};

            virtual int iBotCopy() const {std::cout << "Called virtual function" << std::endl;};

            virtual int iTopCopyId() const {std::cout << "Called virtual function" << std::endl;};

            virtual int iBotCopyId() const {std::cout << "Called virtual function" << std::endl;};

            virtual std::vector<int, std::allocator<int> > motherList() const {std::cout << "Called virtual function" << std::endl;};

            virtual std::vector<int, std::allocator<int> > daughterList() const {std::cout << "Called virtual function" << std::endl;};

            virtual std::vector<int, std::allocator<int> > sisterList(bool traceTopBot) const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isAncestor(int iAncestor) const {std::cout << "Called virtual function" << std::endl;};

            virtual bool undoDecay() {std::cout << "Called virtual function" << std::endl;};

            virtual std::string name() const {std::cout << "Called virtual function" << std::endl;};

            virtual std::string nameWithStatus(int maxLen) const {std::cout << "Called virtual function" << std::endl;};

            virtual int spinType() const {std::cout << "Called virtual function" << std::endl;};

            virtual int chargeType() const {std::cout << "Called virtual function" << std::endl;};

            virtual double charge() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isCharged() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isNeutral() const {std::cout << "Called virtual function" << std::endl;};

            virtual int colType() const {std::cout << "Called virtual function" << std::endl;};

            virtual double m0() const {std::cout << "Called virtual function" << std::endl;};

            virtual double mWidth() const {std::cout << "Called virtual function" << std::endl;};

            virtual double mMin() const {std::cout << "Called virtual function" << std::endl;};

            virtual double mMax() const {std::cout << "Called virtual function" << std::endl;};

            virtual double mSel() const {std::cout << "Called virtual function" << std::endl;};

            virtual double constituentMass() const {std::cout << "Called virtual function" << std::endl;};

            virtual double tau0() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool mayDecay() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool canDecay() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool doExternalDecay() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isResonance() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isVisible() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isLepton() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isQuark() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isGluon() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isDiquark() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isParton() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isHadron() const {std::cout << "Called virtual function" << std::endl;};

            virtual void rescale3(double fac) {std::cout << "Called virtual function" << std::endl;};

            virtual void rescale4(double fac) {std::cout << "Called virtual function" << std::endl;};

            virtual void rescale5(double fac) {std::cout << "Called virtual function" << std::endl;};

            virtual void rot(double thetaIn, double phiIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void bst(double betaX, double betaY, double betaZ) {std::cout << "Called virtual function" << std::endl;};

            virtual void bst(double betaX, double betaY, double betaZ, double gamma) {std::cout << "Called virtual function" << std::endl;};

            virtual void bst_gambit(const Pythia8::Abstract_Vec4& pBst) {std::cout << "Called virtual function" << std::endl;};
            void bst(const Pythia8::Abstract_Vec4& pBst)
            {
                bst_gambit(pBst);
            }

            virtual void bst_gambit(const Pythia8::Abstract_Vec4& pBst, double mBst) {std::cout << "Called virtual function" << std::endl;};
            void bst(const Pythia8::Abstract_Vec4& pBst, double mBst)
            {
                bst_gambit(pBst, mBst);
            }

            virtual void bstback_gambit(const Pythia8::Abstract_Vec4& pBst) {std::cout << "Called virtual function" << std::endl;};
            void bstback(const Pythia8::Abstract_Vec4& pBst)
            {
                bstback_gambit(pBst);
            }

            virtual void bstback_gambit(const Pythia8::Abstract_Vec4& pBst, double mBst) {std::cout << "Called virtual function" << std::endl;};
            void bstback(const Pythia8::Abstract_Vec4& pBst, double mBst)
            {
                bstback_gambit(pBst, mBst);
            }

            virtual void offsetHistory(int minMother, int addMother, int minDaughter, int addDaughter) {std::cout << "Called virtual function" << std::endl;};

            virtual void offsetCol(int addCol) {std::cout << "Called virtual function" << std::endl;};

        public:
            virtual Abstract_Particle* pointerCopy_gambit() {std::cout << "Called virtual function" << std::endl;};
            virtual ~Abstract_Particle() {};
    };
}
#pragma GCC diagnostic pop


#endif /* __ABSTRACT_PARTICLE_H__ */
