// Event.h is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the Particle and Event classes.
// Particle: information on an instance of a particle.
// Junction: information on a junction between three colours.
// Event: list of particles in the current event.

#ifndef Pythia8_Event_H
#define Pythia8_Event_H

#include "Pythia8/Basics.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// Forward references to ParticleDataEntry and ResonanceWidths classes.
class ParticleDataEntry;
class ResonanceWidths;
class Event;

//==========================================================================

// Particle class.
// This class holds info on a particle in general.

} 
#include "Pythia8/abstract_Particle.h"
namespace Pythia8 { 
class Particle : public virtual Abstract__Particle {

public:

  // Constructors.
  Particle() : idSave(0), statusSave(0), mother1Save(0), mother2Save(0),
    daughter1Save(0), daughter2Save(0), colSave(0), acolSave(0),
    pSave(Vec4(0.,0.,0.,0.)), mSave(0.), scaleSave(0.), polSave(9.),
    hasVertexSave(false), vProdSave(Vec4(0.,0.,0.,0.)), tauSave(0.),
    pdePtr(0), evtPtr(0) { }
  Particle(int idIn, int statusIn = 0, int mother1In = 0,
    int mother2In = 0, int daughter1In = 0, int daughter2In = 0,
    int colIn = 0, int acolIn = 0, double pxIn = 0., double pyIn = 0.,
    double pzIn = 0., double eIn = 0., double mIn = 0.,
    double scaleIn = 0., double polIn = 9.)
    : idSave(idIn), statusSave(statusIn), mother1Save(mother1In),
    mother2Save(mother2In), daughter1Save(daughter1In),
    daughter2Save(daughter2In), colSave(colIn), acolSave(acolIn),
    pSave(Vec4(pxIn, pyIn, pzIn, eIn)), mSave(mIn), scaleSave(scaleIn),
    polSave(polIn), hasVertexSave(false), vProdSave(Vec4(0.,0.,0.,0.)),
    tauSave(0.), pdePtr(0), evtPtr(0) { }
  Particle(int idIn, int statusIn, int mother1In, int mother2In,
    int daughter1In, int daughter2In, int colIn, int acolIn,
    Vec4 pIn, double mIn = 0., double scaleIn = 0., double polIn = 9.)
    : idSave(idIn), statusSave(statusIn), mother1Save(mother1In),
    mother2Save(mother2In), daughter1Save(daughter1In),
    daughter2Save(daughter2In), colSave(colIn), acolSave(acolIn),
    pSave(pIn), mSave(mIn), scaleSave(scaleIn), polSave(polIn),
    hasVertexSave(false), vProdSave(Vec4(0.,0.,0.,0.)), tauSave(0.),
    pdePtr(0), evtPtr(0) { }
  Particle(const Particle& pt) : idSave(pt.idSave),
    statusSave(pt.statusSave), mother1Save(pt.mother1Save),
    mother2Save(pt.mother2Save), daughter1Save(pt.daughter1Save),
    daughter2Save(pt.daughter2Save), colSave(pt.colSave),
    acolSave(pt.acolSave), pSave(pt.pSave), mSave(pt.mSave),
    scaleSave(pt.scaleSave), polSave(pt.polSave),
    hasVertexSave(pt.hasVertexSave), vProdSave(pt.vProdSave),
    tauSave(pt.tauSave), pdePtr(pt.pdePtr), evtPtr(pt.evtPtr) { }
  Particle& operator=(const Particle& pt) {if (this != &pt) {
    idSave = pt.idSave; statusSave = pt.statusSave;
    mother1Save = pt.mother1Save; mother2Save = pt.mother2Save;
    daughter1Save = pt.daughter1Save; daughter2Save = pt.daughter2Save;
    colSave = pt.colSave; acolSave = pt.acolSave; pSave = pt.pSave;
    mSave = pt.mSave; scaleSave = pt.scaleSave; polSave = pt.polSave;
    hasVertexSave = pt.hasVertexSave; vProdSave = pt.vProdSave;
    tauSave = pt.tauSave; pdePtr = pt.pdePtr; evtPtr = pt.evtPtr; }
    return *this; }

  // Member functions to set the Event and ParticleDataEntry pointers.
  void setEvtPtr(Event* evtPtrIn) { evtPtr = evtPtrIn; setPDEPtr();}
  void setPDEPtr(ParticleDataEntry* pdePtrIn = 0);
      
  // Member functions for input.
  void id(int idIn) {idSave = idIn; setPDEPtr();}
  void status(int statusIn) {statusSave = statusIn;}
  void statusPos() {statusSave = abs(statusSave);}
  void statusNeg() {statusSave = -abs(statusSave);}
  void statusCode(int statusIn) {statusSave =
    (statusSave > 0) ? abs(statusIn) : -abs(statusIn);}
  void mother1(int mother1In) {mother1Save = mother1In;}
  void mother2(int mother2In) {mother2Save = mother2In;}
  void mothers(int mother1In = 0, int mother2In = 0)
    {mother1Save = mother1In; mother2Save = mother2In;}
  void daughter1(int daughter1In) {daughter1Save = daughter1In;}
  void daughter2(int daughter2In) {daughter2Save = daughter2In;}
  void daughters(int daughter1In = 0, int daughter2In = 0)
    {daughter1Save = daughter1In; daughter2Save = daughter2In;}
  void col(int colIn) {colSave = colIn;}
  void acol(int acolIn) {acolSave = acolIn;}
  void cols(int colIn = 0,int acolIn = 0) {colSave = colIn;
    acolSave = acolIn;}
  void p(Vec4 pIn) {pSave = pIn;}
  void p(double pxIn, double pyIn, double pzIn, double eIn)
    {pSave.p(pxIn, pyIn, pzIn, eIn);}
  void px(double pxIn) {pSave.px(pxIn);}
  void py(double pyIn) {pSave.py(pyIn);}
  void pz(double pzIn) {pSave.pz(pzIn);}
  void e(double eIn) {pSave.e(eIn);}
  void m(double mIn) {mSave = mIn;}
  void scale(double scaleIn) {scaleSave = scaleIn;}
  void pol(double polIn) {polSave = polIn;}
  void vProd(Vec4 vProdIn) {vProdSave = vProdIn; hasVertexSave = true;}
  void vProd(double xProdIn, double yProdIn, double zProdIn, double tProdIn)
    {vProdSave.p(xProdIn, yProdIn, zProdIn, tProdIn); hasVertexSave = true;}
  void xProd(double xProdIn) {vProdSave.px(xProdIn); hasVertexSave = true;}
  void yProd(double yProdIn) {vProdSave.py(yProdIn); hasVertexSave = true;}
  void zProd(double zProdIn) {vProdSave.pz(zProdIn); hasVertexSave = true;}
  void tProd(double tProdIn) {vProdSave.e(tProdIn); hasVertexSave = true;}
  void tau(double tauIn) {tauSave = tauIn;}

  // Member functions for output.
  int    id()        const {return idSave;}
  int    status()    const {return statusSave;}
  int    mother1()   const {return mother1Save;}
  int    mother2()   const {return mother2Save;}
  int    daughter1() const {return daughter1Save;}
  int    daughter2() const {return daughter2Save;}
  int    col()       const {return colSave;}
  int    acol()      const {return acolSave;}
  Vec4   p()         const {return pSave;}
  double px()        const {return pSave.px();}
  double py()        const {return pSave.py();}
  double pz()        const {return pSave.pz();}
  double e()         const {return pSave.e();}
  double m()         const {return mSave;}
  double scale()     const {return scaleSave;}
  double pol()       const {return polSave;}
  bool   hasVertex() const {return hasVertexSave;}
  Vec4   vProd()     const {return vProdSave;}
  double xProd()     const {return vProdSave.px();}
  double yProd()     const {return vProdSave.py();}
  double zProd()     const {return vProdSave.pz();}
  double tProd()     const {return vProdSave.e();}
  double tau()       const {return tauSave;}

  // Member functions for output; derived int and bool quantities.
  int    idAbs()     const {return abs(idSave);}
  int    statusAbs() const {return abs(statusSave);}
  bool   isFinal()   const {return (statusSave > 0);}
  bool   isRescatteredIncoming() const {return
    (statusSave == -34 || statusSave == -45 ||
     statusSave == -46 || statusSave == -54);}

  // Member functions for output; derived double quantities.
  double m2()        const {return (mSave >= 0.) ?  mSave*mSave
                                                 : -mSave*mSave;}
  double mCalc()     const {return pSave.mCalc();}
  double m2Calc()    const {return pSave.m2Calc();}
  double eCalc()     const {return sqrt(abs(m2() + pSave.pAbs2()));}
  double pT()        const {return pSave.pT();}
  double pT2()       const {return pSave.pT2();}
  double mT()        const {double temp = m2() + pSave.pT2();
    return (temp >= 0.) ? sqrt(temp) : -sqrt(-temp);}
  double mT2()       const {return m2() + pSave.pT2();}
  double pAbs()      const {return pSave.pAbs();}
  double pAbs2()     const {return pSave.pAbs2();}
  double eT()        const {return pSave.eT();}
  double eT2()       const {return pSave.eT2();}
  double theta()     const {return pSave.theta();}
  double phi()       const {return pSave.phi();}
  double thetaXZ()   const {return pSave.thetaXZ();}
  double pPos()      const {return pSave.pPos();}
  double pNeg()      const {return pSave.pNeg();}
  double y()         const;
  double eta()       const;
  Vec4   vDec()      const {return (tauSave > 0. && mSave > 0.)
    ? vProdSave + tauSave * pSave / mSave : vProdSave;}
  double xDec()      const {return (tauSave > 0. && mSave > 0.)
    ? vProdSave.px() + tauSave * pSave.px() / mSave : vProdSave.px();}
  double yDec()      const {return (tauSave > 0. && mSave > 0.)
    ? vProdSave.py() + tauSave * pSave.py() / mSave : vProdSave.py();}
  double zDec()      const {return (tauSave > 0. && mSave > 0.)
    ? vProdSave.pz() + tauSave * pSave.pz() / mSave : vProdSave.pz();}
  double tDec()      const {return (tauSave > 0. && mSave > 0.)
    ? vProdSave.e()  + tauSave * pSave.e()  / mSave : vProdSave.e();}

  // Methods that can refer back to the event the particle belongs to.
  int index()        const;
  int statusHepMC()  const;
  int iTopCopy()     const;
  int iBotCopy()     const;
  int iTopCopyId()   const;
  int iBotCopyId()   const;
  vector<int> motherList()   const;
  vector<int> daughterList() const;
  vector<int> sisterList(bool traceTopBot = false) const;
  bool isAncestor(int iAncestor) const;
  bool undoDecay();

  // Further output, based on a pointer to a ParticleDataEntry object.
  string name()      const {
    return (pdePtr != 0) ? pdePtr->name(idSave) : " ";}
  string nameWithStatus(int maxLen = 20) const;
  int    spinType()  const {
    return (pdePtr != 0) ? pdePtr->spinType() : 0;}
  int    chargeType() const {
    return (pdePtr != 0) ? pdePtr->chargeType(idSave) : 0;}
  double charge()    const {
    return (pdePtr != 0) ?  pdePtr->charge(idSave) : 0;}
  bool   isCharged() const {
    return (pdePtr != 0) ? (pdePtr->chargeType(idSave) != 0) : false;}
  bool   isNeutral() const {
    return (pdePtr != 0) ? (pdePtr->chargeType(idSave) == 0) : false;}
  int    colType()   const {
    return (pdePtr != 0) ? pdePtr->colType(idSave) : 0;}
  double m0()        const {
    return (pdePtr != 0) ? pdePtr->m0() : 0.;}
  double mWidth()    const {
    return (pdePtr != 0) ? pdePtr->mWidth() : 0.;}
  double mMin()      const {
    return (pdePtr != 0) ? pdePtr->mMin() : 0.;}
  double mMax()      const {
    return (pdePtr != 0) ? pdePtr->mMax() : 0.;}
  double mSel()      const {
    return (pdePtr != 0) ? pdePtr->mSel() : 0.;}
  double constituentMass() const {
    return (pdePtr != 0) ? pdePtr->constituentMass() : 0.;}
  double tau0()      const {
    return (pdePtr != 0) ? pdePtr->tau0() : 0.;}
  bool   mayDecay()  const {
    return (pdePtr != 0) ? pdePtr->mayDecay() : false;}
  bool   canDecay()  const {
    return (pdePtr != 0) ? pdePtr->canDecay() : false;}
  bool   doExternalDecay() const {
    return (pdePtr != 0) ? pdePtr->doExternalDecay() : false;}
  bool   isResonance() const {
    return (pdePtr != 0) ? pdePtr->isResonance() : false;}
  bool   isVisible() const {
    return (pdePtr != 0) ? pdePtr->isVisible() : false;}
  bool   isLepton()  const {
    return (pdePtr != 0) ? pdePtr->isLepton() : false;}
  bool   isQuark()   const {
    return (pdePtr != 0) ? pdePtr->isQuark() : false;}
  bool   isGluon()   const {
    return (pdePtr != 0) ? pdePtr->isGluon() : false;}
  bool   isDiquark()   const {
    return (pdePtr != 0) ? pdePtr->isDiquark() : false;}
  bool   isParton()   const {
    return (pdePtr != 0) ? pdePtr->isParton() : false;}
  bool   isHadron()  const {
    return (pdePtr != 0) ? pdePtr->isHadron() : false;}
  ParticleDataEntry& particleDataEntry() const {return *pdePtr;}

  // Member functions that perform operations.
  void rescale3(double fac) {pSave.rescale3(fac);}
  void rescale4(double fac) {pSave.rescale4(fac);}
  void rescale5(double fac) {pSave.rescale4(fac); mSave *= fac;}
  void rot(double thetaIn, double phiIn) {pSave.rot(thetaIn, phiIn);
    if (hasVertexSave) vProdSave.rot(thetaIn, phiIn);}
  void bst(double betaX, double betaY, double betaZ) {
    pSave.bst(betaX, betaY, betaZ);
    if (hasVertexSave) vProdSave.bst(betaX, betaY, betaZ);}
  void bst(double betaX, double betaY, double betaZ, double gamma) {
    pSave.bst(betaX, betaY, betaZ, gamma);
    if (hasVertexSave) vProdSave.bst(betaX, betaY, betaZ, gamma);}
  void bst(const Vec4& pBst) {pSave.bst(pBst);
    if (hasVertexSave) vProdSave.bst(pBst);}
  void bst(const Vec4& pBst, double mBst) {pSave.bst(pBst, mBst);
    if (hasVertexSave) vProdSave.bst(pBst, mBst);}
  void bstback(const Vec4& pBst) {pSave.bstback(pBst);
    if (hasVertexSave) vProdSave.bstback(pBst);}
  void bstback(const Vec4& pBst, double mBst) {pSave.bstback(pBst, mBst);
    if (hasVertexSave) vProdSave.bstback(pBst, mBst);}
  void rotbst(const RotBstMatrix& M) {pSave.rotbst(M);
    if (hasVertexSave) vProdSave.rotbst(M);}
  void offsetHistory( int minMother, int addMother, int minDaughter,
    int addDaughter);
  void offsetCol( int addCol);

private:

  // Constants: could only be changed in the code itself.
  static const double TINY;

  // Properties of the current particle.
  int    idSave, statusSave, mother1Save, mother2Save, daughter1Save,
         daughter2Save, colSave, acolSave;
  Vec4   pSave;
  double mSave, scaleSave, polSave;
  bool   hasVertexSave;
  Vec4   vProdSave;
  double tauSave;

  // Pointer to properties of the particle species.
  // Should no be saved in a persistent copy of the event record.
  // The //! below is ROOT notation that this member should not be saved.
  // Event::restorePtrs() can be called to restore the missing information.
  ParticleDataEntry* pdePtr;  //!

  // Pointer to the whole event record to which the particle belongs (if any).
  // As above it should not be saved.
  Event*             evtPtr;  //!


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
};

// Invariant mass of a pair and its square.
// (Not part of class proper, but tightly linked.)
double m(const Particle&, const Particle&);
double m2(const Particle&, const Particle&);

//==========================================================================

// The junction class stores what kind of junction it is, the colour indices
// of the legs at the junction and as far out as legs have been traced,
// and the status codes assigned for fragmentation of each leg.

class Junction {

public:

  // Constructors.
  Junction() : remainsSave(true), kindSave(0) {
    for (int j = 0; j < 3; ++j) {
    colSave[j] = 0; endColSave[j] = 0; statusSave[j] = 0; } }
  Junction( int kindIn, int col0In, int col1In, int col2In)
    : remainsSave(true), kindSave(kindIn) {colSave[0] = col0In;
    colSave[1] = col1In; colSave[2] = col2In;
    for (int j = 0; j < 3; ++j) {
    endColSave[j] = colSave[j]; statusSave[j] = 0; } }
  Junction(const Junction& ju) : remainsSave(ju.remainsSave),
    kindSave(ju.kindSave) { for (int j = 0; j < 3; ++j) {
    colSave[j] = ju.colSave[j]; endColSave[j] = ju.endColSave[j];
    statusSave[j] = ju.statusSave[j]; } }
  Junction& operator=(const Junction& ju) {if (this != &ju) {
    remainsSave = ju.remainsSave; kindSave =  ju.kindSave;
    for (int j = 0; j < 3; ++j) { colSave[j] = ju.colSave[j];
    endColSave[j] = ju.endColSave[j]; statusSave[j] = ju.statusSave[j]; } }
    return *this; }

  // Set values.
  void remains(bool remainsIn) {remainsSave = remainsIn;}
  void col(int j, int colIn) {colSave[j] = colIn; endColSave[j] = colIn;}
  void cols(int j, int colIn, int endColIn) {colSave[j] = colIn;
    endColSave[j] = endColIn;}
  void endCol(int j, int endColIn) {endColSave[j] = endColIn;}
  void status(int j, int statusIn) {statusSave[j] = statusIn;}

  // Read out value.
  bool   remains()     const {return remainsSave;}
  int    kind()        const {return kindSave;}
  int    col(int j)    const {return colSave[j];}
  int    endCol(int j) const {return endColSave[j];}
  int    status(int j) const {return statusSave[j];}
 
private:

  // Kind, positions of the three ends and their status codes.
  bool remainsSave;
  int kindSave, colSave[3], endColSave[3], statusSave[3];

};

//==========================================================================

// The Event class holds all info on the generated event.

} 
#include "Pythia8/abstract_Event.h"
namespace Pythia8 { 
class Event : public virtual Abstract__Event {
    
public:

  // Constructors.
  Event(int capacity = 100) : startColTag(100), maxColTag(100),
    savedSize(0), savedJunctionSize(0), scaleSave(0.), scaleSecondSave(0.),
    headerList("----------------------------------------"),
    particleDataPtr(0) { entry.reserve(capacity); }
  Event& operator=(const Event& oldEvent);
  Event(const Event& oldEvent) {*this = oldEvent;}

  // Initialize header for event listing, particle data table, and colour.
  void init( string headerIn = "", ParticleData* particleDataPtrIn = 0,
    int startColTagIn = 100) {
    headerList.replace(0, headerIn.length() + 2, headerIn + "  ");
     particleDataPtr = particleDataPtrIn; startColTag = startColTagIn;}

  // Clear event record.
  void clear() {entry.resize(0); maxColTag = startColTag; scaleSave = 0.;
    scaleSecondSave = 0.; clearJunctions();}

  // Clear event record, and set first particle empty.
  void reset() {clear(); append(90, -11, 0, 0, 0., 0., 0., 0., 0.);}

  // Overload index operator to access element of event record.
  Particle& operator[](int i) {return entry[i];}
  const Particle& operator[](int i) const {return entry[i];}

  // Implement standard references to elements in the particle array.
  Particle& front()   {return entry.front();}
  Particle& at(int i) {return entry.at(i);}
  Particle& back()    {return entry.back();}

  // Event record size.
  int size() const {return entry.size();}

  // Put a new particle at the end of the event record; return index.
  int append(Particle entryIn) {
    entry.push_back(entryIn); setEvtPtr();
    if (entryIn.col() > maxColTag) maxColTag = entryIn.col();
    if (entryIn.acol() > maxColTag) maxColTag = entryIn.acol();
    return entry.size() - 1;
  }
  int append(int id, int status, int mother1, int mother2, int daughter1,
    int daughter2, int col, int acol, double px, double py, double pz,
    double e, double m = 0., double scaleIn = 0., double polIn = 9.) {
    entry.push_back( Particle(id, status, mother1, mother2, daughter1,
    daughter2, col, acol, px, py, pz, e, m, scaleIn, polIn) ); setEvtPtr();
    if (col > maxColTag) maxColTag = col;
    if (acol > maxColTag) maxColTag = acol;
    return entry.size() - 1;
  }
  int append(int id, int status, int mother1, int mother2, int daughter1,
    int daughter2, int col, int acol, Vec4 p, double m = 0.,
    double scaleIn = 0., double polIn = 9.) {
    entry.push_back( Particle(id, status, mother1, mother2, daughter1,
    daughter2, col, acol, p, m, scaleIn, polIn) ); setEvtPtr();
    if (col > maxColTag) maxColTag = col;
    if (acol > maxColTag) maxColTag = acol;
    return entry.size() - 1;
  }

  // Brief versions of append: no mothers and no daughters.
  int append(int id, int status, int col, int acol, double px, double py,
    double pz, double e, double m = 0., double scaleIn = 0.,
    double polIn = 9.) { entry.push_back( Particle(id, status, 0, 0, 0, 0,
    col, acol, px, py, pz, e, m, scaleIn, polIn) ); setEvtPtr();
    if (col > maxColTag) maxColTag = col;
    if (acol > maxColTag) maxColTag = acol;
    return entry.size() - 1;
  }
  int append(int id, int status, int col, int acol, Vec4 p, double m = 0.,
    double scaleIn = 0., double polIn = 9.) {entry.push_back( Particle(id,
    status, 0, 0, 0, 0, col, acol, p, m, scaleIn, polIn) ); setEvtPtr();
    if (col > maxColTag) maxColTag = col;
    if (acol > maxColTag) maxColTag = acol;
    return entry.size() - 1;
  }

  // Set pointer to the event for a particle, by default latest one.
  void setEvtPtr(int iSet = -1) {if (iSet < 0) iSet = entry.size() - 1;
    entry[iSet].setEvtPtr( this);}

  // Add a copy of an existing particle at the end of the event record.
  int copy(int iCopy, int newStatus = 0);

  // List the particles in an event.
  void list() const;
  void list(ostream& os) const;
  void list(bool showScaleAndVertex, bool showMothersAndDaughters = false)
    const;
  void list(bool showScaleAndVertex, bool showMothersAndDaughters,
    ostream& os) const;

  // Remove last n entries.
  void popBack(int nRemove = 1) { if (nRemove ==1) entry.pop_back();
    else {int newSize = max( 0, size() - nRemove);
    entry.resize(newSize);} }

  // Remove entries from iFirst to iLast, including endpoints.
  void remove(int iFirst, int iLast) {
    if (iFirst < 0 || iLast >= int(entry.size()) || iLast < iFirst) return;
    entry.erase( entry.begin() + iFirst, entry.begin() + iLast + 1);
  }

  // Undo the decay of a single particle (where daughters well-defined).
  bool undoDecay(int i);

  // Restore all ParticleDataEntry* pointers in the Particle vector.
  // Useful when a persistent copy of the event record is read back in.
  void restorePtrs() { for (int i = 0; i < size(); ++i) setEvtPtr(i); }

  // Save or restore the size of the event record (throwing at the end).
  void saveSize() {savedSize = entry.size();}
  void restoreSize() {entry.resize(savedSize);}
  int  savedSizeValue() {return savedSize;}

  // Initialize and access colour tag information.
  void initColTag(int colTag = 0) {maxColTag = max( colTag,startColTag);}
  int lastColTag() const {return maxColTag;}
  int nextColTag() {return ++maxColTag;}

  // Access scale for which event as a whole is defined.
  void scale( double scaleIn) {scaleSave = scaleIn;}
  double scale() const {return scaleSave;}

  // Need a second scale if two hard interactions in event.
  void scaleSecond( double scaleSecondIn) {scaleSecondSave = scaleSecondIn;}
  double scaleSecond() const {return scaleSecondSave;}

  // Find complete list of daughters and mothers.
  vector<int> motherList(int i) const;
  vector<int> daughterList(int i) const;

  // Convert to HepMC status code conventions.
  int statusHepMC(int i) const;

  // Trace the first and last copy of one and the same particle.
  int iTopCopy(int i) const;
  int iBotCopy(int i) const;

  // Trace the first and last copy of a particle, using flavour match.
  int iTopCopyId(int i) const;
  int iBotCopyId(int i) const;

  // Find list of sisters, also tracking up and down identical copies.
  vector<int> sisterList(int i) const;
  vector<int> sisterListTopBot(int i, bool widenSearch = true) const;

  // Check whether two particles have a direct mother-daughter relation.
  bool isAncestor(int i, int iAncestor) const;

  // Member functions for rotations and boosts of an event.
  void rot(double theta, double phi)
    {for (int i = 0; i < size(); ++i) entry[i].rot(theta, phi);}
  void bst(double betaX, double betaY, double betaZ)
    {for (int i = 0; i < size(); ++i) entry[i].bst(betaX, betaY, betaZ);}
  void bst(double betaX, double betaY, double betaZ, double gamma)
    {for (int i = 0; i < size(); ++i) entry[i].bst(betaX, betaY, betaZ,
    gamma);}
  void bst(const Vec4& vec)
    {for (int i = 0; i < size(); ++i) entry[i].bst(vec);}
  void rotbst(const RotBstMatrix& M)
    {for (int i = 0; i < size(); ++i) entry[i].rotbst(M);}

  // Clear the list of junctions.
  void clearJunctions() {junction.resize(0);}
 
  // Add a junction to the list, study it or extra input.
  void appendJunction( int kind, int col0, int col1, int col2)
    { junction.push_back( Junction( kind, col0, col1, col2) );}
  void appendJunction(Junction junctionIn) {junction.push_back(junctionIn);}
  int sizeJunction() const {return junction.size();}
  bool remainsJunction(int i) const {return junction[i].remains();}
  void remainsJunction(int i, bool remainsIn) {junction[i].remains(remainsIn);}
  int kindJunction(int i) const {return junction[i].kind();}
  int colJunction( int i, int j) const {return junction[i].col(j);}
  void colJunction( int i, int j, int colIn) {junction[i].col(j, colIn);}
  int endColJunction( int i, int j) const {return junction[i].endCol(j);}
  void endColJunction( int i, int j, int endColIn)
    {junction[i].endCol(j, endColIn);}
  int statusJunction( int i, int j) const {return junction[i].status(j);}
  void statusJunction( int i, int j, int statusIn)
    {junction[i].status(j, statusIn);}
  Junction& getJunction(int i) {return junction[i];}
  const Junction& getJunction(int i) const {return junction[i];}
  void eraseJunction(int i);

  // Save or restore the size of the junction list (throwing at the end).
  void saveJunctionSize() {savedJunctionSize = junction.size();}
  void restoreJunctionSize() {junction.resize(savedJunctionSize);}

  // List any junctions in the event; for debug mainly.
  void listJunctions(ostream& os = cout) const;

  // Operator overloading allows to append one event to an existing one.
  // Warning: particles should be OK, but some other information unreliable.
  Event& operator+=(const Event& addEvent);

private:

  // The Particle class needs to access particle data.
  friend class Particle;

  // Constants: could only be changed in the code itself.
  static const int IPERLINE;

  // Initialization data, normally only set once.
  int startColTag;

  // The event: a vector containing all particles (entries).
  // The explicit use of Pythia8:: qualifier patches a limitation in ROOT.
  vector<Pythia8::Particle> entry;

  // The list of junctions.
  // The explicit use of Pythia8:: qualifier patches a limitation in ROOT.
  vector<Pythia8::Junction> junction;

  // The maximum colour tag of the event so far.
  int maxColTag;

  // Saved entry and junction list sizes, for simple restoration.
  int savedSize, savedJunctionSize;

  // The scale of the event; linear quantity in GeV.
  double scaleSave, scaleSecondSave;

  // Header specification in event listing (at most 40 characters wide).
  string headerList;

  // Pointer to the particle data table.
  // The //! below is ROOT notation that this member should not be saved.
  ParticleData* particleDataPtr;  //!
  

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
                return temp_pr;
            }
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
};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_Event_H
