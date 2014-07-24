#ifndef __GAMBIT_WRAPPER_EVENT_H__
#define __GAMBIT_WRAPPER_EVENT_H__

#include "GAMBIT_wrapper_WrapperBase.h"
#include <string>
#include "GAMBIT_wrapper_Particle.h"
#include "abstract_Event.h"
#include "GAMBIT_wrapper_Vec4.h"
#include <ostream>
#include <vector>


// Factory function pointers to be filled by dynamic loading
Pythia8::Abstract_Event* (*Factory_Event_0)(int) = NULL;

class Event_gambit : public WrapperBase<Pythia8::Abstract_Event>
{
    public:
        // Member variables: 

        // Member functions: 
        Event_gambit operator=(const WrapperBase< Pythia8::Abstract_Event > oldEvent)        {
            return BEptr->operator=(*oldEvent.BEptr);
        }

        void clear()        {
            BEptr->clear();
        }

        void reset()        {
            BEptr->reset();
        }

        Particle_gambit operator[](int i)        {
            return BEptr->operator[](i);
        }

        const Particle_gambit operator[](int i) const
        {
            return BEptr->operator[](i);
        }

        Particle_gambit front()        {
            return BEptr->front();
        }

        Particle_gambit at(int i)        {
            return BEptr->at(i);
        }

        Particle_gambit back()        {
            return BEptr->back();
        }

        const int size() const
        {
            return BEptr->size();
        }

        int append(WrapperBase< Pythia8::Abstract_Particle > entryIn)        {
            return BEptr->append(*entryIn.BEptr);
        }

        int append(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn, double polIn)        {
            return BEptr->append(id, status, mother1, mother2, daughter1, daughter2, col, acol, px, py, pz, e, m, scaleIn, polIn);
        }

        int append(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, WrapperBase< Pythia8::Abstract_Vec4 > p, double m, double scaleIn, double polIn)        {
            return BEptr->append(id, status, mother1, mother2, daughter1, daughter2, col, acol, *p.BEptr, m, scaleIn, polIn);
        }

        int append(int id, int status, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn, double polIn)        {
            return BEptr->append(id, status, col, acol, px, py, pz, e, m, scaleIn, polIn);
        }

        int append(int id, int status, int col, int acol, WrapperBase< Pythia8::Abstract_Vec4 > p, double m, double scaleIn, double polIn)        {
            return BEptr->append(id, status, col, acol, *p.BEptr, m, scaleIn, polIn);
        }

        void setEvtPtr(int iSet)        {
            BEptr->setEvtPtr(iSet);
        }

        int copy(int iCopy, int newStatus)        {
            return BEptr->copy(iCopy, newStatus);
        }

        const void list() const
        {
            BEptr->list();
        }

        const void list(std::ostream& os) const
        {
            BEptr->list(os);
        }

        const void list(bool showScaleAndVertex, bool showMothersAndDaughters) const
        {
            BEptr->list(showScaleAndVertex, showMothersAndDaughters);
        }

        const void list(bool showScaleAndVertex, bool showMothersAndDaughters, std::ostream& os) const
        {
            BEptr->list(showScaleAndVertex, showMothersAndDaughters, os);
        }

        void popBack(int nRemove)        {
            BEptr->popBack(nRemove);
        }

        void remove(int iFirst, int iLast)        {
            BEptr->remove(iFirst, iLast);
        }

        bool undoDecay(int i)        {
            return BEptr->undoDecay(i);
        }

        void restorePtrs()        {
            BEptr->restorePtrs();
        }

        void saveSize()        {
            BEptr->saveSize();
        }

        void restoreSize()        {
            BEptr->restoreSize();
        }

        int savedSizeValue()        {
            return BEptr->savedSizeValue();
        }

        void initColTag(int colTag)        {
            BEptr->initColTag(colTag);
        }

        const int lastColTag() const
        {
            return BEptr->lastColTag();
        }

        int nextColTag()        {
            return BEptr->nextColTag();
        }

        void scale(double scaleIn)        {
            BEptr->scale(scaleIn);
        }

        const double scale() const
        {
            return BEptr->scale();
        }

        void scaleSecond(double scaleSecondIn)        {
            BEptr->scaleSecond(scaleSecondIn);
        }

        const double scaleSecond() const
        {
            return BEptr->scaleSecond();
        }

        const std::vector<int, std::allocator<int> > motherList(int i) const
        {
            return BEptr->motherList(i);
        }

        const std::vector<int, std::allocator<int> > daughterList(int i) const
        {
            return BEptr->daughterList(i);
        }

        const int statusHepMC(int i) const
        {
            return BEptr->statusHepMC(i);
        }

        const int iTopCopy(int i) const
        {
            return BEptr->iTopCopy(i);
        }

        const int iBotCopy(int i) const
        {
            return BEptr->iBotCopy(i);
        }

        const int iTopCopyId(int i) const
        {
            return BEptr->iTopCopyId(i);
        }

        const int iBotCopyId(int i) const
        {
            return BEptr->iBotCopyId(i);
        }

        const std::vector<int, std::allocator<int> > sisterList(int i) const
        {
            return BEptr->sisterList(i);
        }

        const std::vector<int, std::allocator<int> > sisterListTopBot(int i, bool widenSearch) const
        {
            return BEptr->sisterListTopBot(i, widenSearch);
        }

        const bool isAncestor(int i, int iAncestor) const
        {
            return BEptr->isAncestor(i, iAncestor);
        }

        void rot(double theta, double phi)        {
            BEptr->rot(theta, phi);
        }

        void bst(double betaX, double betaY, double betaZ)        {
            BEptr->bst(betaX, betaY, betaZ);
        }

        void bst(double betaX, double betaY, double betaZ, double gamma)        {
            BEptr->bst(betaX, betaY, betaZ, gamma);
        }

        void bst(const WrapperBase< Pythia8::Abstract_Vec4 > vec)        {
            BEptr->bst(*vec.BEptr);
        }

        void clearJunctions()        {
            BEptr->clearJunctions();
        }

        int appendJunction(int kind, int col0, int col1, int col2)        {
            return BEptr->appendJunction(kind, col0, col1, col2);
        }

        const int sizeJunction() const
        {
            return BEptr->sizeJunction();
        }

        const bool remainsJunction(int i) const
        {
            return BEptr->remainsJunction(i);
        }

        void remainsJunction(int i, bool remainsIn)        {
            BEptr->remainsJunction(i, remainsIn);
        }

        const int kindJunction(int i) const
        {
            return BEptr->kindJunction(i);
        }

        const int colJunction(int i, int j) const
        {
            return BEptr->colJunction(i, j);
        }

        void colJunction(int i, int j, int colIn)        {
            BEptr->colJunction(i, j, colIn);
        }

        const int endColJunction(int i, int j) const
        {
            return BEptr->endColJunction(i, j);
        }

        void endColJunction(int i, int j, int endColIn)        {
            BEptr->endColJunction(i, j, endColIn);
        }

        const int statusJunction(int i, int j) const
        {
            return BEptr->statusJunction(i, j);
        }

        void statusJunction(int i, int j, int statusIn)        {
            BEptr->statusJunction(i, j, statusIn);
        }

        void eraseJunction(int i)        {
            BEptr->eraseJunction(i);
        }

        void saveJunctionSize()        {
            BEptr->saveJunctionSize();
        }

        void restoreJunctionSize()        {
            BEptr->restoreJunctionSize();
        }

        const void listJunctions(std::ostream& os) const
        {
            BEptr->listJunctions(os);
        }

        Event_gambit operator+=(const WrapperBase< Pythia8::Abstract_Event > addEvent)        {
            return BEptr->operator+=(*addEvent.BEptr);
        }


        // Wrappers for original constructors: 
        Event_gambit(int capacity) :
            WrapperBase<Pythia8::Abstract_Event>( Factory_Event_0(capacity), false )
        {
        }

        // Special pointer-based constructor: 
        Event_gambit(Pythia8::Abstract_Event* in) :
            WrapperBase<Pythia8::Abstract_Event>( in, false )
        {
        }

        // Copy constructor: 
        Event_gambit(const Event_gambit& in) :
            WrapperBase<Pythia8::Abstract_Event>( in.BEptr->pointerCopy_gambit(), false )
        {
        }
};

#endif /* __GAMBIT_WRAPPER_EVENT_H__ */
