#ifndef __ABSTRACT_EVENT_H__
#define __ABSTRACT_EVENT_H__

#include <iostream>  // FOR DEBUG: Allow virtual member functions to print a warning if executed.

#include "forward_decls_abstract_classes.h"
#include "abstract_Particle.h"
#include <ostream>
#include <vector>
#include "abstract_Vec4.h"
#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wreturn-type"
namespace Pythia8
{
    class Abstract_Event
    {
        private:
            // IGNORED: Variable  -- Name: IPERLINE  -- XML id: _24569
            // IGNORED: Field  -- Name: startColTag  -- XML id: _24570
            // IGNORED: Field  -- Name: entry  -- XML id: _24571
            // IGNORED: Field  -- Name: junction  -- XML id: _24572
            // IGNORED: Field  -- Name: maxColTag  -- XML id: _24573
            // IGNORED: Field  -- Name: savedSize  -- XML id: _24574
            // IGNORED: Field  -- Name: savedJunctionSize  -- XML id: _24575
            // IGNORED: Field  -- Name: scaleSave  -- XML id: _24576
            // IGNORED: Field  -- Name: scaleSecondSave  -- XML id: _24577
            // IGNORED: Field  -- Name: headerList  -- XML id: _24578
            // IGNORED: Field  -- Name: particleDataPtr  -- XML id: _24579
        public:

            virtual Pythia8::Abstract_Event* operator_assignment_gambit(const Pythia8::Abstract_Event& oldEvent) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Event* operator=(const Pythia8::Abstract_Event& oldEvent)
            {
                return operator_assignment_gambit(oldEvent);
            }

            virtual void clear() {std::cout << "Called virtual function" << std::endl;};

            virtual void reset() {std::cout << "Called virtual function" << std::endl;};

            virtual Pythia8::Abstract_Particle* operator_array_subscript_gambit(int i) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Particle* operator[](int i)
            {
                return operator_array_subscript_gambit(i);
            }

            virtual const Pythia8::Abstract_Particle* operator_array_subscript_gambit(int i) const {std::cout << "Called virtual function" << std::endl;};
            const Pythia8::Abstract_Particle* operator[](int i) const
            {
                return operator_array_subscript_gambit(i);
            }

            virtual Pythia8::Abstract_Particle* front_gambit() {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Particle* front()
            {
                return front_gambit();
            }

            virtual Pythia8::Abstract_Particle* at_gambit(int i) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Particle* at(int i)
            {
                return at_gambit(i);
            }

            virtual Pythia8::Abstract_Particle* back_gambit() {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Particle* back()
            {
                return back_gambit();
            }

            virtual int size() const {std::cout << "Called virtual function" << std::endl;};

            virtual int append_gambit(Pythia8::Abstract_Particle& entryIn) {std::cout << "Called virtual function" << std::endl;};
            int append(Pythia8::Abstract_Particle& entryIn)
            {
                return append_gambit(entryIn);
            }

            virtual int append(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn, double polIn) {std::cout << "Called virtual function" << std::endl;};

            virtual int append_gambit(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract_Vec4& p, double m, double scaleIn, double polIn) {std::cout << "Called virtual function" << std::endl;};
            int append(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract_Vec4& p, double m, double scaleIn, double polIn)
            {
                return append_gambit(id, status, mother1, mother2, daughter1, daughter2, col, acol, p, m, scaleIn, polIn);
            }

            virtual int append(int id, int status, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn, double polIn) {std::cout << "Called virtual function" << std::endl;};

            virtual int append_gambit(int id, int status, int col, int acol, Pythia8::Abstract_Vec4& p, double m, double scaleIn, double polIn) {std::cout << "Called virtual function" << std::endl;};
            int append(int id, int status, int col, int acol, Pythia8::Abstract_Vec4& p, double m, double scaleIn, double polIn)
            {
                return append_gambit(id, status, col, acol, p, m, scaleIn, polIn);
            }

            virtual void setEvtPtr(int iSet) {std::cout << "Called virtual function" << std::endl;};

            virtual int copy(int iCopy, int newStatus) {std::cout << "Called virtual function" << std::endl;};

            virtual void list() const {std::cout << "Called virtual function" << std::endl;};

            virtual void list(std::ostream& os) const {std::cout << "Called virtual function" << std::endl;};

            virtual void list(bool showScaleAndVertex, bool showMothersAndDaughters) const {std::cout << "Called virtual function" << std::endl;};

            virtual void list(bool showScaleAndVertex, bool showMothersAndDaughters, std::ostream& os) const {std::cout << "Called virtual function" << std::endl;};

            virtual void popBack(int nRemove) {std::cout << "Called virtual function" << std::endl;};

            virtual void remove(int iFirst, int iLast) {std::cout << "Called virtual function" << std::endl;};

            virtual bool undoDecay(int i) {std::cout << "Called virtual function" << std::endl;};

            virtual void restorePtrs() {std::cout << "Called virtual function" << std::endl;};

            virtual void saveSize() {std::cout << "Called virtual function" << std::endl;};

            virtual void restoreSize() {std::cout << "Called virtual function" << std::endl;};

            virtual int savedSizeValue() {std::cout << "Called virtual function" << std::endl;};

            virtual void initColTag(int colTag) {std::cout << "Called virtual function" << std::endl;};

            virtual int lastColTag() const {std::cout << "Called virtual function" << std::endl;};

            virtual int nextColTag() {std::cout << "Called virtual function" << std::endl;};

            virtual void scale(double scaleIn) {std::cout << "Called virtual function" << std::endl;};

            virtual double scale() const {std::cout << "Called virtual function" << std::endl;};

            virtual void scaleSecond(double scaleSecondIn) {std::cout << "Called virtual function" << std::endl;};

            virtual double scaleSecond() const {std::cout << "Called virtual function" << std::endl;};

            virtual std::vector<int, std::allocator<int> > motherList(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual std::vector<int, std::allocator<int> > daughterList(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int statusHepMC(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int iTopCopy(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int iBotCopy(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int iTopCopyId(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int iBotCopyId(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual std::vector<int, std::allocator<int> > sisterList(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual std::vector<int, std::allocator<int> > sisterListTopBot(int i, bool widenSearch) const {std::cout << "Called virtual function" << std::endl;};

            virtual bool isAncestor(int i, int iAncestor) const {std::cout << "Called virtual function" << std::endl;};

            virtual void rot(double theta, double phi) {std::cout << "Called virtual function" << std::endl;};

            virtual void bst(double betaX, double betaY, double betaZ) {std::cout << "Called virtual function" << std::endl;};

            virtual void bst(double betaX, double betaY, double betaZ, double gamma) {std::cout << "Called virtual function" << std::endl;};

            virtual void bst_gambit(const Pythia8::Abstract_Vec4& vec) {std::cout << "Called virtual function" << std::endl;};
            void bst(const Pythia8::Abstract_Vec4& vec)
            {
                bst_gambit(vec);
            }

            virtual void clearJunctions() {std::cout << "Called virtual function" << std::endl;};

            virtual int appendJunction(int kind, int col0, int col1, int col2) {std::cout << "Called virtual function" << std::endl;};

            virtual int sizeJunction() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool remainsJunction(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual void remainsJunction(int i, bool remainsIn) {std::cout << "Called virtual function" << std::endl;};

            virtual int kindJunction(int i) const {std::cout << "Called virtual function" << std::endl;};

            virtual int colJunction(int i, int j) const {std::cout << "Called virtual function" << std::endl;};

            virtual void colJunction(int i, int j, int colIn) {std::cout << "Called virtual function" << std::endl;};

            virtual int endColJunction(int i, int j) const {std::cout << "Called virtual function" << std::endl;};

            virtual void endColJunction(int i, int j, int endColIn) {std::cout << "Called virtual function" << std::endl;};

            virtual int statusJunction(int i, int j) const {std::cout << "Called virtual function" << std::endl;};

            virtual void statusJunction(int i, int j, int statusIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void eraseJunction(int i) {std::cout << "Called virtual function" << std::endl;};

            virtual void saveJunctionSize() {std::cout << "Called virtual function" << std::endl;};

            virtual void restoreJunctionSize() {std::cout << "Called virtual function" << std::endl;};

            virtual void listJunctions(std::ostream& os) const {std::cout << "Called virtual function" << std::endl;};

            virtual Pythia8::Abstract_Event* operator_addition_assignment_gambit(const Pythia8::Abstract_Event& addEvent) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Event* operator+=(const Pythia8::Abstract_Event& addEvent)
            {
                return operator_addition_assignment_gambit(addEvent);
            }

        public:
            virtual Abstract_Event* pointerCopy_gambit() {std::cout << "Called virtual function" << std::endl;};
            virtual ~Abstract_Event() {};
    };
}
#pragma GCC diagnostic pop


#endif /* __ABSTRACT_EVENT_H__ */
