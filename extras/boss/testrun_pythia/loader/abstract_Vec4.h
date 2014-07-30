#ifndef __ABSTRACT_VEC4_H__
#define __ABSTRACT_VEC4_H__

#include <iostream>  // FOR DEBUG: Allow virtual member functions to print a warning if executed.

#include "forward_decls_abstract_classes.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wreturn-type"
namespace Pythia8
{
    class Abstract_Vec4
    {
        private:
            // IGNORED: Variable  -- Name: TINY  -- XML id: _21851
            // IGNORED: Field  -- Name: xx  -- XML id: _21852
            // IGNORED: Field  -- Name: yy  -- XML id: _21853
            // IGNORED: Field  -- Name: zz  -- XML id: _21854
            // IGNORED: Field  -- Name: tt  -- XML id: _21855
        public:

            virtual Pythia8::Abstract_Vec4* operator_assignment_gambit(const Pythia8::Abstract_Vec4& v) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Vec4* operator=(const Pythia8::Abstract_Vec4& v)
            {
                return operator_assignment_gambit(v);
            }

            virtual Pythia8::Abstract_Vec4* operator_assignment_gambit(double value) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Vec4* operator=(double value)
            {
                return operator_assignment_gambit(value);
            }

            virtual void reset() {std::cout << "Called virtual function" << std::endl;};

            virtual void p(double xIn, double yIn, double zIn, double tIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void p_gambit(Pythia8::Abstract_Vec4& pIn) {std::cout << "Called virtual function" << std::endl;};
            void p(Pythia8::Abstract_Vec4& pIn)
            {
                p_gambit(pIn);
            }

            virtual void px(double xIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void py(double yIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void pz(double zIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void e(double tIn) {std::cout << "Called virtual function" << std::endl;};

            virtual double px() const {std::cout << "Called virtual function" << std::endl;};

            virtual double py() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pz() const {std::cout << "Called virtual function" << std::endl;};

            virtual double e() const {std::cout << "Called virtual function" << std::endl;};

            virtual double mCalc() const {std::cout << "Called virtual function" << std::endl;};

            virtual double m2Calc() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pT() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pT2() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pAbs() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pAbs2() const {std::cout << "Called virtual function" << std::endl;};

            virtual double eT() const {std::cout << "Called virtual function" << std::endl;};

            virtual double eT2() const {std::cout << "Called virtual function" << std::endl;};

            virtual double theta() const {std::cout << "Called virtual function" << std::endl;};

            virtual double phi() const {std::cout << "Called virtual function" << std::endl;};

            virtual double thetaXZ() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pPos() const {std::cout << "Called virtual function" << std::endl;};

            virtual double pNeg() const {std::cout << "Called virtual function" << std::endl;};

            virtual double rap() const {std::cout << "Called virtual function" << std::endl;};

            virtual double eta() const {std::cout << "Called virtual function" << std::endl;};

            virtual void rescale3(double fac) {std::cout << "Called virtual function" << std::endl;};

            virtual void rescale4(double fac) {std::cout << "Called virtual function" << std::endl;};

            virtual void flip3() {std::cout << "Called virtual function" << std::endl;};

            virtual void flip4() {std::cout << "Called virtual function" << std::endl;};

            virtual void rot(double thetaIn, double phiIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void rotaxis(double phiIn, double nx, double ny, double nz) {std::cout << "Called virtual function" << std::endl;};

            virtual void rotaxis_gambit(double phiIn, const Pythia8::Abstract_Vec4& n) {std::cout << "Called virtual function" << std::endl;};
            void rotaxis(double phiIn, const Pythia8::Abstract_Vec4& n)
            {
                rotaxis_gambit(phiIn, n);
            }

            virtual void bst(double betaX, double betaY, double betaZ) {std::cout << "Called virtual function" << std::endl;};

            virtual void bst(double betaX, double betaY, double betaZ, double gamma) {std::cout << "Called virtual function" << std::endl;};

            virtual void bst_gambit(const Pythia8::Abstract_Vec4& pIn) {std::cout << "Called virtual function" << std::endl;};
            void bst(const Pythia8::Abstract_Vec4& pIn)
            {
                bst_gambit(pIn);
            }

            virtual void bst_gambit(const Pythia8::Abstract_Vec4& pIn, double mIn) {std::cout << "Called virtual function" << std::endl;};
            void bst(const Pythia8::Abstract_Vec4& pIn, double mIn)
            {
                bst_gambit(pIn, mIn);
            }

            virtual void bstback_gambit(const Pythia8::Abstract_Vec4& pIn) {std::cout << "Called virtual function" << std::endl;};
            void bstback(const Pythia8::Abstract_Vec4& pIn)
            {
                bstback_gambit(pIn);
            }

            virtual void bstback_gambit(const Pythia8::Abstract_Vec4& pIn, double mIn) {std::cout << "Called virtual function" << std::endl;};
            void bstback(const Pythia8::Abstract_Vec4& pIn, double mIn)
            {
                bstback_gambit(pIn, mIn);
            }

            virtual Pythia8::Abstract_Vec4* operator_subtraction_gambit() {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Vec4* operator-()
            {
                return operator_subtraction_gambit();
            }

            virtual Pythia8::Abstract_Vec4* operator_addition_assignment_gambit(const Pythia8::Abstract_Vec4& v) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Vec4* operator+=(const Pythia8::Abstract_Vec4& v)
            {
                return operator_addition_assignment_gambit(v);
            }

            virtual Pythia8::Abstract_Vec4* operator_subtraction_assignment_gambit(const Pythia8::Abstract_Vec4& v) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Vec4* operator-=(const Pythia8::Abstract_Vec4& v)
            {
                return operator_subtraction_assignment_gambit(v);
            }

            virtual Pythia8::Abstract_Vec4* operator_multiplication_assignment_gambit(double f) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Vec4* operator*=(double f)
            {
                return operator_multiplication_assignment_gambit(f);
            }

            virtual Pythia8::Abstract_Vec4* operator_division_assignment_gambit(double f) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Vec4* operator/=(double f)
            {
                return operator_division_assignment_gambit(f);
            }

        public:
            virtual Abstract_Vec4* pointerCopy_gambit() {std::cout << "Called virtual function" << std::endl;};
            virtual ~Abstract_Vec4() {};
    };
}
#pragma GCC diagnostic pop


#endif /* __ABSTRACT_VEC4_H__ */
