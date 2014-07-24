#ifndef __ABSTRACT_HIST_H__
#define __ABSTRACT_HIST_H__

#include <iostream>  // FOR DEBUG: Allow virtual member functions to print a warning if executed.

#include "forward_decls_abstract_classes.h"
#include <vector>
#include <ostream>
#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wreturn-type"
namespace Pythia8
{
    class Abstract_Hist
    {
        private:
            // IGNORED: Variable  -- Name: NBINMAX  -- XML id: _23079
            // IGNORED: Variable  -- Name: NCOLMAX  -- XML id: _23080
            // IGNORED: Variable  -- Name: NLINES  -- XML id: _23081
            // IGNORED: Variable  -- Name: TOLERANCE  -- XML id: _23082
            // IGNORED: Variable  -- Name: TINY  -- XML id: _23083
            // IGNORED: Variable  -- Name: LARGE  -- XML id: _23084
            // IGNORED: Variable  -- Name: SMALLFRAC  -- XML id: _23085
            // IGNORED: Variable  -- Name: DYAC  -- XML id: _23086
            // IGNORED: Variable  -- Name: NUMBER  -- XML id: _23087
            // IGNORED: Field  -- Name: title  -- XML id: _23088
            // IGNORED: Field  -- Name: nBin  -- XML id: _23089
            // IGNORED: Field  -- Name: nFill  -- XML id: _23090
            // IGNORED: Field  -- Name: xMin  -- XML id: _23091
            // IGNORED: Field  -- Name: xMax  -- XML id: _23092
            // IGNORED: Field  -- Name: dx  -- XML id: _23093
            // IGNORED: Field  -- Name: under  -- XML id: _23094
            // IGNORED: Field  -- Name: inside  -- XML id: _23095
            // IGNORED: Field  -- Name: over  -- XML id: _23096
            // IGNORED: Field  -- Name: res  -- XML id: _23097
        public:

            virtual Pythia8::Abstract_Hist* operator_assignment_gambit(const Pythia8::Abstract_Hist& h) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Hist* operator=(const Pythia8::Abstract_Hist& h)
            {
                return operator_assignment_gambit(h);
            }

            virtual void book(std::string titleIn, int nBinIn, double xMinIn, double xMaxIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void name(std::string titleIn) {std::cout << "Called virtual function" << std::endl;};

            virtual void null() {std::cout << "Called virtual function" << std::endl;};

            virtual void fill(double x, double w) {std::cout << "Called virtual function" << std::endl;};

            virtual void table(std::ostream& os, bool printOverUnder, bool xMidBin) const {std::cout << "Called virtual function" << std::endl;};

            virtual void table(std::string fileName, bool printOverUnder, bool xMidBin) const {std::cout << "Called virtual function" << std::endl;};

            virtual double getBinContent(int iBin) const {std::cout << "Called virtual function" << std::endl;};

            virtual int getEntries() const {std::cout << "Called virtual function" << std::endl;};

            virtual bool sameSize_gambit(const Pythia8::Abstract_Hist& h) const {std::cout << "Called virtual function" << std::endl;};
            bool sameSize(const Pythia8::Abstract_Hist& h) const
            {
                return sameSize_gambit(h);
            }

            virtual void takeLog(bool tenLog) {std::cout << "Called virtual function" << std::endl;};

            virtual void takeSqrt() {std::cout << "Called virtual function" << std::endl;};

            virtual Pythia8::Abstract_Hist* operator_addition_assignment_gambit(const Pythia8::Abstract_Hist& h) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Hist* operator+=(const Pythia8::Abstract_Hist& h)
            {
                return operator_addition_assignment_gambit(h);
            }

            virtual Pythia8::Abstract_Hist* operator_subtraction_assignment_gambit(const Pythia8::Abstract_Hist& h) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Hist* operator-=(const Pythia8::Abstract_Hist& h)
            {
                return operator_subtraction_assignment_gambit(h);
            }

            virtual Pythia8::Abstract_Hist* operator_multiplication_assignment_gambit(const Pythia8::Abstract_Hist& h) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Hist* operator*=(const Pythia8::Abstract_Hist& h)
            {
                return operator_multiplication_assignment_gambit(h);
            }

            virtual Pythia8::Abstract_Hist* operator_division_assignment_gambit(const Pythia8::Abstract_Hist& h) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Hist* operator/=(const Pythia8::Abstract_Hist& h)
            {
                return operator_division_assignment_gambit(h);
            }

            virtual Pythia8::Abstract_Hist* operator_addition_assignment_gambit(double f) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Hist* operator+=(double f)
            {
                return operator_addition_assignment_gambit(f);
            }

            virtual Pythia8::Abstract_Hist* operator_subtraction_assignment_gambit(double f) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Hist* operator-=(double f)
            {
                return operator_subtraction_assignment_gambit(f);
            }

            virtual Pythia8::Abstract_Hist* operator_multiplication_assignment_gambit(double f) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Hist* operator*=(double f)
            {
                return operator_multiplication_assignment_gambit(f);
            }

            virtual Pythia8::Abstract_Hist* operator_division_assignment_gambit(double f) {std::cout << "Called virtual function" << std::endl;};
            Pythia8::Abstract_Hist* operator/=(double f)
            {
                return operator_division_assignment_gambit(f);
            }

        public:
            virtual Abstract_Hist* pointerCopy_gambit() {std::cout << "Called virtual function" << std::endl;};
            virtual ~Abstract_Hist() {};
    };
}
#pragma GCC diagnostic pop


#endif /* __ABSTRACT_HIST_H__ */
