#include <string>
#include <vector>
#include <ostream>
#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"
#include "Pythia8/Basics.h"

void Pythia8::Hist::book_GAMBIT(std::string titleIn, int nBinIn, double xMinIn)
{
    book(titleIn, nBinIn, xMinIn);
}


void Pythia8::Hist::book_GAMBIT(std::string titleIn, int nBinIn)
{
    book(titleIn, nBinIn);
}


void Pythia8::Hist::book_GAMBIT(std::string titleIn)
{
    book(titleIn);
}


void Pythia8::Hist::book_GAMBIT()
{
    book();
}


void Pythia8::Hist::name_GAMBIT()
{
    name();
}


void Pythia8::Hist::fill_GAMBIT(double x)
{
    fill(x);
}


void Pythia8::Hist::table_GAMBIT(std::ostream& os, bool printOverUnder) const
{
    table(os, printOverUnder);
}


void Pythia8::Hist::table_GAMBIT(std::ostream& os) const
{
    table(os);
}


void Pythia8::Hist::table_GAMBIT() const
{
    table();
}


void Pythia8::Hist::table_GAMBIT(std::string fileName, bool printOverUnder) const
{
    table(fileName, printOverUnder);
}


void Pythia8::Hist::table_GAMBIT(std::string fileName) const
{
    table(fileName);
}


bool Pythia8::Hist::sameSize_GAMBIT(const Pythia8::Abstract_Hist& h) const
{
    return sameSize(dynamic_cast< const Pythia8::Hist& >(h));
}


void Pythia8::Hist::takeLog_GAMBIT()
{
    takeLog();
}



Pythia8::Abstract_Hist* Pythia8::Hist::operator_equal_GAMBIT(const Pythia8::Abstract_Hist& h)
{
    return &(operator=(dynamic_cast< const Pythia8::Hist& >(h)));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_plus_equal_GAMBIT(const Pythia8::Abstract_Hist& h)
{
    return &(operator+=(dynamic_cast< const Pythia8::Hist& >(h)));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_minus_equal_GAMBIT(const Pythia8::Abstract_Hist& h)
{
    return &(operator-=(dynamic_cast< const Pythia8::Hist& >(h)));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_asterix_equal_GAMBIT(const Pythia8::Abstract_Hist& h)
{
    return &(operator*=(dynamic_cast< const Pythia8::Hist& >(h)));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_slash_equal_GAMBIT(const Pythia8::Abstract_Hist& h)
{
    return &(operator/=(dynamic_cast< const Pythia8::Hist& >(h)));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_plus_equal_GAMBIT(double f)
{
    return &(operator+=(f));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_minus_equal_GAMBIT(double f)
{
    return &(operator-=(f));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_asterix_equal_GAMBIT(double f)
{
    return &(operator*=(f));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_slash_equal_GAMBIT(double f)
{
    return &(operator/=(f));
}



Pythia8::Abstract_Hist* Pythia8::Hist::pointerCopy_GAMBIT() { return new Pythia8::Hist(*this); }
void Pythia8::Hist::pointerAssign_GAMBIT(Pythia8::Abstract_Hist* in) { *this = *dynamic_cast<Hist*>(in); }
