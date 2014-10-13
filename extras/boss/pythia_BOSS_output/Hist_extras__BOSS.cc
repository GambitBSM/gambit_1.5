#include <string>
#include <vector>
#include <ostream>
#include "abstracttypedefs.h"
#include "wrappertypedefs.h"
#include "Pythia8/Basics.h"

void Pythia8::Hist::book__BOSS(std::string titleIn, int nBinIn, double xMinIn)
{
    book(titleIn, nBinIn, xMinIn);
}


void Pythia8::Hist::book__BOSS(std::string titleIn, int nBinIn)
{
    book(titleIn, nBinIn);
}


void Pythia8::Hist::book__BOSS(std::string titleIn)
{
    book(titleIn);
}


void Pythia8::Hist::book__BOSS()
{
    book();
}


void Pythia8::Hist::name__BOSS()
{
    name();
}


void Pythia8::Hist::fill__BOSS(double x)
{
    fill(x);
}


void Pythia8::Hist::table__BOSS(std::ostream& os, bool printOverUnder) const
{
    table(os, printOverUnder);
}


void Pythia8::Hist::table__BOSS(std::ostream& os) const
{
    table(os);
}


void Pythia8::Hist::table__BOSS() const
{
    table();
}


void Pythia8::Hist::table__BOSS(std::string fileName, bool printOverUnder) const
{
    table(fileName, printOverUnder);
}


void Pythia8::Hist::table__BOSS(std::string fileName) const
{
    table(fileName);
}


bool Pythia8::Hist::sameSize__BOSS(const Pythia8::Abstract_Hist& h) const
{
    return sameSize(dynamic_cast< const Pythia8::Hist& >(h));
}


void Pythia8::Hist::takeLog__BOSS()
{
    takeLog();
}



Pythia8::Abstract_Hist* Pythia8::Hist::operator_equal__BOSS(const Pythia8::Abstract_Hist& h)
{
    return &(operator=(dynamic_cast< const Pythia8::Hist& >(h)));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_plus_equal__BOSS(const Pythia8::Abstract_Hist& h)
{
    return &(operator+=(dynamic_cast< const Pythia8::Hist& >(h)));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_minus_equal__BOSS(const Pythia8::Abstract_Hist& h)
{
    return &(operator-=(dynamic_cast< const Pythia8::Hist& >(h)));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_asterix_equal__BOSS(const Pythia8::Abstract_Hist& h)
{
    return &(operator*=(dynamic_cast< const Pythia8::Hist& >(h)));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_slash_equal__BOSS(const Pythia8::Abstract_Hist& h)
{
    return &(operator/=(dynamic_cast< const Pythia8::Hist& >(h)));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_plus_equal__BOSS(double f)
{
    return &(operator+=(f));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_minus_equal__BOSS(double f)
{
    return &(operator-=(f));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_asterix_equal__BOSS(double f)
{
    return &(operator*=(f));
}


Pythia8::Abstract_Hist* Pythia8::Hist::operator_slash_equal__BOSS(double f)
{
    return &(operator/=(f));
}



Pythia8::Abstract_Hist* Pythia8::Hist::pointerCopy__BOSS() { return new Pythia8::Hist(*this); }
void Pythia8::Hist::pointerAssign__BOSS(Pythia8::Abstract_Hist* in) { *this = *dynamic_cast<Hist*>(in); }
