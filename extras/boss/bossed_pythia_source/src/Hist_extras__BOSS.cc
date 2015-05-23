#include <string>
#include <vector>
#include "Pythia8/Basics.h"
#include <ostream>
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

void Pythia8::Hist::book__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > titleIn, int nBinIn, double xMinIn)
{
    book(titleIn, nBinIn, xMinIn);
}


void Pythia8::Hist::book__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > titleIn, int nBinIn)
{
    book(titleIn, nBinIn);
}


void Pythia8::Hist::book__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > titleIn)
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


void Pythia8::Hist::table__BOSS(std::basic_ostream<char,std::char_traits<char> >& os, bool printOverUnder) const
{
    table(os, printOverUnder);
}


void Pythia8::Hist::table__BOSS(std::basic_ostream<char,std::char_traits<char> >& os) const
{
    table(os);
}


void Pythia8::Hist::table__BOSS() const
{
    table();
}


void Pythia8::Hist::table__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > fileName, bool printOverUnder) const
{
    table(fileName, printOverUnder);
}


void Pythia8::Hist::table__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > fileName) const
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



#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_Hist* Pythia8::Hist::pointerCopy__BOSS()
{
    Pythia8::Abstract_Hist* new_ptr = new Pythia8::Hist(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::Hist::pointerAssign__BOSS(Pythia8::Abstract_Hist* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Hist* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<Hist*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
