#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"
#include "Pythia8/Basics.h"

void Pythia8::Vec4::p_GAMBIT(Pythia8::Abstract_Vec4& pIn)
{
    p(dynamic_cast< Pythia8::Vec4& >(pIn));
}


void Pythia8::Vec4::rotaxis_GAMBIT(double phiIn, const Pythia8::Abstract_Vec4& n)
{
    rotaxis(phiIn, dynamic_cast< const Pythia8::Vec4& >(n));
}


void Pythia8::Vec4::bst_GAMBIT(const Pythia8::Abstract_Vec4& pIn)
{
    bst(dynamic_cast< const Pythia8::Vec4& >(pIn));
}


void Pythia8::Vec4::bst_GAMBIT(const Pythia8::Abstract_Vec4& pIn, double mIn)
{
    bst(dynamic_cast< const Pythia8::Vec4& >(pIn), mIn);
}


void Pythia8::Vec4::bstback_GAMBIT(const Pythia8::Abstract_Vec4& pIn)
{
    bstback(dynamic_cast< const Pythia8::Vec4& >(pIn));
}


void Pythia8::Vec4::bstback_GAMBIT(const Pythia8::Abstract_Vec4& pIn, double mIn)
{
    bstback(dynamic_cast< const Pythia8::Vec4& >(pIn), mIn);
}



Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_equal_GAMBIT(const Pythia8::Abstract_Vec4& v)
{
    return &(operator=(dynamic_cast< const Pythia8::Vec4& >(v)));
}


Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_equal_GAMBIT(double value)
{
    return &(operator=(value));
}


Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_minus_GAMBIT()
{
    return new Pythia8::Vec4(operator-());
}


Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_plus_equal_GAMBIT(const Pythia8::Abstract_Vec4& v)
{
    return &(operator+=(dynamic_cast< const Pythia8::Vec4& >(v)));
}


Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_minus_equal_GAMBIT(const Pythia8::Abstract_Vec4& v)
{
    return &(operator-=(dynamic_cast< const Pythia8::Vec4& >(v)));
}


Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_asterix_equal_GAMBIT(double f)
{
    return &(operator*=(f));
}


Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_slash_equal_GAMBIT(double f)
{
    return &(operator/=(f));
}



Pythia8::Abstract_Vec4* Pythia8::Vec4::pointerCopy_GAMBIT() { return new Pythia8::Vec4(*this); }
void Pythia8::Vec4::pointerAssign_GAMBIT(Pythia8::Abstract_Vec4* in) { *this = *dynamic_cast<Vec4*>(in); }
