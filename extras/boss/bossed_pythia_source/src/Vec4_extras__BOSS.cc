#include "Pythia8/Basics.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

void Pythia8::Vec4::p__BOSS(Pythia8::Abstract_Vec4& pIn)
{
    p(dynamic_cast< Pythia8::Vec4& >(pIn));
}


void Pythia8::Vec4::rotaxis__BOSS(double phiIn, const Pythia8::Abstract_Vec4& n)
{
    rotaxis(phiIn, dynamic_cast< const Pythia8::Vec4& >(n));
}


void Pythia8::Vec4::bst__BOSS(const Pythia8::Abstract_Vec4& pIn)
{
    bst(dynamic_cast< const Pythia8::Vec4& >(pIn));
}


void Pythia8::Vec4::bst__BOSS(const Pythia8::Abstract_Vec4& pIn, double mIn)
{
    bst(dynamic_cast< const Pythia8::Vec4& >(pIn), mIn);
}


void Pythia8::Vec4::bstback__BOSS(const Pythia8::Abstract_Vec4& pIn)
{
    bstback(dynamic_cast< const Pythia8::Vec4& >(pIn));
}


void Pythia8::Vec4::bstback__BOSS(const Pythia8::Abstract_Vec4& pIn, double mIn)
{
    bstback(dynamic_cast< const Pythia8::Vec4& >(pIn), mIn);
}



Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_equal__BOSS(const Pythia8::Abstract_Vec4& v)
{
    return &(operator=(dynamic_cast< const Pythia8::Vec4& >(v)));
}


Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_equal__BOSS(double value)
{
    return &(operator=(value));
}


Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_minus__BOSS()
{
    return new Pythia8::Vec4(operator-());
}


Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_plus_equal__BOSS(const Pythia8::Abstract_Vec4& v)
{
    return &(operator+=(dynamic_cast< const Pythia8::Vec4& >(v)));
}


Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_minus_equal__BOSS(const Pythia8::Abstract_Vec4& v)
{
    return &(operator-=(dynamic_cast< const Pythia8::Vec4& >(v)));
}


Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_asterix_equal__BOSS(double f)
{
    return &(operator*=(f));
}


Pythia8::Abstract_Vec4* Pythia8::Vec4::operator_slash_equal__BOSS(double f)
{
    return &(operator/=(f));
}



#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_Vec4* Pythia8::Vec4::pointerCopy__BOSS()
{
    Pythia8::Abstract_Vec4* new_ptr = new Pythia8::Vec4(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::Vec4::pointerAssign__BOSS(Pythia8::Abstract_Vec4* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Vec4* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<Vec4*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
