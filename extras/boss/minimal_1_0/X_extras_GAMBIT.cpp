#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"
#include "classes.hpp"

nspace1::nspace2::Abstract_X* nspace1::nspace2::X::return_ref_this_GAMBIT()
{
    return &(return_ref_this());
}


nspace1::nspace2::Abstract_X* nspace1::nspace2::X::return_ptr_this_GAMBIT()
{
    return return_ptr_this();
}



nspace1::nspace2::Abstract_X* nspace1::nspace2::X::operator_plus_GAMBIT(nspace1::nspace2::Abstract_X& x_rhs)
{
    return new nspace1::nspace2::X(operator+(dynamic_cast< nspace1::nspace2::X& >(x_rhs)));
}


int& nspace1::nspace2::X::i_ref_GAMBIT() { return i; }


nspace1::nspace2::Abstract_X* nspace1::nspace2::X::pointerCopy_GAMBIT() { return new nspace1::nspace2::X(*this); }
void nspace1::nspace2::X::pointerAssign_GAMBIT(nspace1::nspace2::Abstract_X* in) { *this = *dynamic_cast<X*>(in); }
