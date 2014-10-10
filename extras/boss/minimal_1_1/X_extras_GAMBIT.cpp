#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"
#include "classes.hpp"

Abstract_X* X::return_ref_this_GAMBIT()
{
    return &(return_ref_this());
}


Abstract_X* X::return_ptr_this_GAMBIT()
{
    return return_ptr_this();
}



Abstract_X* X::operator_plus_GAMBIT(Abstract_X& x_rhs)
{
    return new X(operator+(dynamic_cast< X& >(x_rhs)));
}


int& X::i_ref_GAMBIT() { return i; }


Abstract_X* X::pointerCopy_GAMBIT() { return new X(*this); }
void X::pointerAssign_GAMBIT(Abstract_X* in) { *this = *dynamic_cast<X*>(in); }
