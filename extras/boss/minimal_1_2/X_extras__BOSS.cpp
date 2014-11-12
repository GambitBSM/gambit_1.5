#include "abstracttypedefs.hpp"
#include "wrappertypedefs.hpp"
#include "classes.hpp"

Abstract_X* X::return_ref_this__BOSS()
{
    return &(return_ref_this());
}


Abstract_X* X::return_ptr_this__BOSS()
{
    return return_ptr_this();
}



Abstract_X* X::operator_plus__BOSS(Abstract_X& x_rhs)
{
    return new X(operator+(dynamic_cast< X& >(x_rhs)));
}


int& X::i_ref__BOSS() { return i; }


Abstract_X* X::pointerCopy__BOSS() { return new X(*this); }
void X::pointerAssign__BOSS(Abstract_X* in) { *this = *dynamic_cast<X*>(in); }
