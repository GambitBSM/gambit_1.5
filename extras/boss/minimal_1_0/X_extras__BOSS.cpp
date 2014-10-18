#include "abstracttypedefs.hpp"
#include "wrappertypedefs.hpp"
#include "classes.hpp"

nspace1::nspace2::Abstract_X* nspace1::nspace2::X::return_ref_this__BOSS()
{
    return &(return_ref_this());
}


nspace1::nspace2::Abstract_X* nspace1::nspace2::X::return_ptr_this__BOSS()
{
    return return_ptr_this();
}



nspace1::nspace2::Abstract_X* nspace1::nspace2::X::operator_plus__BOSS(nspace1::nspace2::Abstract_X& x_rhs)
{
    return new nspace1::nspace2::X(operator+(dynamic_cast< nspace1::nspace2::X& >(x_rhs)));
}


int& nspace1::nspace2::X::i_ref__BOSS() { return i; }


nspace1::nspace2::Abstract_X* nspace1::nspace2::X::pointerCopy__BOSS() { return new nspace1::nspace2::X(*this); }
void nspace1::nspace2::X::pointerAssign__BOSS(nspace1::nspace2::Abstract_X* in) { *this = *dynamic_cast<X*>(in); }
