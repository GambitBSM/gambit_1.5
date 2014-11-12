#include "backend_types/BOSSMinimalExample_1_0/abstract_X.hpp"
#include "abstracttypedefs.hpp"
#include "wrappertypedefs.hpp"
#include "classes.hpp"

nspace1::nspace2::Abstract_X* nspace3::Y::get_x__BOSS()
{
    return new nspace1::nspace2::X(get_x());
}


void nspace3::Y::set_x__BOSS(nspace1::nspace2::Abstract_X& x_in)
{
    set_x(dynamic_cast< nspace1::nspace2::X& >(x_in));
}


void nspace3::Y::set_x_ptr__BOSS(nspace1::nspace2::Abstract_X* x_in)
{
    set_x_ptr(dynamic_cast< nspace1::nspace2::X* >(x_in));
}


nspace1::nspace2::Abstract_X& nspace3::Y::x_ref__BOSS() { return x; }


nspace3::Abstract_Y* nspace3::Y::pointerCopy__BOSS() { return new nspace3::Y(*this); }
void nspace3::Y::pointerAssign__BOSS(nspace3::Abstract_Y* in) { *this = *dynamic_cast<Y*>(in); }
