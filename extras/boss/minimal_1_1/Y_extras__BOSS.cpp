#include "backend_types/BOSSMinimalExample_1_1/abstract_X.hpp"
#include "abstracttypedefs.hpp"
#include "wrappertypedefs.hpp"
#include "classes.hpp"

Abstract_X* Y::get_x__BOSS()
{
    return new X(get_x());
}


void Y::set_x__BOSS(Abstract_X& x_in)
{
    set_x(dynamic_cast< X& >(x_in));
}


void Y::set_x_ptr__BOSS(Abstract_X* x_in)
{
    set_x_ptr(dynamic_cast< X* >(x_in));
}


Abstract_X& Y::x_ref__BOSS() { return x; }


Abstract_Y* Y::pointerCopy__BOSS() { return new Y(*this); }
void Y::pointerAssign__BOSS(Abstract_Y* in) { *this = *dynamic_cast<Y*>(in); }
