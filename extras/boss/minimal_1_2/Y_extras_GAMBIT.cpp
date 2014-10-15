#include "backend_types/BOSSMinimalExample_1_2/abstract_X.hpp"
#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"
#include "classes.hpp"

Abstract_X* Y::get_x_GAMBIT()
{
    return new X(get_x());
}


void Y::set_x_GAMBIT(Abstract_X& x_in)
{
    set_x(dynamic_cast< X& >(x_in));
}


void Y::set_x_ptr_GAMBIT(Abstract_X* x_in)
{
    set_x_ptr(dynamic_cast< X* >(x_in));
}


Abstract_X& Y::x_ref_GAMBIT() { return x; }


Abstract_Y* Y::pointerCopy_GAMBIT() { return new Y(*this); }
void Y::pointerAssign_GAMBIT(Abstract_Y* in) { *this = *dynamic_cast<Y*>(in); }
