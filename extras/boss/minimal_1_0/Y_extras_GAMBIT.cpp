#include "backend_types/BOSSMinimalExample_1_0/abstract_X.hpp"
#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"
#include "classes.hpp"

nspace1::nspace2::Abstract_X* nspace3::Y::get_x_GAMBIT()
{
    return new nspace1::nspace2::X(get_x());
}


void nspace3::Y::set_x_GAMBIT(nspace1::nspace2::Abstract_X& x_in)
{
    set_x(dynamic_cast< nspace1::nspace2::X& >(x_in));
}


void nspace3::Y::set_x_ptr_GAMBIT(nspace1::nspace2::Abstract_X* x_in)
{
    set_x_ptr(dynamic_cast< nspace1::nspace2::X* >(x_in));
}


nspace1::nspace2::Abstract_X& nspace3::Y::x_ref_GAMBIT() { return x; }


nspace3::Abstract_Y* nspace3::Y::pointerCopy_GAMBIT() { return new nspace3::Y(*this); }
void nspace3::Y::pointerAssign_GAMBIT(nspace3::Abstract_Y* in) { *this = *dynamic_cast<Y*>(in); }
