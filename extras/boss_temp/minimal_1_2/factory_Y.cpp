#include "classes.hpp"
#include "wrappers.hpp"

Abstract_Y* Factory_Y()
{
    return new Y();
}

Abstract_Y* Factory_Y(wrapper_X& x_in)
{
    return new Y(reinterpret_cast< X& >(*x_in.BEptr));
}

#include "backend_undefs.hpp"
