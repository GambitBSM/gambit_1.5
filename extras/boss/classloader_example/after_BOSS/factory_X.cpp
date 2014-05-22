#include "classes.hpp"

Abstract_X* Factory_X()
{
    return new X();
}

Abstract_X* Factory_X(Abstract_T& t_in)
{
    return new X( dynamic_cast<T&>( t_in ) );
}
