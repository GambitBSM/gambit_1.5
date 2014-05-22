#include "classes.hpp"

Abstract_T* Factory_T()
{
    return new T();
}

Abstract_T* Factory_T(int i_in, double d_in)
{
    return new T(i_in, d_in);
}
