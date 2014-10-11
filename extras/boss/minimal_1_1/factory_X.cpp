#include "classes.hpp"
#include "backend_types/BOSSMinimalExample_1_1/wrapper_X_decl.hpp"
#include "backend_types/BOSSMinimalExample_1_1/wrapper_X_def.hpp"
#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"

Abstract_X* Factory_X()
{
    return new X();
}

Abstract_X* Factory_X(int i_in)
{
    return new X(i_in);
}


