#include "classes.hpp"
#include "backend_types/BOSSMinimalExample_1_1/wrapper_Y_decl.hpp"
#include "backend_types/BOSSMinimalExample_1_1/wrapper_Y_def.hpp"
#include "backend_types/BOSSMinimalExample_1_1/wrapper_X_decl.hpp"
#include "backend_types/BOSSMinimalExample_1_1/wrapper_X_def.hpp"
#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"

Abstract_Y* Factory_Y()
{
    return new Y();
}

Abstract_Y* Factory_Y(X_GAMBIT& x_in)
{
    return new Y(dynamic_cast< X& >(*x_in.BEptr));
}


