#include "classes.hpp"
#include "backend_types/BOSSMinimalExample_1_1/wrapper_Y_decl.hpp"
#include "backend_types/BOSSMinimalExample_1_1/wrapper_Y_def.hpp"
#include "backend_types/BOSSMinimalExample_1_1/wrapper_X_decl.hpp"
#include "backend_types/BOSSMinimalExample_1_1/wrapper_X_def.hpp"
#include "abstracttypedefs.hpp"
#include "wrappertypedefs.hpp"

// FACTORY_SIGNATURES_ORDER: ##()##(X__BOSS&)##

Abstract_Y* Factory_Y()
{
    return new Y();
}

Abstract_Y* Factory_Y(X__BOSS& x_in)
{
    return new Y(dynamic_cast< X& >(*x_in.BEptr));
}


