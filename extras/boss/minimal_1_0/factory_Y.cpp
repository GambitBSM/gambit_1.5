#include "classes.hpp"
#include "backend_types/BOSSMinimalExample_1_0/wrapper_Y_decl.hpp"
#include "backend_types/BOSSMinimalExample_1_0/wrapper_Y_def.hpp"
#include "backend_types/BOSSMinimalExample_1_0/wrapper_X_decl.hpp"
#include "backend_types/BOSSMinimalExample_1_0/wrapper_X_def.hpp"
#include "abstracttypedefs.hpp"
#include "wrappertypedefs.hpp"

// FACTORY_SIGNATURES_ORDER: ##()##(nspace1::nspace2::X__BOSS&)##

namespace nspace3
{
    Abstract_Y* Factory_Y()
    {
        return new Y();
    }
    
    Abstract_Y* Factory_Y(nspace1::nspace2::X__BOSS& x_in)
    {
        return new Y(dynamic_cast< nspace1::nspace2::X& >(*x_in.BEptr));
    }
    
}

