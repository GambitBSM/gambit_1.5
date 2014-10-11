#include "classes.hpp"
#include "backend_types/BOSSMinimalExample_1_0/wrapper_Y_decl.hpp"
#include "backend_types/BOSSMinimalExample_1_0/wrapper_Y_def.hpp"
#include "backend_types/BOSSMinimalExample_1_0/wrapper_X_decl.hpp"
#include "backend_types/BOSSMinimalExample_1_0/wrapper_X_def.hpp"
#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"

namespace nspace3
{
    Abstract_Y* Factory_Y()
    {
        return new Y();
    }
    
    Abstract_Y* Factory_Y(nspace1::nspace2::X_GAMBIT& x_in)
    {
        return new Y(dynamic_cast< nspace1::nspace2::X& >(*x_in.BEptr));
    }
    
}

