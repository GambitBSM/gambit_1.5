#include "backend_types/BOSSMinimalExample_1_0/forward_decls_wrapper_classes.hpp"
#include "backend_types/BOSSMinimalExample_1_0/identification.hpp"

namespace nspace3
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::nspace3::Y Y_GAMBIT;
}

namespace nspace1
{
    namespace nspace2
    {
        typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::nspace1::nspace2::X X_GAMBIT;
    }
}

#include "backend_undefs.hpp"
