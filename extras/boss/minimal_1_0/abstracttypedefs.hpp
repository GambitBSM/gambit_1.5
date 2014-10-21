#include "backend_types/BOSSMinimalExample_1_0/forward_decls_abstract_classes.hpp"
#include "backend_types/BOSSMinimalExample_1_0/identification.hpp"

namespace nspace3
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::nspace3::Abstract_Y Abstract_Y;
}

namespace nspace1
{
    namespace nspace2
    {
        typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::nspace1::nspace2::Abstract_X Abstract_X;
    }
}

#include "backend_undefs.hpp"
