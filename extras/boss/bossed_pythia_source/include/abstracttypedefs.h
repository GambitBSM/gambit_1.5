#include "backend_types/Pythia_8_186/forward_decls_abstract_classes.h"
#include "backend_types/Pythia_8_186/identification.hpp"

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Abstract_Particle Abstract_Particle;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Abstract_Info Abstract_Info;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Abstract_Vec4 Abstract_Vec4;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Abstract_Hist Abstract_Hist;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Abstract_Event Abstract_Event;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Abstract_Pythia Abstract_Pythia;
}

#include "backend_undefs.hpp"
