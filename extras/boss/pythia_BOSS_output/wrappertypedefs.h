#include "backend_types/Pythia_8_186/forward_decls_wrapper_classes.h"
#include "backend_types/Pythia_8_186/identification.hpp"

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Particle Particle__BOSS;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Info Info__BOSS;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Vec4 Vec4__BOSS;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Hist Hist__BOSS;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Event Event__BOSS;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Pythia Pythia__BOSS;
}

#include "backend_undefs.hpp"
