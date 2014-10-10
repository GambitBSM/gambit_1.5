#include "backend_types/BOSSedPythia_1_0/forward_decls_wrapper_classes.h"
#include "backend_types/BOSSedPythia_1_0/identification.hpp"

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Particle Particle_GAMBIT;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Info Info_GAMBIT;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Vec4 Vec4_GAMBIT;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Hist Hist_GAMBIT;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Event Event_GAMBIT;
}

namespace Pythia8
{
    typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Pythia Pythia_GAMBIT;
}

#include "backend_undefs.hpp"
