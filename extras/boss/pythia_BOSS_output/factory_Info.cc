#include "Pythia8/Info.h"
#include "backend_types/Pythia_8_186/wrapper_Info_decl.h"
#include "backend_types/Pythia_8_186/wrapper_Info_def.h"
#include "abstracttypedefs.h"
#include "wrappertypedefs.h"

// FACTORY_SIGNATURES_ORDER: ##()##

namespace Pythia8
{
    Abstract_Info* Factory_Info()
    {
        return new Info();
    }
    
}

