#include "Pythia8/Info.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Info_decl.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Info_def.h"
#include "abstracttypedefs.h"
#include "wrappertypedefs.h"

namespace Pythia8
{
    Abstract_Info* Factory_Info()
    {
        return new Info();
    }
    
}

