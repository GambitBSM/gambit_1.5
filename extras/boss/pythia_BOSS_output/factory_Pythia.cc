#include <string>
#include "Pythia8/Pythia.h"
#include "backend_types/Pythia_8_186/wrapper_Pythia_decl.h"
#include "backend_types/Pythia_8_186/wrapper_Pythia_def.h"
#include "abstracttypedefs.h"
#include "wrappertypedefs.h"

// FACTORY_SIGNATURES_ORDER: ##(std::string, bool)##(std::string)##()##

namespace Pythia8
{
    Abstract_Pythia* Factory_Pythia(std::string xmlDir, bool printBanner)
    {
        return new Pythia(xmlDir, printBanner);
    }
    
    Abstract_Pythia* Factory_Pythia(std::string xmlDir)
    {
        return new Pythia(xmlDir);
    }
    
    Abstract_Pythia* Factory_Pythia()
    {
        return new Pythia();
    }
    
}

