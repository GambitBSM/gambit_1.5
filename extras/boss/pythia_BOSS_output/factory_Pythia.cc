#include <string>
#include "Pythia8/Pythia.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Pythia_decl.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Pythia_def.h"
#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"

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

