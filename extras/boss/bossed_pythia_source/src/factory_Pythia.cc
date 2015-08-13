#include <string>
#include "Pythia8/Pythia.h"
#include "backend_types/Pythia_8_209/wrapper_Pythia.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_Pythia* Factory_Pythia_0(std::basic_string<char,std::char_traits<char>,std::allocator<char> > xmlDir, bool printBanner)
    {
        return new Pythia(xmlDir, printBanner);
    }
    
    Abstract_Pythia* Factory_Pythia_1(std::basic_string<char,std::char_traits<char>,std::allocator<char> > xmlDir)
    {
        return new Pythia(xmlDir);
    }
    
    Abstract_Pythia* Factory_Pythia_2()
    {
        return new Pythia();
    }
    
}

