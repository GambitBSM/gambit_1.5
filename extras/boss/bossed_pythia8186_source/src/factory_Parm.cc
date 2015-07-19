#include <string>
#include "Pythia8/Settings.h"
#include "backend_types/Pythia_8_186/wrapper_Parm.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_Parm* Factory_Parm_0(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, double defaultIn, bool hasMinIn, bool hasMaxIn, double minIn, double maxIn)
    {
        return new Parm(nameIn, defaultIn, hasMinIn, hasMaxIn, minIn, maxIn);
    }
    
    Abstract_Parm* Factory_Parm_1(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, double defaultIn, bool hasMinIn, bool hasMaxIn, double minIn)
    {
        return new Parm(nameIn, defaultIn, hasMinIn, hasMaxIn, minIn);
    }
    
    Abstract_Parm* Factory_Parm_2(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, double defaultIn, bool hasMinIn, bool hasMaxIn)
    {
        return new Parm(nameIn, defaultIn, hasMinIn, hasMaxIn);
    }
    
    Abstract_Parm* Factory_Parm_3(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, double defaultIn, bool hasMinIn)
    {
        return new Parm(nameIn, defaultIn, hasMinIn);
    }
    
    Abstract_Parm* Factory_Parm_4(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, double defaultIn)
    {
        return new Parm(nameIn, defaultIn);
    }
    
    Abstract_Parm* Factory_Parm_5(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn)
    {
        return new Parm(nameIn);
    }
    
    Abstract_Parm* Factory_Parm_6()
    {
        return new Parm();
    }
    
}

