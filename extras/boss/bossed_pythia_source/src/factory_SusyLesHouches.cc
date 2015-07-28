#include "Pythia8/SusyLesHouches.h"
#include "backend_types/Pythia_8_209/wrapper_SusyLesHouches.h"
#include <string>
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_SusyLesHouches* Factory_SusyLesHouches_0(int verboseIn)
    {
        return new SusyLesHouches(verboseIn);
    }
    
    Abstract_SusyLesHouches* Factory_SusyLesHouches_1()
    {
        return new SusyLesHouches();
    }
    
    Abstract_SusyLesHouches* Factory_SusyLesHouches_2(std::basic_string<char,std::char_traits<char>,std::allocator<char> > filename, int verboseIn)
    {
        return new SusyLesHouches(filename, verboseIn);
    }
    
    Abstract_SusyLesHouches* Factory_SusyLesHouches_3(std::basic_string<char,std::char_traits<char>,std::allocator<char> > filename)
    {
        return new SusyLesHouches(filename);
    }
    
}

