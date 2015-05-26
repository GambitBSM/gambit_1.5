#include "Pythia8/SusyLesHouches.h"
#include "backend_types/Pythia_8_209/wrapper_LHdecayChannel.h"
#include <vector>
#include <string>
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_LHdecayChannel* Factory_LHdecayChannel_0()
    {
        return new LHdecayChannel();
    }
    
    Abstract_LHdecayChannel* Factory_LHdecayChannel_1(double bratIn, int nDaIn, std::vector<int,std::allocator<int> > idDaIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > cIn)
    {
        return new LHdecayChannel(bratIn, nDaIn, idDaIn, cIn);
    }
    
    Abstract_LHdecayChannel* Factory_LHdecayChannel_2(double bratIn, int nDaIn, std::vector<int,std::allocator<int> > idDaIn)
    {
        return new LHdecayChannel(bratIn, nDaIn, idDaIn);
    }
    
}

