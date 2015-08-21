#include "Pythia8/SusyLesHouches.h"
#include "backend_types/Pythia_8_186/wrapper_LHdecayTable.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_LHdecayTable* Factory_LHdecayTable_0()
    {
        return new LHdecayTable();
    }
    
    Abstract_LHdecayTable* Factory_LHdecayTable_1(int idIn)
    {
        return new LHdecayTable(idIn);
    }
    
    Abstract_LHdecayTable* Factory_LHdecayTable_2(int idIn, double widthIn)
    {
        return new LHdecayTable(idIn, widthIn);
    }
    
}

