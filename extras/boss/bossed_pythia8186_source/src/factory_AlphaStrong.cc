#include "Pythia8/StandardModel.h"
#include "backend_types/Pythia_8_186/wrapper_AlphaStrong.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_AlphaStrong* Factory_AlphaStrong_0()
    {
        return new AlphaStrong();
    }
    
    Abstract_AlphaStrong* Factory_AlphaStrong_1(double valueIn, int orderIn)
    {
        return new AlphaStrong(valueIn, orderIn);
    }
    
    Abstract_AlphaStrong* Factory_AlphaStrong_2(double valueIn)
    {
        return new AlphaStrong(valueIn);
    }
    
}

