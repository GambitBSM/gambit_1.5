#include "Pythia8/ParticleData.h"
#include "backend_types/Pythia_8_209/wrapper_DecayChannel.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_DecayChannel* Factory_DecayChannel_0(int onModeIn, double bRatioIn, int meModeIn, int prod0, int prod1, int prod2, int prod3, int prod4, int prod5, int prod6, int prod7)
    {
        return new DecayChannel(onModeIn, bRatioIn, meModeIn, prod0, prod1, prod2, prod3, prod4, prod5, prod6, prod7);
    }
    
    Abstract_DecayChannel* Factory_DecayChannel_1(int onModeIn, double bRatioIn, int meModeIn, int prod0, int prod1, int prod2, int prod3, int prod4, int prod5, int prod6)
    {
        return new DecayChannel(onModeIn, bRatioIn, meModeIn, prod0, prod1, prod2, prod3, prod4, prod5, prod6);
    }
    
    Abstract_DecayChannel* Factory_DecayChannel_2(int onModeIn, double bRatioIn, int meModeIn, int prod0, int prod1, int prod2, int prod3, int prod4, int prod5)
    {
        return new DecayChannel(onModeIn, bRatioIn, meModeIn, prod0, prod1, prod2, prod3, prod4, prod5);
    }
    
    Abstract_DecayChannel* Factory_DecayChannel_3(int onModeIn, double bRatioIn, int meModeIn, int prod0, int prod1, int prod2, int prod3, int prod4)
    {
        return new DecayChannel(onModeIn, bRatioIn, meModeIn, prod0, prod1, prod2, prod3, prod4);
    }
    
    Abstract_DecayChannel* Factory_DecayChannel_4(int onModeIn, double bRatioIn, int meModeIn, int prod0, int prod1, int prod2, int prod3)
    {
        return new DecayChannel(onModeIn, bRatioIn, meModeIn, prod0, prod1, prod2, prod3);
    }
    
    Abstract_DecayChannel* Factory_DecayChannel_5(int onModeIn, double bRatioIn, int meModeIn, int prod0, int prod1, int prod2)
    {
        return new DecayChannel(onModeIn, bRatioIn, meModeIn, prod0, prod1, prod2);
    }
    
    Abstract_DecayChannel* Factory_DecayChannel_6(int onModeIn, double bRatioIn, int meModeIn, int prod0, int prod1)
    {
        return new DecayChannel(onModeIn, bRatioIn, meModeIn, prod0, prod1);
    }
    
    Abstract_DecayChannel* Factory_DecayChannel_7(int onModeIn, double bRatioIn, int meModeIn, int prod0)
    {
        return new DecayChannel(onModeIn, bRatioIn, meModeIn, prod0);
    }
    
    Abstract_DecayChannel* Factory_DecayChannel_8(int onModeIn, double bRatioIn, int meModeIn)
    {
        return new DecayChannel(onModeIn, bRatioIn, meModeIn);
    }
    
    Abstract_DecayChannel* Factory_DecayChannel_9(int onModeIn, double bRatioIn)
    {
        return new DecayChannel(onModeIn, bRatioIn);
    }
    
    Abstract_DecayChannel* Factory_DecayChannel_10(int onModeIn)
    {
        return new DecayChannel(onModeIn);
    }
    
    Abstract_DecayChannel* Factory_DecayChannel_11()
    {
        return new DecayChannel();
    }
    
}

