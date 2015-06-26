#include "Pythia8/Basics.h"
#include "backend_types/Pythia_8_209/wrapper_Vec4.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_Vec4* Factory_Vec4_0(double xIn, double yIn, double zIn, double tIn)
    {
        return new Vec4(xIn, yIn, zIn, tIn);
    }
    
    Abstract_Vec4* Factory_Vec4_1(double xIn, double yIn, double zIn)
    {
        return new Vec4(xIn, yIn, zIn);
    }
    
    Abstract_Vec4* Factory_Vec4_2(double xIn, double yIn)
    {
        return new Vec4(xIn, yIn);
    }
    
    Abstract_Vec4* Factory_Vec4_3(double xIn)
    {
        return new Vec4(xIn);
    }
    
    Abstract_Vec4* Factory_Vec4_4()
    {
        return new Vec4();
    }
    
}

