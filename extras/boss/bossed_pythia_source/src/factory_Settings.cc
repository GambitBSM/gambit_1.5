#include "Pythia8/Settings.h"
#include "backend_types/Pythia_8_209/wrapper_Settings.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_Settings* Factory_Settings_0()
    {
        return new Settings();
    }
    
}

