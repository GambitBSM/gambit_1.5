#include "Pythia8/Analysis.h"
#include "backend_types/Pythia_8_209/wrapper_SlowJet.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_SlowJet* Factory_SlowJet_0(int powerIn, double Rin, double pTjetMinIn, double etaMaxIn, int selectIn, int massSetIn)
    {
        return new SlowJet(powerIn, Rin, pTjetMinIn, etaMaxIn, selectIn, massSetIn);
    }
    
    Abstract_SlowJet* Factory_SlowJet_1(int powerIn, double Rin, double pTjetMinIn, double etaMaxIn, int selectIn)
    {
        return new SlowJet(powerIn, Rin, pTjetMinIn, etaMaxIn, selectIn);
    }
    
    Abstract_SlowJet* Factory_SlowJet_2(int powerIn, double Rin, double pTjetMinIn, double etaMaxIn)
    {
        return new SlowJet(powerIn, Rin, pTjetMinIn, etaMaxIn);
    }
    
    Abstract_SlowJet* Factory_SlowJet_3(int powerIn, double Rin, double pTjetMinIn)
    {
        return new SlowJet(powerIn, Rin, pTjetMinIn);
    }
    
    Abstract_SlowJet* Factory_SlowJet_4(int powerIn, double Rin)
    {
        return new SlowJet(powerIn, Rin);
    }
    
}

