#include "Pythia8/Basics.h"
#include "backend_types/Pythia_8_186/wrapper_Hist_decl.h"
#include "backend_types/Pythia_8_186/wrapper_Hist_def.h"
#include <string>
#include "abstracttypedefs.h"
#include "wrappertypedefs.h"

namespace Pythia8
{
    Abstract_Hist* Factory_Hist_0()
    {
        return new Hist();
    }
    
    Abstract_Hist* Factory_Hist_1(std::string titleIn, int nBinIn, double xMinIn, double xMaxIn)
    {
        return new Hist(titleIn, nBinIn, xMinIn, xMaxIn);
    }
    
    Abstract_Hist* Factory_Hist_2(std::string titleIn, int nBinIn, double xMinIn)
    {
        return new Hist(titleIn, nBinIn, xMinIn);
    }
    
    Abstract_Hist* Factory_Hist_3(std::string titleIn, int nBinIn)
    {
        return new Hist(titleIn, nBinIn);
    }
    
    Abstract_Hist* Factory_Hist_4(std::string titleIn)
    {
        return new Hist(titleIn);
    }
    
    Abstract_Hist* Factory_Hist_5(std::string titleIn, const Pythia8::Hist__BOSS& h)
    {
        return new Hist(titleIn, dynamic_cast< const Pythia8::Hist& >(*h.BEptr));
    }
    
}

