#include "Pythia8/Basics.h"
#include "backend_types/Pythia_8_186/wrapper_Hist_decl.h"
#include "backend_types/Pythia_8_186/wrapper_Hist_def.h"
#include <string>
#include "abstracttypedefs.h"
#include "wrappertypedefs.h"

// FACTORY_SIGNATURES_ORDER: ##()##(std::string, int, double, double)##(std::string, int, double)##(std::string, int)##(std::string)##(std::string, const Pythia8::Hist__BOSS&)##

namespace Pythia8
{
    Abstract_Hist* Factory_Hist()
    {
        return new Hist();
    }
    
    Abstract_Hist* Factory_Hist(std::string titleIn, int nBinIn, double xMinIn, double xMaxIn)
    {
        return new Hist(titleIn, nBinIn, xMinIn, xMaxIn);
    }
    
    Abstract_Hist* Factory_Hist(std::string titleIn, int nBinIn, double xMinIn)
    {
        return new Hist(titleIn, nBinIn, xMinIn);
    }
    
    Abstract_Hist* Factory_Hist(std::string titleIn, int nBinIn)
    {
        return new Hist(titleIn, nBinIn);
    }
    
    Abstract_Hist* Factory_Hist(std::string titleIn)
    {
        return new Hist(titleIn);
    }
    
    Abstract_Hist* Factory_Hist(std::string titleIn, const Pythia8::Hist__BOSS& h)
    {
        return new Hist(titleIn, dynamic_cast< const Pythia8::Hist& >(*h.BEptr));
    }
    
}

