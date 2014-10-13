#include "Pythia8/Basics.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Hist_decl.h"
#include "backend_types/BOSSedPythia_1_0/wrapper_Hist_def.h"
#include <string>
#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"

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

