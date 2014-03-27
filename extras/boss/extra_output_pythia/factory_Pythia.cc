#include "Pythia8/Pythia.h"

Pythia8::Abstract__Pythia* Factory_Pythia(std::string xmlDir, bool printBanner)
{
    return new Pythia8::Pythia( xmlDir,  printBanner);
}


