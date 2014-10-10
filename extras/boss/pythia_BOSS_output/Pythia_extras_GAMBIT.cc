#include "backend_types/BOSSedPythia_1_0/abstract_Event.h"
#include "backend_types/BOSSedPythia_1_0/abstract_Info.h"
#include <string>
#include "backend_types/BOSSedPythia_1_0/abstract_Vec4.h"
#include <vector>
#include <istream>
#include <ostream>
#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"
#include "Pythia8/Pythia.h"

bool Pythia8::Pythia::readString_GAMBIT(std::string arg_1)
{
    return readString(arg_1);
}


bool Pythia8::Pythia::readFile_GAMBIT(std::string fileName, bool warn)
{
    return readFile(fileName, warn);
}


bool Pythia8::Pythia::readFile_GAMBIT(std::string fileName)
{
    return readFile(fileName);
}


bool Pythia8::Pythia::readFile_GAMBIT(std::istream& is, bool warn)
{
    return readFile(is, warn);
}


bool Pythia8::Pythia::readFile_GAMBIT(std::istream& is)
{
    return readFile(is);
}


bool Pythia8::Pythia::readFile_GAMBIT()
{
    return readFile();
}


bool Pythia8::Pythia::init_GAMBIT(std::string LesHouchesEventFile)
{
    return init(LesHouchesEventFile);
}


int Pythia8::Pythia::forceTimeShower_GAMBIT(int iBeg, int iEnd, double pTmax)
{
    return forceTimeShower(iBeg, iEnd, pTmax);
}


bool Pythia8::Pythia::forceHadronLevel_GAMBIT()
{
    return forceHadronLevel();
}


void Pythia8::Pythia::LHAeventList_GAMBIT()
{
    LHAeventList();
}


void Pythia8::Pythia::statistics_GAMBIT(bool all)
{
    statistics(all);
}


void Pythia8::Pythia::statistics_GAMBIT()
{
    statistics();
}


void Pythia8::Pythia::banner_GAMBIT()
{
    banner();
}


int Pythia8::Pythia::readSubrun_GAMBIT(std::string line, bool warn)
{
    return readSubrun(line, warn);
}


int Pythia8::Pythia::readSubrun_GAMBIT(std::string line)
{
    return readSubrun(line);
}


bool Pythia8::Pythia::check_GAMBIT()
{
    return check();
}



Pythia8::Abstract_Pythia* Pythia8::Pythia::operator_equal_GAMBIT(const Pythia8::Abstract_Pythia& arg_1)
{
    return &(operator=(dynamic_cast< const Pythia8::Pythia& >(arg_1)));
}


Pythia8::Abstract_Event& Pythia8::Pythia::process_ref_GAMBIT() { return process; }

Pythia8::Abstract_Event& Pythia8::Pythia::event_ref_GAMBIT() { return event; }

Pythia8::Abstract_Info& Pythia8::Pythia::info_ref_GAMBIT() { return info; }


Pythia8::Abstract_Pythia* Pythia8::Pythia::pointerCopy_GAMBIT() { return new Pythia8::Pythia(*this); }
void Pythia8::Pythia::pointerAssign_GAMBIT(Pythia8::Abstract_Pythia* in) { *this = *dynamic_cast<Pythia*>(in); }
