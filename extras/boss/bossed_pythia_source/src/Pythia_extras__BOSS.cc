#include "backend_types/Pythia_8_186/abstract_Event.h"
#include "backend_types/Pythia_8_186/abstract_Info.h"
#include <string>
#include "backend_types/Pythia_8_186/abstract_Vec4.h"
#include <vector>
#include <istream>
#include <ostream>
#include "abstracttypedefs.h"
#include "wrappertypedefs.h"
#include "Pythia8/Pythia.h"

bool Pythia8::Pythia::readString__BOSS(std::string arg_1)
{
    return readString(arg_1);
}


bool Pythia8::Pythia::readFile__BOSS(std::string fileName, bool warn)
{
    return readFile(fileName, warn);
}


bool Pythia8::Pythia::readFile__BOSS(std::string fileName)
{
    return readFile(fileName);
}


bool Pythia8::Pythia::readFile__BOSS(std::istream& is, bool warn)
{
    return readFile(is, warn);
}


bool Pythia8::Pythia::readFile__BOSS(std::istream& is)
{
    return readFile(is);
}


bool Pythia8::Pythia::readFile__BOSS()
{
    return readFile();
}


bool Pythia8::Pythia::init__BOSS(std::string LesHouchesEventFile)
{
    return init(LesHouchesEventFile);
}


int Pythia8::Pythia::forceTimeShower__BOSS(int iBeg, int iEnd, double pTmax)
{
    return forceTimeShower(iBeg, iEnd, pTmax);
}


bool Pythia8::Pythia::forceHadronLevel__BOSS()
{
    return forceHadronLevel();
}


void Pythia8::Pythia::LHAeventList__BOSS()
{
    LHAeventList();
}


void Pythia8::Pythia::statistics__BOSS(bool all)
{
    statistics(all);
}


void Pythia8::Pythia::statistics__BOSS()
{
    statistics();
}


void Pythia8::Pythia::banner__BOSS()
{
    banner();
}


int Pythia8::Pythia::readSubrun__BOSS(std::string line, bool warn)
{
    return readSubrun(line, warn);
}


int Pythia8::Pythia::readSubrun__BOSS(std::string line)
{
    return readSubrun(line);
}


bool Pythia8::Pythia::check__BOSS()
{
    return check();
}



Pythia8::Abstract_Pythia* Pythia8::Pythia::operator_equal__BOSS(const Pythia8::Abstract_Pythia& arg_1)
{
    return &(operator=(dynamic_cast< const Pythia8::Pythia& >(arg_1)));
}


Pythia8::Abstract_Event& Pythia8::Pythia::process_ref__BOSS() { return process; }

Pythia8::Abstract_Event& Pythia8::Pythia::event_ref__BOSS() { return event; }

Pythia8::Abstract_Info& Pythia8::Pythia::info_ref__BOSS() { return info; }


Pythia8::Abstract_Pythia* Pythia8::Pythia::pointerCopy__BOSS() { return new Pythia8::Pythia(*this); }
void Pythia8::Pythia::pointerAssign__BOSS(Pythia8::Abstract_Pythia* in) { *this = *dynamic_cast<Pythia*>(in); }
