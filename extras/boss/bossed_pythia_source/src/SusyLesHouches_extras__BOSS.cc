#include <string>
#include <map>
#include <vector>
#include <istream>
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/SusyLesHouches.h"

int Pythia8::SusyLesHouches::readFile__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > slhaFileIn, int verboseIn)
{
    return readFile(slhaFileIn, verboseIn);
}


int Pythia8::SusyLesHouches::readFile__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > slhaFileIn)
{
    return readFile(slhaFileIn);
}


int Pythia8::SusyLesHouches::readFile__BOSS()
{
    return readFile();
}


int Pythia8::SusyLesHouches::readFile__BOSS(std::basic_istream<char,std::char_traits<char> >& arg_1, int verboseIn)
{
    return readFile(arg_1, verboseIn);
}


int Pythia8::SusyLesHouches::readFile__BOSS(std::basic_istream<char,std::char_traits<char> >& arg_1)
{
    return readFile(arg_1);
}


int Pythia8::SusyLesHouches::readSLHAea__BOSS(int verboseIn)
{
    return readSLHAea(verboseIn);
}


int Pythia8::SusyLesHouches::readSLHAea__BOSS()
{
    return readSLHAea();
}


void Pythia8::SusyLesHouches::printSpectrum__BOSS()
{
    printSpectrum();
}


void Pythia8::SusyLesHouches::message__BOSS(int arg_1, std::basic_string<char,std::char_traits<char>,std::allocator<char> > arg_2, std::basic_string<char,std::char_traits<char>,std::allocator<char> > arg_3)
{
    message(arg_1, arg_2, arg_3);
}



std::basic_string<char,std::char_traits<char>,std::allocator<char> >& Pythia8::SusyLesHouches::slhaFile_ref__BOSS() { return slhaFile; }

std::map<int,int,std::less<int>,std::allocator<std::pair<const int, int> > >& Pythia8::SusyLesHouches::decayIndices_ref__BOSS() { return decayIndices; }

std::vector<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::allocator<std::basic_string<char, std::char_traits<char>, std::allocator<char> > > >& Pythia8::SusyLesHouches::qnumbersName_ref__BOSS() { return qnumbersName; }

std::vector<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::allocator<std::basic_string<char, std::char_traits<char>, std::allocator<char> > > >& Pythia8::SusyLesHouches::qnumbersAntiName_ref__BOSS() { return qnumbersAntiName; }


#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_SusyLesHouches* Pythia8::SusyLesHouches::pointerCopy__BOSS()
{
    Pythia8::Abstract_SusyLesHouches* new_ptr = new Pythia8::SusyLesHouches(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::SusyLesHouches::pointerAssign__BOSS(Pythia8::Abstract_SusyLesHouches* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::SusyLesHouches* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<SusyLesHouches*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
