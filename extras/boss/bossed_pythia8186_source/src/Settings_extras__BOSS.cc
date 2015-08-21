#include <string>
#include <ostream>
#include <vector>
#include "backend_types/Pythia_8_186/wrapper_Info.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/Settings.h"

void Pythia8::Settings::initPtr__BOSS(Pythia8::Abstract_Info* infoPtrIn)
{
    initPtr(dynamic_cast< Pythia8::Info* >(infoPtrIn));
}


bool Pythia8::Settings::init__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > startFile, bool append)
{
    return init(startFile, append);
}


bool Pythia8::Settings::init__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > startFile)
{
    return init(startFile);
}


bool Pythia8::Settings::init__BOSS()
{
    return init();
}


bool Pythia8::Settings::reInit__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > startFile)
{
    return reInit(startFile);
}


bool Pythia8::Settings::reInit__BOSS()
{
    return reInit();
}


bool Pythia8::Settings::readString__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > line, bool warn)
{
    return readString(line, warn);
}


bool Pythia8::Settings::readString__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > line)
{
    return readString(line);
}


bool Pythia8::Settings::writeFile__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > toFile)
{
    return writeFile(toFile);
}


bool Pythia8::Settings::writeFile__BOSS(std::basic_ostream<char,std::char_traits<char> >& os)
{
    return writeFile(os);
}


bool Pythia8::Settings::writeFile__BOSS()
{
    return writeFile();
}


void Pythia8::Settings::listAll__BOSS()
{
    listAll();
}


void Pythia8::Settings::listChanged__BOSS()
{
    listChanged();
}


void Pythia8::Settings::list__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > match)
{
    list(match);
}


void Pythia8::Settings::list__BOSS(bool doListAll, bool doListString, std::basic_string<char,std::char_traits<char>,std::allocator<char> > match)
{
    list(doListAll, doListString, match);
}




#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_Settings* Pythia8::Settings::pointerCopy__BOSS()
{
    Pythia8::Abstract_Settings* new_ptr = new Pythia8::Settings(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::Settings::pointerAssign__BOSS(Pythia8::Abstract_Settings* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Settings* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<Settings*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
