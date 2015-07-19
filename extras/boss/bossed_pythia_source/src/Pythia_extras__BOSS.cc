#include <string>
#include <vector>
#include <istream>
#include <ostream>
#include "Pythia8/Pythia.h"
#include "backend_types/Pythia_8_209/wrapper_Event.h"
#include "backend_types/Pythia_8_209/wrapper_Info.h"
#include "backend_types/Pythia_8_209/wrapper_Settings.h"
#include "backend_types/Pythia_8_209/wrapper_ParticleData.h"
#include "backend_types/Pythia_8_209/wrapper_Rndm.h"
#include "backend_types/Pythia_8_209/wrapper_Couplings.h"
#include "backend_types/Pythia_8_209/wrapper_SLHAinterface.h"
#include "backend_types/Pythia_8_209/wrapper_Vec4.h"
#include "backend_types/Pythia_8_209/wrapper_BeamParticle.h"
#include "backend_types/Pythia_8_209/wrapper_UserHooks.h"
#include "backend_types/Pythia_8_209/wrapper_PartonLevel.h"
#include "backend_types/Pythia_8_209/wrapper_SigmaTotal.h"
#include "backend_types/Pythia_8_209/wrapper_SigmaProcess.h"
#include "backend_types/Pythia_8_209/wrapper_ResonanceWidths.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

bool Pythia8::Pythia::readString__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > arg_1)
{
    return readString(arg_1);
}


bool Pythia8::Pythia::readFile__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > fileName, bool warn)
{
    return readFile(fileName, warn);
}


bool Pythia8::Pythia::readFile__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > fileName)
{
    return readFile(fileName);
}


bool Pythia8::Pythia::readFile__BOSS(std::basic_istream<char,std::char_traits<char> >& is, bool warn)
{
    return readFile(is, warn);
}


bool Pythia8::Pythia::readFile__BOSS(std::basic_istream<char,std::char_traits<char> >& is)
{
    return readFile(is);
}


bool Pythia8::Pythia::readFile__BOSS()
{
    return readFile();
}


bool Pythia8::Pythia::setUserHooksPtr__BOSS(Pythia8::Abstract_UserHooks* userHooksPtrIn)
{
    return setUserHooksPtr(dynamic_cast< Pythia8::UserHooks* >(userHooksPtrIn));
}


bool Pythia8::Pythia::setSigmaPtr__BOSS(Pythia8::Abstract_SigmaProcess* sigmaPtrIn)
{
    return setSigmaPtr(dynamic_cast< Pythia8::SigmaProcess* >(sigmaPtrIn));
}


bool Pythia8::Pythia::setResonancePtr__BOSS(Pythia8::Abstract_ResonanceWidths* resonancePtrIn)
{
    return setResonancePtr(dynamic_cast< Pythia8::ResonanceWidths* >(resonancePtrIn));
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


void Pythia8::Pythia::banner__BOSS()
{
    banner();
}


int Pythia8::Pythia::readSubrun__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > line, bool warn)
{
    return readSubrun(line, warn);
}


int Pythia8::Pythia::readSubrun__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > line)
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

Pythia8::Abstract_Settings& Pythia8::Pythia::settings_ref__BOSS() { return settings; }

Pythia8::Abstract_ParticleData& Pythia8::Pythia::particleData_ref__BOSS() { return particleData; }

Pythia8::Abstract_Rndm& Pythia8::Pythia::rndm_ref__BOSS() { return rndm; }

Pythia8::Abstract_Couplings& Pythia8::Pythia::couplings_ref__BOSS() { return couplings; }

Pythia8::Abstract_SLHAinterface& Pythia8::Pythia::slhaInterface_ref__BOSS() { return slhaInterface; }


#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_Pythia* Pythia8::Pythia::pointerCopy__BOSS()
{
    Pythia8::Abstract_Pythia* new_ptr = new Pythia8::Pythia(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::Pythia::pointerAssign__BOSS(Pythia8::Abstract_Pythia* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Pythia* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<Pythia*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
