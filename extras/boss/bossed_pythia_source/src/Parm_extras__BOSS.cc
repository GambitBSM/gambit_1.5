#include <string>
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/Settings.h"


std::basic_string<char,std::char_traits<char>,std::allocator<char> >& Pythia8::Parm::name_ref__BOSS() { return name; }

double& Pythia8::Parm::valNow_ref__BOSS() { return valNow; }

double& Pythia8::Parm::valDefault_ref__BOSS() { return valDefault; }

bool& Pythia8::Parm::hasMin_ref__BOSS() { return hasMin; }

bool& Pythia8::Parm::hasMax_ref__BOSS() { return hasMax; }

double& Pythia8::Parm::valMin_ref__BOSS() { return valMin; }

double& Pythia8::Parm::valMax_ref__BOSS() { return valMax; }


#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_Parm* Pythia8::Parm::pointerCopy__BOSS()
{
    Pythia8::Abstract_Parm* new_ptr = new Pythia8::Parm(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::Parm::pointerAssign__BOSS(Pythia8::Abstract_Parm* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Parm* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<Parm*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
