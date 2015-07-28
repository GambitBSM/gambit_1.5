#include "backend_types/Pythia_8_186/wrapper_Settings.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/StandardModel.h"

void Pythia8::AlphaEM::init__BOSS(int orderIn, Pythia8::Abstract_Settings* settingsPtr)
{
    init(orderIn, dynamic_cast< Pythia8::Settings* >(settingsPtr));
}




#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_AlphaEM* Pythia8::AlphaEM::pointerCopy__BOSS()
{
    Pythia8::Abstract_AlphaEM* new_ptr = new Pythia8::AlphaEM(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::AlphaEM::pointerAssign__BOSS(Pythia8::Abstract_AlphaEM* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::AlphaEM* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<AlphaEM*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
