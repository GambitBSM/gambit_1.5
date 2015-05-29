#include "backend_types/Pythia_8_186/wrapper_Rndm.h"
#include "backend_types/Pythia_8_186/wrapper_Settings.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/StandardModel.h"

void Pythia8::CoupSM::init__BOSS(Pythia8::Abstract_Settings& settings, Pythia8::Abstract_Rndm* rndmPtrIn)
{
    init(dynamic_cast< Pythia8::Settings& >(settings), dynamic_cast< Pythia8::Rndm* >(rndmPtrIn));
}




#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_CoupSM* Pythia8::CoupSM::pointerCopy__BOSS()
{
    Pythia8::Abstract_CoupSM* new_ptr = new Pythia8::CoupSM(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::CoupSM::pointerAssign__BOSS(Pythia8::Abstract_CoupSM* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::CoupSM* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<CoupSM*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
