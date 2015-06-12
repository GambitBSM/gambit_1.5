#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/StandardModel.h"


bool& Pythia8::Couplings::isSUSY_ref__BOSS() { return isSUSY; }


#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_Couplings* Pythia8::Couplings::pointerCopy__BOSS()
{
    Pythia8::Abstract_Couplings* new_ptr = new Pythia8::Couplings(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::Couplings::pointerAssign__BOSS(Pythia8::Abstract_Couplings* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Couplings* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<Couplings*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
