#include <vector>
#include <string>
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/Basics.h"

void Pythia8::Rndm::init__BOSS()
{
    init();
}




#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_Rndm* Pythia8::Rndm::pointerCopy__BOSS()
{
    Pythia8::Abstract_Rndm* new_ptr = new Pythia8::Rndm(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::Rndm::pointerAssign__BOSS(Pythia8::Abstract_Rndm* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Rndm* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<Rndm*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
