#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/ResonanceWidths.h"

void Pythia8::ResonanceGmZ::calcPreFac__BOSS()
{
    calcPreFac();
}


void Pythia8::ResonanceGmZ::calcWidth__BOSS()
{
    calcWidth();
}




#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_ResonanceGmZ* Pythia8::ResonanceGmZ::pointerCopy__BOSS()
{
    Pythia8::Abstract_ResonanceGmZ* new_ptr = new Pythia8::ResonanceGmZ(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::ResonanceGmZ::pointerAssign__BOSS(Pythia8::Abstract_ResonanceGmZ* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::ResonanceGmZ* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<ResonanceGmZ*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
