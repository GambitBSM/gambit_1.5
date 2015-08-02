#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/StandardModel.h"

void Pythia8::AlphaStrong::init__BOSS(double valueIn, int orderIn, int nfmaxIn)
{
    init(valueIn, orderIn, nfmaxIn);
}


void Pythia8::AlphaStrong::init__BOSS(double valueIn, int orderIn)
{
    init(valueIn, orderIn);
}


void Pythia8::AlphaStrong::init__BOSS(double valueIn)
{
    init(valueIn);
}


void Pythia8::AlphaStrong::init__BOSS()
{
    init();
}




#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_AlphaStrong* Pythia8::AlphaStrong::pointerCopy__BOSS()
{
    Pythia8::Abstract_AlphaStrong* new_ptr = new Pythia8::AlphaStrong(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::AlphaStrong::pointerAssign__BOSS(Pythia8::Abstract_AlphaStrong* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::AlphaStrong* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<AlphaStrong*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
