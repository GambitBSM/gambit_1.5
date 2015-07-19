#include <vector>
#include <string>
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/SusyLesHouches.h"

void Pythia8::LHdecayChannel::setChannel__BOSS(double bratIn, int nDaIn, std::vector<int,std::allocator<int> > idDaIn)
{
    setChannel(bratIn, nDaIn, idDaIn);
}




#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_LHdecayChannel* Pythia8::LHdecayChannel::pointerCopy__BOSS()
{
    Pythia8::Abstract_LHdecayChannel* new_ptr = new Pythia8::LHdecayChannel(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::LHdecayChannel::pointerAssign__BOSS(Pythia8::Abstract_LHdecayChannel* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::LHdecayChannel* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<LHdecayChannel*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
