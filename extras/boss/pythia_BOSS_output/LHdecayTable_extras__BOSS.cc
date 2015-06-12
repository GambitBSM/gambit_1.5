#include <vector>
#include <string>
#include "backend_types/Pythia_8_209/wrapper_LHdecayChannel.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/SusyLesHouches.h"

void Pythia8::LHdecayTable::reset__BOSS()
{
    reset();
}


void Pythia8::LHdecayTable::addChannel__BOSS(Pythia8::Abstract_LHdecayChannel& channelIn)
{
    addChannel(dynamic_cast< Pythia8::LHdecayChannel& >(channelIn));
}


void Pythia8::LHdecayTable::addChannel__BOSS(double bratIn, int nDaIn, std::vector<int,std::allocator<int> > idDaIn)
{
    addChannel(bratIn, nDaIn, idDaIn);
}


Pythia8::Abstract_LHdecayChannel* Pythia8::LHdecayTable::getChannel__BOSS(int iChannel)
{
    return new Pythia8::LHdecayChannel(getChannel(iChannel));
}




#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_LHdecayTable* Pythia8::LHdecayTable::pointerCopy__BOSS()
{
    Pythia8::Abstract_LHdecayTable* new_ptr = new Pythia8::LHdecayTable(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::LHdecayTable::pointerAssign__BOSS(Pythia8::Abstract_LHdecayTable* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::LHdecayTable* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<LHdecayTable*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
