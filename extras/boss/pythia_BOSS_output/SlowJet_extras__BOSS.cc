#include <vector>
#include <ostream>
#include "backend_types/Pythia_8_209/wrapper_Event.h"
#include "backend_types/Pythia_8_209/wrapper_Vec4.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/Analysis.h"

bool Pythia8::SlowJet::analyze__BOSS(const Pythia8::Abstract_Event& event)
{
    return analyze(dynamic_cast< const Pythia8::Event& >(event));
}


bool Pythia8::SlowJet::setup__BOSS(const Pythia8::Abstract_Event& event)
{
    return setup(dynamic_cast< const Pythia8::Event& >(event));
}


Pythia8::Abstract_Vec4* Pythia8::SlowJet::p__BOSS(int i) const
{
    return new Pythia8::Vec4(p(i));
}


void Pythia8::SlowJet::list__BOSS(bool listAll) const
{
    list(listAll);
}


void Pythia8::SlowJet::list__BOSS() const
{
    list();
}




#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_SlowJet* Pythia8::SlowJet::pointerCopy__BOSS()
{
    Pythia8::Abstract_SlowJet* new_ptr = new Pythia8::SlowJet(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::SlowJet::pointerAssign__BOSS(Pythia8::Abstract_SlowJet* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::SlowJet* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<SlowJet*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
