#include <ostream>
#include <vector>
#include "backend_types/Pythia_8_186/wrapper_Info.h"
#include "backend_types/Pythia_8_186/wrapper_ParticleData.h"
#include "backend_types/Pythia_8_186/wrapper_Rndm.h"
#include "backend_types/Pythia_8_186/wrapper_Vec4.h"
#include "backend_types/Pythia_8_186/wrapper_Settings.h"
#include "backend_types/Pythia_8_186/wrapper_Event.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/BeamParticle.h"

Pythia8::Abstract_Vec4* Pythia8::BeamParticle::p__BOSS() const
{
    return new Pythia8::Vec4(p());
}


double Pythia8::BeamParticle::xMax__BOSS()
{
    return xMax();
}


int Pythia8::BeamParticle::append__BOSS(int iPos, int idIn, double x)
{
    return append(iPos, idIn, x);
}


void Pythia8::BeamParticle::list__BOSS() const
{
    list();
}


bool Pythia8::BeamParticle::remnantFlavours__BOSS(Pythia8::Abstract_Event& event)
{
    return remnantFlavours(dynamic_cast< Pythia8::Event& >(event));
}


bool Pythia8::BeamParticle::remnantColours__BOSS(Pythia8::Abstract_Event& event, std::vector<int,std::allocator<int> >& colFrom, std::vector<int,std::allocator<int> >& colTo)
{
    return remnantColours(dynamic_cast< Pythia8::Event& >(event), colFrom, colTo);
}



Pythia8::ResolvedParton& Pythia8::BeamParticle::operator_square_bracket_pair__BOSS(int i)
{
    return operator[](i);
}


const Pythia8::ResolvedParton& Pythia8::BeamParticle::operator_square_bracket_pair__BOSS(int i) const
{
    return operator[](i);
}



#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_BeamParticle* Pythia8::BeamParticle::pointerCopy__BOSS()
{
    Pythia8::Abstract_BeamParticle* new_ptr = new Pythia8::BeamParticle(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::BeamParticle::pointerAssign__BOSS(Pythia8::Abstract_BeamParticle* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::BeamParticle* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<BeamParticle*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
