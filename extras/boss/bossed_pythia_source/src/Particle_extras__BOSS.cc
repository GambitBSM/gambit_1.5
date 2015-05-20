#include "Pythia8/Event.h"
#include <vector>
#include <string>
#include "backend_types/Pythia_8_186/wrapper_Vec4.h"
#include "backend_types/Pythia_8_186/wrapper_ParticleDataEntry.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

void Pythia8::Particle::setEvtPtr__BOSS(Pythia8::Abstract_Event* evtPtrIn)
{
    setEvtPtr(dynamic_cast< Pythia8::Event* >(evtPtrIn));
}


void Pythia8::Particle::setPDEPtr__BOSS(Pythia8::Abstract_ParticleDataEntry* pdePtrIn)
{
    setPDEPtr(dynamic_cast< Pythia8::ParticleDataEntry* >(pdePtrIn));
}


void Pythia8::Particle::setPDEPtr__BOSS()
{
    setPDEPtr();
}


void Pythia8::Particle::mothers__BOSS(int mother1In)
{
    mothers(mother1In);
}


void Pythia8::Particle::mothers__BOSS()
{
    mothers();
}


void Pythia8::Particle::daughters__BOSS(int daughter1In)
{
    daughters(daughter1In);
}


void Pythia8::Particle::daughters__BOSS()
{
    daughters();
}


void Pythia8::Particle::cols__BOSS(int colIn)
{
    cols(colIn);
}


void Pythia8::Particle::cols__BOSS()
{
    cols();
}


void Pythia8::Particle::p__BOSS(Pythia8::Abstract_Vec4& pIn)
{
    p(dynamic_cast< Pythia8::Vec4& >(pIn));
}


void Pythia8::Particle::vProd__BOSS(Pythia8::Abstract_Vec4& vProdIn)
{
    vProd(dynamic_cast< Pythia8::Vec4& >(vProdIn));
}


Pythia8::Abstract_Vec4* Pythia8::Particle::p__BOSS() const
{
    return new Pythia8::Vec4(p());
}


Pythia8::Abstract_Vec4* Pythia8::Particle::vProd__BOSS() const
{
    return new Pythia8::Vec4(vProd());
}


Pythia8::Abstract_Vec4* Pythia8::Particle::vDec__BOSS() const
{
    return new Pythia8::Vec4(vDec());
}


std::vector<int,std::allocator<int> > Pythia8::Particle::sisterList__BOSS() const
{
    return sisterList();
}


std::basic_string<char,std::char_traits<char>,std::allocator<char> > Pythia8::Particle::nameWithStatus__BOSS() const
{
    return nameWithStatus();
}


Pythia8::Abstract_ParticleDataEntry* Pythia8::Particle::particleDataEntry__BOSS() const
{
    return &(particleDataEntry());
}


void Pythia8::Particle::bst__BOSS(const Pythia8::Abstract_Vec4& pBst)
{
    bst(dynamic_cast< const Pythia8::Vec4& >(pBst));
}


void Pythia8::Particle::bst__BOSS(const Pythia8::Abstract_Vec4& pBst, double mBst)
{
    bst(dynamic_cast< const Pythia8::Vec4& >(pBst), mBst);
}


void Pythia8::Particle::bstback__BOSS(const Pythia8::Abstract_Vec4& pBst)
{
    bstback(dynamic_cast< const Pythia8::Vec4& >(pBst));
}


void Pythia8::Particle::bstback__BOSS(const Pythia8::Abstract_Vec4& pBst, double mBst)
{
    bstback(dynamic_cast< const Pythia8::Vec4& >(pBst), mBst);
}



Pythia8::Abstract_Particle* Pythia8::Particle::operator_equal__BOSS(const Pythia8::Abstract_Particle& pt)
{
    return &(operator=(dynamic_cast< const Pythia8::Particle& >(pt)));
}



#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_Particle* Pythia8::Particle::pointerCopy__BOSS()
{
    Pythia8::Abstract_Particle* new_ptr = new Pythia8::Particle(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::Particle::pointerAssign__BOSS(Pythia8::Abstract_Particle* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Particle* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<Particle*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
