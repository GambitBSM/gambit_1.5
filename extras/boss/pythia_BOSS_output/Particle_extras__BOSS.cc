#include "backend_types/BOSSedPythia_1_0/abstract_Vec4.h"
#include <vector>
#include <string>
#include "abstracttypedefs.h"
#include "wrappertypedefs.h"
#include "Pythia8/Event.h"

void Pythia8::Particle::setEvtPtr__BOSS(Pythia8::Abstract_Event* evtPtrIn)
{
    setEvtPtr(dynamic_cast< Pythia8::Event* >(evtPtrIn));
}


Pythia8::Abstract_Event* Pythia8::Particle::getEvtPtr__BOSS()
{
    return getEvtPtr();
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


std::vector<int, std::allocator<int> > Pythia8::Particle::sisterList__BOSS() const
{
    return sisterList();
}


std::string Pythia8::Particle::nameWithStatus__BOSS() const
{
    return nameWithStatus();
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



Pythia8::Abstract_Particle* Pythia8::Particle::pointerCopy__BOSS() { return new Pythia8::Particle(*this); }
void Pythia8::Particle::pointerAssign__BOSS(Pythia8::Abstract_Particle* in) { *this = *dynamic_cast<Particle*>(in); }
