#include "backend_types/BOSSedPythia_1_0/abstract_Vec4.h"
#include <vector>
#include <string>
#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"
#include "Pythia8/Event.h"

void Pythia8::Particle::setEvtPtr_GAMBIT(Pythia8::Abstract_Event* evtPtrIn)
{
    setEvtPtr(dynamic_cast< Pythia8::Event* >(evtPtrIn));
}


Pythia8::Abstract_Event* Pythia8::Particle::getEvtPtr_GAMBIT()
{
    return getEvtPtr();
}


void Pythia8::Particle::mothers_GAMBIT(int mother1In)
{
    mothers(mother1In);
}


void Pythia8::Particle::mothers_GAMBIT()
{
    mothers();
}


void Pythia8::Particle::daughters_GAMBIT(int daughter1In)
{
    daughters(daughter1In);
}


void Pythia8::Particle::daughters_GAMBIT()
{
    daughters();
}


void Pythia8::Particle::cols_GAMBIT(int colIn)
{
    cols(colIn);
}


void Pythia8::Particle::cols_GAMBIT()
{
    cols();
}


void Pythia8::Particle::p_GAMBIT(Pythia8::Abstract_Vec4& pIn)
{
    p(dynamic_cast< Pythia8::Vec4& >(pIn));
}


void Pythia8::Particle::vProd_GAMBIT(Pythia8::Abstract_Vec4& vProdIn)
{
    vProd(dynamic_cast< Pythia8::Vec4& >(vProdIn));
}


Pythia8::Abstract_Vec4* Pythia8::Particle::p_GAMBIT() const
{
    return new Pythia8::Vec4(p());
}


Pythia8::Abstract_Vec4* Pythia8::Particle::vProd_GAMBIT() const
{
    return new Pythia8::Vec4(vProd());
}


Pythia8::Abstract_Vec4* Pythia8::Particle::vDec_GAMBIT() const
{
    return new Pythia8::Vec4(vDec());
}


std::vector<int, std::allocator<int> > Pythia8::Particle::sisterList_GAMBIT() const
{
    return sisterList();
}


std::string Pythia8::Particle::nameWithStatus_GAMBIT() const
{
    return nameWithStatus();
}


void Pythia8::Particle::bst_GAMBIT(const Pythia8::Abstract_Vec4& pBst)
{
    bst(dynamic_cast< const Pythia8::Vec4& >(pBst));
}


void Pythia8::Particle::bst_GAMBIT(const Pythia8::Abstract_Vec4& pBst, double mBst)
{
    bst(dynamic_cast< const Pythia8::Vec4& >(pBst), mBst);
}


void Pythia8::Particle::bstback_GAMBIT(const Pythia8::Abstract_Vec4& pBst)
{
    bstback(dynamic_cast< const Pythia8::Vec4& >(pBst));
}


void Pythia8::Particle::bstback_GAMBIT(const Pythia8::Abstract_Vec4& pBst, double mBst)
{
    bstback(dynamic_cast< const Pythia8::Vec4& >(pBst), mBst);
}



Pythia8::Abstract_Particle* Pythia8::Particle::operator_equal_GAMBIT(const Pythia8::Abstract_Particle& pt)
{
    return &(operator=(dynamic_cast< const Pythia8::Particle& >(pt)));
}



Pythia8::Abstract_Particle* Pythia8::Particle::pointerCopy_GAMBIT() { return new Pythia8::Particle(*this); }
void Pythia8::Particle::pointerAssign_GAMBIT(Pythia8::Abstract_Particle* in) { *this = *dynamic_cast<Particle*>(in); }
