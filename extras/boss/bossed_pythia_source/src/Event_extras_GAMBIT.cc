#include <vector>
#include <string>
#include "backend_types/BOSSedPythia_1_0/abstract_Particle.h"
#include "backend_types/BOSSedPythia_1_0/abstract_Vec4.h"
#include <ostream>
#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"
#include "Pythia8/Event.h"

Pythia8::Abstract_Particle* Pythia8::Event::front_GAMBIT()
{
    return &(front());
}


Pythia8::Abstract_Particle* Pythia8::Event::at_GAMBIT(int i)
{
    return &(at(i));
}


Pythia8::Abstract_Particle* Pythia8::Event::back_GAMBIT()
{
    return &(back());
}


int Pythia8::Event::append_GAMBIT(Pythia8::Abstract_Particle& entryIn)
{
    return append(dynamic_cast< Pythia8::Particle& >(entryIn));
}


int Pythia8::Event::append_GAMBIT(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, px, py, pz, e, m, scaleIn);
}


int Pythia8::Event::append_GAMBIT(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz, double e, double m)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, px, py, pz, e, m);
}


int Pythia8::Event::append_GAMBIT(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz, double e)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, px, py, pz, e);
}


int Pythia8::Event::append_GAMBIT(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract_Vec4& p, double m, double scaleIn, double polIn)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, dynamic_cast< Pythia8::Vec4& >(p), m, scaleIn, polIn);
}


int Pythia8::Event::append_GAMBIT(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract_Vec4& p, double m, double scaleIn)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, dynamic_cast< Pythia8::Vec4& >(p), m, scaleIn);
}


int Pythia8::Event::append_GAMBIT(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract_Vec4& p, double m)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, dynamic_cast< Pythia8::Vec4& >(p), m);
}


int Pythia8::Event::append_GAMBIT(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract_Vec4& p)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, dynamic_cast< Pythia8::Vec4& >(p));
}


int Pythia8::Event::append_GAMBIT(int id, int status, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn)
{
    return append(id, status, col, acol, px, py, pz, e, m, scaleIn);
}


int Pythia8::Event::append_GAMBIT(int id, int status, int col, int acol, double px, double py, double pz, double e, double m)
{
    return append(id, status, col, acol, px, py, pz, e, m);
}


int Pythia8::Event::append_GAMBIT(int id, int status, int col, int acol, double px, double py, double pz, double e)
{
    return append(id, status, col, acol, px, py, pz, e);
}


int Pythia8::Event::append_GAMBIT(int id, int status, int col, int acol, Pythia8::Abstract_Vec4& p, double m, double scaleIn, double polIn)
{
    return append(id, status, col, acol, dynamic_cast< Pythia8::Vec4& >(p), m, scaleIn, polIn);
}


int Pythia8::Event::append_GAMBIT(int id, int status, int col, int acol, Pythia8::Abstract_Vec4& p, double m, double scaleIn)
{
    return append(id, status, col, acol, dynamic_cast< Pythia8::Vec4& >(p), m, scaleIn);
}


int Pythia8::Event::append_GAMBIT(int id, int status, int col, int acol, Pythia8::Abstract_Vec4& p, double m)
{
    return append(id, status, col, acol, dynamic_cast< Pythia8::Vec4& >(p), m);
}


int Pythia8::Event::append_GAMBIT(int id, int status, int col, int acol, Pythia8::Abstract_Vec4& p)
{
    return append(id, status, col, acol, dynamic_cast< Pythia8::Vec4& >(p));
}


void Pythia8::Event::setEvtPtr_GAMBIT()
{
    setEvtPtr();
}


int Pythia8::Event::copy_GAMBIT(int iCopy)
{
    return copy(iCopy);
}


void Pythia8::Event::list_GAMBIT(bool showScaleAndVertex) const
{
    list(showScaleAndVertex);
}


void Pythia8::Event::popBack_GAMBIT()
{
    popBack();
}


void Pythia8::Event::initColTag_GAMBIT()
{
    initColTag();
}


std::vector<int, std::allocator<int> > Pythia8::Event::sisterListTopBot_GAMBIT(int i) const
{
    return sisterListTopBot(i);
}


void Pythia8::Event::bst_GAMBIT(const Pythia8::Abstract_Vec4& vec)
{
    bst(dynamic_cast< const Pythia8::Vec4& >(vec));
}


void Pythia8::Event::listJunctions_GAMBIT() const
{
    listJunctions();
}



Pythia8::Abstract_Event* Pythia8::Event::operator_equal_GAMBIT(const Pythia8::Abstract_Event& oldEvent)
{
    return &(operator=(dynamic_cast< const Pythia8::Event& >(oldEvent)));
}


Pythia8::Abstract_Particle* Pythia8::Event::operator_square_bracket_pair_GAMBIT(int i)
{
    return &(operator[](i));
}


const Pythia8::Abstract_Particle* Pythia8::Event::operator_square_bracket_pair_GAMBIT(int i) const
{
    return &(operator[](i));
}


Pythia8::Abstract_Event* Pythia8::Event::operator_plus_equal_GAMBIT(const Pythia8::Abstract_Event& addEvent)
{
    return &(operator+=(dynamic_cast< const Pythia8::Event& >(addEvent)));
}



Pythia8::Abstract_Event* Pythia8::Event::pointerCopy_GAMBIT() { return new Pythia8::Event(*this); }
void Pythia8::Event::pointerAssign_GAMBIT(Pythia8::Abstract_Event* in) { *this = *dynamic_cast<Event*>(in); }
