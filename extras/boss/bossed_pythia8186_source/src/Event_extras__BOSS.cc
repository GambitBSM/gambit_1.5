#include <string>
#include "Pythia8/Event.h"
#include <ostream>
#include <vector>
#include "backend_types/Pythia_8_186/wrapper_ParticleData.h"
#include "backend_types/Pythia_8_186/wrapper_Particle.h"
#include "backend_types/Pythia_8_186/wrapper_Vec4.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

void Pythia8::Event::init__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > headerIn, Pythia8::Abstract_ParticleData* particleDataPtrIn, int startColTagIn)
{
    init(headerIn, dynamic_cast< Pythia8::ParticleData* >(particleDataPtrIn), startColTagIn);
}


void Pythia8::Event::init__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > headerIn, Pythia8::Abstract_ParticleData* particleDataPtrIn)
{
    init(headerIn, dynamic_cast< Pythia8::ParticleData* >(particleDataPtrIn));
}


void Pythia8::Event::init__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > headerIn)
{
    init(headerIn);
}


void Pythia8::Event::init__BOSS()
{
    init();
}


Pythia8::Abstract_Particle* Pythia8::Event::front__BOSS()
{
    return &(front());
}


Pythia8::Abstract_Particle* Pythia8::Event::at__BOSS(int i)
{
    return &(at(i));
}


Pythia8::Abstract_Particle* Pythia8::Event::back__BOSS()
{
    return &(back());
}


int Pythia8::Event::append__BOSS(Pythia8::Abstract_Particle& entryIn)
{
    return append(dynamic_cast< Pythia8::Particle& >(entryIn));
}


int Pythia8::Event::append__BOSS(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, px, py, pz, e, m, scaleIn);
}


int Pythia8::Event::append__BOSS(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz, double e, double m)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, px, py, pz, e, m);
}


int Pythia8::Event::append__BOSS(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz, double e)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, px, py, pz, e);
}


int Pythia8::Event::append__BOSS(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract_Vec4& p, double m, double scaleIn, double polIn)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, dynamic_cast< Pythia8::Vec4& >(p), m, scaleIn, polIn);
}


int Pythia8::Event::append__BOSS(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract_Vec4& p, double m, double scaleIn)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, dynamic_cast< Pythia8::Vec4& >(p), m, scaleIn);
}


int Pythia8::Event::append__BOSS(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract_Vec4& p, double m)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, dynamic_cast< Pythia8::Vec4& >(p), m);
}


int Pythia8::Event::append__BOSS(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract_Vec4& p)
{
    return append(id, status, mother1, mother2, daughter1, daughter2, col, acol, dynamic_cast< Pythia8::Vec4& >(p));
}


int Pythia8::Event::append__BOSS(int id, int status, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn)
{
    return append(id, status, col, acol, px, py, pz, e, m, scaleIn);
}


int Pythia8::Event::append__BOSS(int id, int status, int col, int acol, double px, double py, double pz, double e, double m)
{
    return append(id, status, col, acol, px, py, pz, e, m);
}


int Pythia8::Event::append__BOSS(int id, int status, int col, int acol, double px, double py, double pz, double e)
{
    return append(id, status, col, acol, px, py, pz, e);
}


int Pythia8::Event::append__BOSS(int id, int status, int col, int acol, Pythia8::Abstract_Vec4& p, double m, double scaleIn, double polIn)
{
    return append(id, status, col, acol, dynamic_cast< Pythia8::Vec4& >(p), m, scaleIn, polIn);
}


int Pythia8::Event::append__BOSS(int id, int status, int col, int acol, Pythia8::Abstract_Vec4& p, double m, double scaleIn)
{
    return append(id, status, col, acol, dynamic_cast< Pythia8::Vec4& >(p), m, scaleIn);
}


int Pythia8::Event::append__BOSS(int id, int status, int col, int acol, Pythia8::Abstract_Vec4& p, double m)
{
    return append(id, status, col, acol, dynamic_cast< Pythia8::Vec4& >(p), m);
}


int Pythia8::Event::append__BOSS(int id, int status, int col, int acol, Pythia8::Abstract_Vec4& p)
{
    return append(id, status, col, acol, dynamic_cast< Pythia8::Vec4& >(p));
}


void Pythia8::Event::setEvtPtr__BOSS()
{
    setEvtPtr();
}


int Pythia8::Event::copy__BOSS(int iCopy)
{
    return copy(iCopy);
}


void Pythia8::Event::list__BOSS(bool showScaleAndVertex) const
{
    list(showScaleAndVertex);
}


void Pythia8::Event::popBack__BOSS()
{
    popBack();
}


void Pythia8::Event::initColTag__BOSS()
{
    initColTag();
}


std::vector<int,std::allocator<int> > Pythia8::Event::sisterListTopBot__BOSS(int i) const
{
    return sisterListTopBot(i);
}


void Pythia8::Event::bst__BOSS(const Pythia8::Abstract_Vec4& vec)
{
    bst(dynamic_cast< const Pythia8::Vec4& >(vec));
}


void Pythia8::Event::listJunctions__BOSS() const
{
    listJunctions();
}



Pythia8::Abstract_Event* Pythia8::Event::operator_equal__BOSS(const Pythia8::Abstract_Event& oldEvent)
{
    return &(operator=(dynamic_cast< const Pythia8::Event& >(oldEvent)));
}


Pythia8::Abstract_Particle* Pythia8::Event::operator_square_bracket_pair__BOSS(int i)
{
    return &(operator[](i));
}


const Pythia8::Abstract_Particle* Pythia8::Event::operator_square_bracket_pair__BOSS(int i) const
{
    return &(operator[](i));
}


Pythia8::Abstract_Event* Pythia8::Event::operator_plus_equal__BOSS(const Pythia8::Abstract_Event& addEvent)
{
    return &(operator+=(dynamic_cast< const Pythia8::Event& >(addEvent)));
}



#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_Event* Pythia8::Event::pointerCopy__BOSS()
{
    Pythia8::Abstract_Event* new_ptr = new Pythia8::Event(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::Event::pointerAssign__BOSS(Pythia8::Abstract_Event* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Event* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<Event*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
