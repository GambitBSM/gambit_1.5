#include <string>
#include <ostream>
#include <vector>
#include "backend_types/Pythia_8_186/wrapper_Info.h"
#include "backend_types/Pythia_8_186/wrapper_Settings.h"
#include "backend_types/Pythia_8_186/wrapper_Rndm.h"
#include "backend_types/Pythia_8_186/wrapper_Couplings.h"
#include "backend_types/Pythia_8_186/wrapper_ParticleDataEntry.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/ParticleData.h"

void Pythia8::ParticleData::initPtr__BOSS(Pythia8::Abstract_Info* infoPtrIn, Pythia8::Abstract_Settings* settingsPtrIn, Pythia8::Abstract_Rndm* rndmPtrIn, Pythia8::Abstract_Couplings* couplingsPtrIn)
{
    initPtr(dynamic_cast< Pythia8::Info* >(infoPtrIn), dynamic_cast< Pythia8::Settings* >(settingsPtrIn), dynamic_cast< Pythia8::Rndm* >(rndmPtrIn), dynamic_cast< Pythia8::Couplings* >(couplingsPtrIn));
}


bool Pythia8::ParticleData::init__BOSS()
{
    return init();
}


bool Pythia8::ParticleData::reInit__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > startFile)
{
    return reInit(startFile);
}


bool Pythia8::ParticleData::readXML__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > inFile)
{
    return readXML(inFile);
}


bool Pythia8::ParticleData::readFF__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > inFile)
{
    return readFF(inFile);
}


bool Pythia8::ParticleData::readString__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > lineIn, bool warn)
{
    return readString(lineIn, warn);
}


bool Pythia8::ParticleData::readString__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > lineIn)
{
    return readString(lineIn);
}


void Pythia8::ParticleData::listAll__BOSS()
{
    listAll();
}


void Pythia8::ParticleData::listChanged__BOSS()
{
    listChanged();
}


void Pythia8::ParticleData::listChanged__BOSS(bool changedRes)
{
    listChanged(changedRes);
}


void Pythia8::ParticleData::list__BOSS(bool changedOnly, bool changedRes)
{
    list(changedOnly, changedRes);
}


void Pythia8::ParticleData::list__BOSS(bool changedOnly)
{
    list(changedOnly);
}


void Pythia8::ParticleData::list__BOSS()
{
    list();
}


void Pythia8::ParticleData::list__BOSS(int idList)
{
    list(idList);
}


void Pythia8::ParticleData::list__BOSS(std::vector<int,std::allocator<int> > idList)
{
    list(idList);
}


void Pythia8::ParticleData::checkTable__BOSS()
{
    checkTable();
}


void Pythia8::ParticleData::checkTable__BOSS(int verbosity)
{
    checkTable(verbosity);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn)
{
    addParticle(idIn, nameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn, mMaxIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn)
{
    addParticle(idIn, nameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn)
{
    addParticle(idIn, nameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In)
{
    addParticle(idIn, nameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn)
{
    addParticle(idIn, nameIn, spinTypeIn, chargeTypeIn, colTypeIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn)
{
    addParticle(idIn, nameIn, spinTypeIn, chargeTypeIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn)
{
    addParticle(idIn, nameIn, spinTypeIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn)
{
    addParticle(idIn, nameIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn)
{
    addParticle(idIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn)
{
    addParticle(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn, mMaxIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn)
{
    addParticle(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn)
{
    addParticle(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In)
{
    addParticle(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn)
{
    addParticle(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn)
{
    addParticle(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn)
{
    addParticle(idIn, nameIn, antiNameIn, spinTypeIn);
}


void Pythia8::ParticleData::addParticle__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn)
{
    addParticle(idIn, nameIn, antiNameIn);
}


void Pythia8::ParticleData::setAll__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn)
{
    setAll(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn, mMaxIn);
}


void Pythia8::ParticleData::setAll__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn)
{
    setAll(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn);
}


void Pythia8::ParticleData::setAll__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn)
{
    setAll(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn);
}


void Pythia8::ParticleData::setAll__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In)
{
    setAll(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In);
}


void Pythia8::ParticleData::setAll__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn)
{
    setAll(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn);
}


void Pythia8::ParticleData::setAll__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn)
{
    setAll(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn);
}


void Pythia8::ParticleData::setAll__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn)
{
    setAll(idIn, nameIn, antiNameIn, spinTypeIn);
}


void Pythia8::ParticleData::setAll__BOSS(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn)
{
    setAll(idIn, nameIn, antiNameIn);
}


void Pythia8::ParticleData::rescaleBR__BOSS(int idIn)
{
    rescaleBR(idIn);
}


double Pythia8::ParticleData::resWidth__BOSS(int idIn, double mHat, int idInFlav, bool openOnly)
{
    return resWidth(idIn, mHat, idInFlav, openOnly);
}


double Pythia8::ParticleData::resWidth__BOSS(int idIn, double mHat, int idInFlav)
{
    return resWidth(idIn, mHat, idInFlav);
}


double Pythia8::ParticleData::resWidth__BOSS(int idIn, double mHat)
{
    return resWidth(idIn, mHat);
}


double Pythia8::ParticleData::resWidthOpen__BOSS(int idIn, double mHat)
{
    return resWidthOpen(idIn, mHat);
}


double Pythia8::ParticleData::resWidthStore__BOSS(int idIn, double mHat)
{
    return resWidthStore(idIn, mHat);
}


double Pythia8::ParticleData::resOpenFrac__BOSS(int id1In, int id2In)
{
    return resOpenFrac(id1In, id2In);
}


double Pythia8::ParticleData::resOpenFrac__BOSS(int id1In)
{
    return resOpenFrac(id1In);
}


double Pythia8::ParticleData::resWidthChan__BOSS(int idIn, double mHat, int idAbs1)
{
    return resWidthChan(idIn, mHat, idAbs1);
}


double Pythia8::ParticleData::resWidthChan__BOSS(int idIn, double mHat)
{
    return resWidthChan(idIn, mHat);
}


Pythia8::Abstract_ParticleDataEntry* Pythia8::ParticleData::particleDataEntryPtr__BOSS(int idIn)
{
    return particleDataEntryPtr(idIn);
}




#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_ParticleData* Pythia8::ParticleData::pointerCopy__BOSS()
{
    Pythia8::Abstract_ParticleData* new_ptr = new Pythia8::ParticleData(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::ParticleData::pointerAssign__BOSS(Pythia8::Abstract_ParticleData* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::ParticleData* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<ParticleData*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
