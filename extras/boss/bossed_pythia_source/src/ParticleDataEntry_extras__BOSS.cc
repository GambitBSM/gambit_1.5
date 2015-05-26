#include <string>
#include "Pythia8/ResonanceWidths.h"
#include "Pythia8/ParticleData.h"
#include "backend_types/Pythia_8_209/wrapper_DecayChannel.h"
#include "backend_types/Pythia_8_209/wrapper_Info.h"
#include "backend_types/Pythia_8_209/wrapper_Settings.h"
#include "backend_types/Pythia_8_209/wrapper_Couplings.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

void Pythia8::ParticleDataEntry::initPtr__BOSS(Pythia8::Abstract_ParticleData* particleDataPtrIn)
{
    initPtr(dynamic_cast< Pythia8::ParticleData* >(particleDataPtrIn));
}


void Pythia8::ParticleDataEntry::setAll__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn)
{
    setAll(nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn, mMaxIn);
}


void Pythia8::ParticleDataEntry::setAll__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn)
{
    setAll(nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn);
}


void Pythia8::ParticleDataEntry::setAll__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn)
{
    setAll(nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn);
}


void Pythia8::ParticleDataEntry::setAll__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In)
{
    setAll(nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In);
}


void Pythia8::ParticleDataEntry::setAll__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn)
{
    setAll(nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn);
}


void Pythia8::ParticleDataEntry::setAll__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn)
{
    setAll(nameIn, antiNameIn, spinTypeIn, chargeTypeIn);
}


void Pythia8::ParticleDataEntry::setAll__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn)
{
    setAll(nameIn, antiNameIn, spinTypeIn);
}


void Pythia8::ParticleDataEntry::setAll__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn)
{
    setAll(nameIn, antiNameIn);
}


void Pythia8::ParticleDataEntry::setMWidth__BOSS(double mWidthIn)
{
    setMWidth(mWidthIn);
}


void Pythia8::ParticleDataEntry::setMayDecay__BOSS(bool mayDecayIn)
{
    setMayDecay(mayDecayIn);
}


std::basic_string<char,std::char_traits<char>,std::allocator<char> > Pythia8::ParticleDataEntry::name__BOSS() const
{
    return name();
}


int Pythia8::ParticleDataEntry::chargeType__BOSS() const
{
    return chargeType();
}


double Pythia8::ParticleDataEntry::charge__BOSS() const
{
    return charge();
}


int Pythia8::ParticleDataEntry::colType__BOSS() const
{
    return colType();
}


int Pythia8::ParticleDataEntry::heaviestQuark__BOSS() const
{
    return heaviestQuark();
}


int Pythia8::ParticleDataEntry::baryonNumberType__BOSS() const
{
    return baryonNumberType();
}


void Pythia8::ParticleDataEntry::addChannel__BOSS(int onMode, double bRatio, int meMode, int prod0, int prod1, int prod2, int prod3, int prod4, int prod5, int prod6)
{
    addChannel(onMode, bRatio, meMode, prod0, prod1, prod2, prod3, prod4, prod5, prod6);
}


void Pythia8::ParticleDataEntry::addChannel__BOSS(int onMode, double bRatio, int meMode, int prod0, int prod1, int prod2, int prod3, int prod4, int prod5)
{
    addChannel(onMode, bRatio, meMode, prod0, prod1, prod2, prod3, prod4, prod5);
}


void Pythia8::ParticleDataEntry::addChannel__BOSS(int onMode, double bRatio, int meMode, int prod0, int prod1, int prod2, int prod3, int prod4)
{
    addChannel(onMode, bRatio, meMode, prod0, prod1, prod2, prod3, prod4);
}


void Pythia8::ParticleDataEntry::addChannel__BOSS(int onMode, double bRatio, int meMode, int prod0, int prod1, int prod2, int prod3)
{
    addChannel(onMode, bRatio, meMode, prod0, prod1, prod2, prod3);
}


void Pythia8::ParticleDataEntry::addChannel__BOSS(int onMode, double bRatio, int meMode, int prod0, int prod1, int prod2)
{
    addChannel(onMode, bRatio, meMode, prod0, prod1, prod2);
}


void Pythia8::ParticleDataEntry::addChannel__BOSS(int onMode, double bRatio, int meMode, int prod0, int prod1)
{
    addChannel(onMode, bRatio, meMode, prod0, prod1);
}


void Pythia8::ParticleDataEntry::addChannel__BOSS(int onMode, double bRatio, int meMode, int prod0)
{
    addChannel(onMode, bRatio, meMode, prod0);
}


void Pythia8::ParticleDataEntry::addChannel__BOSS(int onMode, double bRatio, int meMode)
{
    addChannel(onMode, bRatio, meMode);
}


void Pythia8::ParticleDataEntry::addChannel__BOSS(int onMode, double bRatio)
{
    addChannel(onMode, bRatio);
}


void Pythia8::ParticleDataEntry::addChannel__BOSS(int onMode)
{
    addChannel(onMode);
}


void Pythia8::ParticleDataEntry::addChannel__BOSS()
{
    addChannel();
}


Pythia8::Abstract_DecayChannel* Pythia8::ParticleDataEntry::channel__BOSS(int i)
{
    return &(channel(i));
}


const Pythia8::Abstract_DecayChannel* Pythia8::ParticleDataEntry::channel__BOSS(int i) const
{
    return &(channel(i));
}


void Pythia8::ParticleDataEntry::rescaleBR__BOSS()
{
    rescaleBR();
}


bool Pythia8::ParticleDataEntry::preparePick__BOSS(int idSgn, double mHat)
{
    return preparePick(idSgn, mHat);
}


bool Pythia8::ParticleDataEntry::preparePick__BOSS(int idSgn)
{
    return preparePick(idSgn);
}


Pythia8::Abstract_DecayChannel* Pythia8::ParticleDataEntry::pickChannel__BOSS()
{
    return &(pickChannel());
}


void Pythia8::ParticleDataEntry::setResonancePtr__BOSS(Pythia8::Abstract_ResonanceWidths* resonancePtrIn)
{
    setResonancePtr(dynamic_cast< Pythia8::ResonanceWidths* >(resonancePtrIn));
}


Pythia8::Abstract_ResonanceWidths* Pythia8::ParticleDataEntry::getResonancePtr__BOSS()
{
    return getResonancePtr();
}


void Pythia8::ParticleDataEntry::resInit__BOSS(Pythia8::Abstract_Info* infoPtrIn, Pythia8::Abstract_Settings* settingsPtrIn, Pythia8::Abstract_ParticleData* particleDataPtrIn, Pythia8::Abstract_Couplings* couplingsPtrIn)
{
    resInit(dynamic_cast< Pythia8::Info* >(infoPtrIn), dynamic_cast< Pythia8::Settings* >(settingsPtrIn), dynamic_cast< Pythia8::ParticleData* >(particleDataPtrIn), dynamic_cast< Pythia8::Couplings* >(couplingsPtrIn));
}


double Pythia8::ParticleDataEntry::resWidth__BOSS(int idSgn, double mHat, int idIn, bool openOnly)
{
    return resWidth(idSgn, mHat, idIn, openOnly);
}


double Pythia8::ParticleDataEntry::resWidth__BOSS(int idSgn, double mHat, int idIn)
{
    return resWidth(idSgn, mHat, idIn);
}


double Pythia8::ParticleDataEntry::resWidth__BOSS(int idSgn, double mHat)
{
    return resWidth(idSgn, mHat);
}


double Pythia8::ParticleDataEntry::resWidthOpen__BOSS(int idSgn, double mHat)
{
    return resWidthOpen(idSgn, mHat);
}


double Pythia8::ParticleDataEntry::resWidthStore__BOSS(int idSgn, double mHat)
{
    return resWidthStore(idSgn, mHat);
}


double Pythia8::ParticleDataEntry::resWidthChan__BOSS(double mHat, int idAbs1)
{
    return resWidthChan(mHat, idAbs1);
}


double Pythia8::ParticleDataEntry::resWidthChan__BOSS(double mHat)
{
    return resWidthChan(mHat);
}




#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_ParticleDataEntry* Pythia8::ParticleDataEntry::pointerCopy__BOSS()
{
    Pythia8::Abstract_ParticleDataEntry* new_ptr = new Pythia8::ParticleDataEntry(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::ParticleDataEntry::pointerAssign__BOSS(Pythia8::Abstract_ParticleDataEntry* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::ParticleDataEntry* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<ParticleDataEntry*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
