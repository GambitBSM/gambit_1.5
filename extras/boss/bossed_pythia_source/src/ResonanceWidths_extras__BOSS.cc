#include "Pythia8/ParticleData.h"
#include "backend_types/Pythia_8_209/wrapper_Info.h"
#include "backend_types/Pythia_8_209/wrapper_Settings.h"
#include "backend_types/Pythia_8_209/wrapper_Couplings.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/ResonanceWidths.h"

void Pythia8::ResonanceWidths::initBasic__BOSS(int idResIn)
{
    initBasic(idResIn);
}


bool Pythia8::ResonanceWidths::init__BOSS(Pythia8::Abstract_Info* infoPtrIn, Pythia8::Abstract_Settings* settingsPtrIn, Pythia8::Abstract_ParticleData* particleDataPtrIn, Pythia8::Abstract_Couplings* couplingsPtrIn)
{
    return init(dynamic_cast< Pythia8::Info* >(infoPtrIn), dynamic_cast< Pythia8::Settings* >(settingsPtrIn), dynamic_cast< Pythia8::ParticleData* >(particleDataPtrIn), dynamic_cast< Pythia8::Couplings* >(couplingsPtrIn));
}


double Pythia8::ResonanceWidths::width__BOSS(int idSgn, double mHatIn, int idInFlavIn, bool openOnly, bool setBR, int idOutFlav1)
{
    return width(idSgn, mHatIn, idInFlavIn, openOnly, setBR, idOutFlav1);
}


double Pythia8::ResonanceWidths::width__BOSS(int idSgn, double mHatIn, int idInFlavIn, bool openOnly, bool setBR)
{
    return width(idSgn, mHatIn, idInFlavIn, openOnly, setBR);
}


double Pythia8::ResonanceWidths::width__BOSS(int idSgn, double mHatIn, int idInFlavIn, bool openOnly)
{
    return width(idSgn, mHatIn, idInFlavIn, openOnly);
}


double Pythia8::ResonanceWidths::width__BOSS(int idSgn, double mHatIn, int idInFlavIn)
{
    return width(idSgn, mHatIn, idInFlavIn);
}


double Pythia8::ResonanceWidths::width__BOSS(int idSgn, double mHatIn)
{
    return width(idSgn, mHatIn);
}


double Pythia8::ResonanceWidths::widthOpen__BOSS(int idSgn, double mHatIn)
{
    return widthOpen(idSgn, mHatIn);
}


double Pythia8::ResonanceWidths::widthStore__BOSS(int idSgn, double mHatIn)
{
    return widthStore(idSgn, mHatIn);
}


void Pythia8::ResonanceWidths::calcPreFac__BOSS()
{
    calcPreFac();
}


void Pythia8::ResonanceWidths::calcWidth__BOSS()
{
    calcWidth();
}


double Pythia8::ResonanceWidths::numInt1BW__BOSS(double mHatIn, double m1, double Gamma1, double mMin1, double m2)
{
    return numInt1BW(mHatIn, m1, Gamma1, mMin1, m2);
}


double Pythia8::ResonanceWidths::numInt2BW__BOSS(double mHatIn, double m1, double Gamma1, double mMin1, double m2, double Gamma2, double mMin2)
{
    return numInt2BW(mHatIn, m1, Gamma1, mMin1, m2, Gamma2, mMin2);
}




#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_ResonanceWidths* Pythia8::ResonanceWidths::pointerCopy__BOSS()
{
    Pythia8::Abstract_ResonanceWidths* new_ptr = new Pythia8::ResonanceWidths(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::ResonanceWidths::pointerAssign__BOSS(Pythia8::Abstract_ResonanceWidths* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::ResonanceWidths* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<ResonanceWidths*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
