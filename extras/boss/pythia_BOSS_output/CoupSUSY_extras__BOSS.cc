#include <complex>
#include "Pythia8/ParticleData.h"
#include "backend_types/Pythia_8_209/wrapper_SusyLesHouches.h"
#include "backend_types/Pythia_8_209/wrapper_Info.h"
#include "backend_types/Pythia_8_209/wrapper_Settings.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/SusyCouplings.h"

void Pythia8::CoupSUSY::initSUSY__BOSS(Pythia8::Abstract_SusyLesHouches* slhaPtrIn, Pythia8::Abstract_Info* infoPtrIn, Pythia8::Abstract_ParticleData* particleDataPtrIn, Pythia8::Abstract_Settings* settingsPtrIn)
{
    initSUSY(dynamic_cast< Pythia8::SusyLesHouches* >(slhaPtrIn), dynamic_cast< Pythia8::Info* >(infoPtrIn), dynamic_cast< Pythia8::ParticleData* >(particleDataPtrIn), dynamic_cast< Pythia8::Settings* >(settingsPtrIn));
}



bool& Pythia8::CoupSUSY::isInit_ref__BOSS() { return isInit; }

bool& Pythia8::CoupSUSY::isNMSSM_ref__BOSS() { return isNMSSM; }

double& Pythia8::CoupSUSY::mWpole_ref__BOSS() { return mWpole; }

double& Pythia8::CoupSUSY::wWpole_ref__BOSS() { return wWpole; }

double& Pythia8::CoupSUSY::mZpole_ref__BOSS() { return mZpole; }

double& Pythia8::CoupSUSY::wZpole_ref__BOSS() { return wZpole; }

double& Pythia8::CoupSUSY::mW_ref__BOSS() { return mW; }

double& Pythia8::CoupSUSY::mZ_ref__BOSS() { return mZ; }

double& Pythia8::CoupSUSY::sin2W_ref__BOSS() { return sin2W; }

double& Pythia8::CoupSUSY::sinW_ref__BOSS() { return sinW; }

double& Pythia8::CoupSUSY::cosW_ref__BOSS() { return cosW; }

double& Pythia8::CoupSUSY::tanb_ref__BOSS() { return tanb; }

double& Pythia8::CoupSUSY::cosb_ref__BOSS() { return cosb; }

double& Pythia8::CoupSUSY::sinb_ref__BOSS() { return sinb; }

double& Pythia8::CoupSUSY::muHiggs_ref__BOSS() { return muHiggs; }

double& Pythia8::CoupSUSY::alphaHiggs_ref__BOSS() { return alphaHiggs; }

double& Pythia8::CoupSUSY::mAHiggs_ref__BOSS() { return mAHiggs; }

std::complex<double> (&Pythia8::CoupSUSY::LsddG_ref__BOSS())[7][4] { return LsddG; }

std::complex<double> (&Pythia8::CoupSUSY::RsddG_ref__BOSS())[7][4] { return RsddG; }

std::complex<double> (&Pythia8::CoupSUSY::LsuuG_ref__BOSS())[7][4] { return LsuuG; }

std::complex<double> (&Pythia8::CoupSUSY::RsuuG_ref__BOSS())[7][4] { return RsuuG; }

std::complex<double> (&Pythia8::CoupSUSY::OLpp_ref__BOSS())[6][6] { return OLpp; }

std::complex<double> (&Pythia8::CoupSUSY::ORpp_ref__BOSS())[6][6] { return ORpp; }

std::complex<double> (&Pythia8::CoupSUSY::OLp_ref__BOSS())[3][3] { return OLp; }

std::complex<double> (&Pythia8::CoupSUSY::ORp_ref__BOSS())[3][3] { return ORp; }

std::complex<double> (&Pythia8::CoupSUSY::OL_ref__BOSS())[6][3] { return OL; }

std::complex<double> (&Pythia8::CoupSUSY::OR_ref__BOSS())[6][3] { return OR; }

double (&Pythia8::CoupSUSY::LqqZ_ref__BOSS())[7] { return LqqZ; }

double (&Pythia8::CoupSUSY::RqqZ_ref__BOSS())[7] { return RqqZ; }

std::complex<double> (&Pythia8::CoupSUSY::LsdsdZ_ref__BOSS())[7][7] { return LsdsdZ; }

std::complex<double> (&Pythia8::CoupSUSY::RsdsdZ_ref__BOSS())[7][7] { return RsdsdZ; }

std::complex<double> (&Pythia8::CoupSUSY::LsusuZ_ref__BOSS())[7][7] { return LsusuZ; }

std::complex<double> (&Pythia8::CoupSUSY::RsusuZ_ref__BOSS())[7][7] { return RsusuZ; }

std::complex<double> (&Pythia8::CoupSUSY::LudW_ref__BOSS())[4][4] { return LudW; }

std::complex<double> (&Pythia8::CoupSUSY::RudW_ref__BOSS())[4][4] { return RudW; }

std::complex<double> (&Pythia8::CoupSUSY::LsusdW_ref__BOSS())[7][7] { return LsusdW; }

std::complex<double> (&Pythia8::CoupSUSY::RsusdW_ref__BOSS())[7][7] { return RsusdW; }

std::complex<double> (&Pythia8::CoupSUSY::LsddX_ref__BOSS())[7][4][6] { return LsddX; }

std::complex<double> (&Pythia8::CoupSUSY::RsddX_ref__BOSS())[7][4][6] { return RsddX; }

std::complex<double> (&Pythia8::CoupSUSY::LsuuX_ref__BOSS())[7][4][6] { return LsuuX; }

std::complex<double> (&Pythia8::CoupSUSY::RsuuX_ref__BOSS())[7][4][6] { return RsuuX; }

std::complex<double> (&Pythia8::CoupSUSY::LsduX_ref__BOSS())[7][4][3] { return LsduX; }

std::complex<double> (&Pythia8::CoupSUSY::RsduX_ref__BOSS())[7][4][3] { return RsduX; }

std::complex<double> (&Pythia8::CoupSUSY::LsudX_ref__BOSS())[7][4][3] { return LsudX; }

std::complex<double> (&Pythia8::CoupSUSY::RsudX_ref__BOSS())[7][4][3] { return RsudX; }

double (&Pythia8::CoupSUSY::LllZ_ref__BOSS())[7] { return LllZ; }

double (&Pythia8::CoupSUSY::RllZ_ref__BOSS())[7] { return RllZ; }

std::complex<double> (&Pythia8::CoupSUSY::LlvW_ref__BOSS())[4][4] { return LlvW; }

std::complex<double> (&Pythia8::CoupSUSY::RlvW_ref__BOSS())[4][4] { return RlvW; }

std::complex<double> (&Pythia8::CoupSUSY::LslslZ_ref__BOSS())[7][7] { return LslslZ; }

std::complex<double> (&Pythia8::CoupSUSY::RslslZ_ref__BOSS())[7][7] { return RslslZ; }

std::complex<double> (&Pythia8::CoupSUSY::LsvsvZ_ref__BOSS())[7][7] { return LsvsvZ; }

std::complex<double> (&Pythia8::CoupSUSY::RsvsvZ_ref__BOSS())[7][7] { return RsvsvZ; }

std::complex<double> (&Pythia8::CoupSUSY::LslsvW_ref__BOSS())[7][7] { return LslsvW; }

std::complex<double> (&Pythia8::CoupSUSY::RslsvW_ref__BOSS())[7][7] { return RslsvW; }

std::complex<double> (&Pythia8::CoupSUSY::LsvvX_ref__BOSS())[7][4][6] { return LsvvX; }

std::complex<double> (&Pythia8::CoupSUSY::RsvvX_ref__BOSS())[7][4][6] { return RsvvX; }

std::complex<double> (&Pythia8::CoupSUSY::LsllX_ref__BOSS())[7][4][6] { return LsllX; }

std::complex<double> (&Pythia8::CoupSUSY::RsllX_ref__BOSS())[7][4][6] { return RsllX; }

std::complex<double> (&Pythia8::CoupSUSY::LsvlX_ref__BOSS())[7][4][3] { return LsvlX; }

std::complex<double> (&Pythia8::CoupSUSY::RsvlX_ref__BOSS())[7][4][3] { return RsvlX; }

std::complex<double> (&Pythia8::CoupSUSY::LslvX_ref__BOSS())[7][4][3] { return LslvX; }

std::complex<double> (&Pythia8::CoupSUSY::RslvX_ref__BOSS())[7][4][3] { return RslvX; }

double (&Pythia8::CoupSUSY::rvLLE_ref__BOSS())[4][4][4] { return rvLLE; }

double (&Pythia8::CoupSUSY::rvLQD_ref__BOSS())[4][4][4] { return rvLQD; }

double (&Pythia8::CoupSUSY::rvUDD_ref__BOSS())[4][4][4] { return rvUDD; }

bool& Pythia8::CoupSUSY::isLLE_ref__BOSS() { return isLLE; }

bool& Pythia8::CoupSUSY::isLQD_ref__BOSS() { return isLQD; }

bool& Pythia8::CoupSUSY::isUDD_ref__BOSS() { return isUDD; }

std::complex<double> (&Pythia8::CoupSUSY::Rusq_ref__BOSS())[7][7] { return Rusq; }

std::complex<double> (&Pythia8::CoupSUSY::Rdsq_ref__BOSS())[7][7] { return Rdsq; }

std::complex<double> (&Pythia8::CoupSUSY::Rsl_ref__BOSS())[7][7] { return Rsl; }

std::complex<double> (&Pythia8::CoupSUSY::Rsv_ref__BOSS())[7][7] { return Rsv; }


#include "backend_types/Pythia_8_209/identification.hpp"

Pythia8::Abstract_CoupSUSY* Pythia8::CoupSUSY::pointerCopy__BOSS()
{
    Pythia8::Abstract_CoupSUSY* new_ptr = new Pythia8::CoupSUSY(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::CoupSUSY::pointerAssign__BOSS(Pythia8::Abstract_CoupSUSY* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::CoupSUSY* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<CoupSUSY*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
