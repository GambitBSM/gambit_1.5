//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Models for QCD axions and axion-like particles.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Sebastian Hoof
///          (s.hoof15@imperial.ac.uk)
///  \date 2019 Feb
///
///  *********************************************

#include <cmath>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/numerical_constants.hpp"

#include "gambit/Models/models/Axions.hpp"

#define MODEL GeneralCosmoALP
#define FRIEND DecayingDM_photon
void MODEL_NAMESPACE::GeneralCosmoALP_to_DecayingDM_photon (const ModelParameters &myparams, ModelParameters &friendparams)
{
    USE_MODEL_PIPE(FRIEND) // get pipe for "interpret as friend" function
    logger()<<"Running interpret_as_friend calculations for GeneralCosmoALP -> DecayingDM_photon ..."<<EOM;

    friendparams.setValue("lifetime", *Dep::lifetime);
    friendparams.setValue("mass", 1e-9*myparams["ma0"]); // Convert units from eV (GeneralCosmoALP) to GeV (DecayingDM_photon)
    friendparams.setValue("fraction", *Dep::DM_fraction);
}
#undef FRIEND

#define FRIEND etaBBN_rBBN_rCMB_dNurBBN_dNurCMB
void MODEL_NAMESPACE::GeneralCosmoALP_to_etaBBN_rBBN_rCMB_dNurBBN_dNurCMB (const ModelParameters & /*myparams*/, ModelParameters &friendparams)
{
    USE_MODEL_PIPE(FRIEND) // get pipe for "interpret as friend" function
    logger()<<"Running interpret_as_friend calculations for GeneralCosmoALP -> etaBBN_rBBN_rCMB_dNurBBN_dNurCMB ..."<<EOM;

    // Get the results for the entropy evolution between BBN and CMB
    // and map them onto r_BBN/CMB and dNurBBN/CMB.
    map_str_dbl Neff_results = *Dep::external_dNeff_etaBBN;

    // We assume that the initial ratio r (at BBN) is unchanged
    friendparams.setValue("r_BBN", 1.0);
    friendparams.setValue("r_CMB", pow(Neff_results["Neff_ratio"], 1./4.));
    // No dark raddiation
    friendparams.setValue("dNur_BBN", 0.0);
    friendparams.setValue("dNur_CMB", 0.0);
    friendparams.setValue("eta_BBN", (*Dep::eta0)*Neff_results["eta_ratio"]);
}
#undef FRIEND
#undef MODEL

#define MODEL CosmoALP
void MODEL_NAMESPACE::CosmoALP_to_GeneralCosmoALP (const ModelParameters &myparams, ModelParameters &parentparams)
{
    logger()<<"Running interpret_as_parent calculations for CosmoALP -> GeneralCosmoALP ..."<<EOM;

    const double alpha_red = alpha_EM/(2.0*pi);
    double fa  = myparams["fa"];

    parentparams.setValue("gagg", alpha_red*myparams["Cagg"]/fa);
    parentparams.setValue("gaee", 0);
    parentparams.setValue("fa", fa);
    parentparams.setValue("ma0", myparams["ma0"]);
    parentparams.setValue("Tchi", 1E99);
    parentparams.setValue("beta", 0);
    parentparams.setValue("thetai", myparams["thetai"]);
    parentparams.setValue("f0_thermal", myparams["f0_thermal"]);
    parentparams.setValue("T_R", myparams["T_R"]);
}
#undef MODEL

#define MODEL GeneralALP
void MODEL_NAMESPACE::GeneralALP_to_GeneralCosmoALP (const ModelParameters &myparams, ModelParameters &parentparams)
{
    logger()<<"Running interpret_as_parent calculations for GeneralALP -> GeneralCosmoALP ..."<<EOM;

    parentparams.setValue("gagg", myparams["gagg"]);
    parentparams.setValue("gaee", myparams["gaee"]);
    parentparams.setValue("fa", myparams["fa"]);
    parentparams.setValue("ma0", myparams["ma0"]);
    parentparams.setValue("Tchi", myparams["Tchi"]);
    parentparams.setValue("beta", myparams["beta"]);
    parentparams.setValue("thetai", myparams["thetai"]);
    // Set f0_thermal = 0, i.e. no thermal component.
    parentparams.setValue("f0_thermal", 0);
    // Use extremely high reheating temperature to guarantee that Tosc < T_R
    parentparams.setValue("T_R", 1.0E99);
}
#undef MODEL

#define MODEL QCDAxion
void MODEL_NAMESPACE::QCDAxion_to_GeneralALP (const ModelParameters &myparams, ModelParameters &parentparams)
{
    logger()<<"Running interpret_as_parent calculations for QCDAxion -> GeneralALP ..."<<EOM;

    const double alpha_red = alpha_EM/(2.0*pi);
    double fa  = myparams["fa"];
    double L2  = myparams["LambdaChi"]*myparams["LambdaChi"];
    double EoN = myparams["EoverN"];
    double CG  = myparams["CaggQCD"];

    parentparams.setValue("gagg", alpha_red*std::fabs(EoN-CG)/fa);
    parentparams.setValue("gaee", m_electron*myparams["Caee"]/fa);
    parentparams.setValue("fa", fa);
    parentparams.setValue("ma0", 1E+3*L2/fa);
    parentparams.setValue("Tchi", myparams["Tchi"]);
    parentparams.setValue("beta", myparams["beta"]);
    parentparams.setValue("thetai", myparams["thetai"]);
}
#undef MODEL

#define MODEL KSVZAxion
void MODEL_NAMESPACE::KSVZAxion_to_QCDAxion (const ModelParameters &myparams, ModelParameters &parentparams)
{
    logger()<<"Running interpret_as_parent calculations for KSVZAxion -> QCDAxion ..."<<EOM;

    const double prefactor = 3.0*alpha_EM*alpha_EM/(2.0*pi);
    const double scale     = 1.0;

    double EoN     = myparams["EoverN"];
    double CaggQCD = myparams["CaggQCD"];
    double fa      = myparams["fa"];

    parentparams.setValue("EoverN", EoN);
    parentparams.setValue("CaggQCD", CaggQCD);
    parentparams.setValue("Caee", prefactor*(EoN*std::log(fa/m_electron) - CaggQCD*std::log(scale/m_electron)));
    parentparams.setValue("fa", fa);
    parentparams.setValue("LambdaChi", myparams["LambdaChi"]);
    parentparams.setValue("Tchi", myparams["Tchi"]);
    parentparams.setValue("beta", myparams["beta"]);
    parentparams.setValue("thetai", myparams["thetai"]);
}
#undef MODEL

#define MODEL DFSZAxion_I
void MODEL_NAMESPACE::DFSZAxion_I_to_QCDAxion (const ModelParameters &myparams, ModelParameters &parentparams)
{
    logger()<<"Running interpret_as_parent calculations for DFSZAxion I -> QCDAxion ..."<<EOM;

    double angle = std::atan(myparams["tanbeta"]);
    double S2    = std::sin(angle);
           S2    = S2*S2;

    parentparams.setValue("EoverN", myparams["EoverN"]);
    parentparams.setValue("CaggQCD", myparams["CaggQCD"]);
    parentparams.setValue("Caee", S2/3.0);
    parentparams.setValue("fa", myparams["fa"]);
    parentparams.setValue("LambdaChi", myparams["LambdaChi"]);
    parentparams.setValue("Tchi", myparams["Tchi"]);
    parentparams.setValue("beta", myparams["beta"]);
    parentparams.setValue("thetai", myparams["thetai"]);
}
#undef MODEL

#define MODEL DFSZAxion_II
void MODEL_NAMESPACE::DFSZAxion_II_to_QCDAxion (const ModelParameters &myparams, ModelParameters &parentparams)
{
    logger()<<"Running interpret_as_parent calculations for DFSZAxion II -> QCDAxion ..."<<EOM;

    double angle = std::atan(myparams["tanbeta"]);
    double S2    = std::sin(angle);
           S2    = S2*S2;

    parentparams.setValue("EoverN", myparams["EoverN"]);
    parentparams.setValue("CaggQCD", myparams["CaggQCD"]);
    parentparams.setValue("Caee", (1.0-S2)/3.0);
    parentparams.setValue("fa", myparams["fa"]);
    parentparams.setValue("LambdaChi", myparams["LambdaChi"]);
    parentparams.setValue("Tchi", myparams["Tchi"]);
    parentparams.setValue("beta", myparams["beta"]);
    parentparams.setValue("thetai", myparams["thetai"]);
}
#undef MODEL

#define MODEL ConstantMassALP
void MODEL_NAMESPACE::ConstantMassALP_to_GeneralALP (const ModelParameters &myparams, ModelParameters &parentparams)
{
    logger()<<"Running interpret_as_parent calculations for ConstantMassALP -> GeneralALP ..."<<EOM;

    const double alpha_red = alpha_EM/(2.0*pi);

    double L2 = myparams["Lambda"]*myparams["Lambda"];
    double fa = myparams["fa"];

    parentparams.setValue("gagg", alpha_red*myparams["Cagg"]/fa);
    parentparams.setValue("gaee", m_electron*myparams["Caee"]/fa);
    parentparams.setValue("fa", fa);
    parentparams.setValue("ma0", 1E+3*L2/fa);
    parentparams.setValue("Tchi", 1.0E99);
    parentparams.setValue("beta", 0);
    parentparams.setValue("thetai", myparams["thetai"]);
}
#undef MODEL
