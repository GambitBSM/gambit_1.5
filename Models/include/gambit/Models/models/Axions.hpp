//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  QCD-Axion and axion-like particles.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Sebastian Hoof
///  \date 2016 Oct
///  \date 2017 Feb, May, Jul
///  \date 2018 Feb
///
///  *********************************************

#include <string>
#include <cmath>

#include "gambit/Logs/logging.hpp"
#include "gambit/Utils/numerical_constants.hpp"

#ifndef __GeneralALP_hpp__
#define __GeneralALP_hpp__

// General axion model with parametric temperature-dependent mass
#define MODEL GeneralALP
  START_MODEL
  DEFINEPARS(gagg,gaee,fa,ma0,Tcrit,beta,thetai)
#undef MODEL

// QCD axion
// Parameter values are taken from results in 0910.1066
#define MODEL QCDAxion
#define PARENT GeneralALP
  START_MODEL
  DEFINEPARS(fa,Tcrit,beta,thetai)
  DEFINEPARS(LambdaQCD,EoverN,CaggQCD,Caee)

  INTERPRET_AS_PARENT_FUNCTION(QCDAxion_IAPfunc)

  void MODEL_NAMESPACE::QCDAxion_IAPfunc (const ModelParameters &myparams, ModelParameters &parentparams)
  {
      logger()<<"Running interpret_as_parent calculations for QCDAxion -> GeneralALP ..."<<EOM;

      const double alpha_red = alpha_EM/sqrt(2.0*pi);

      double fa  = myparams["fa"];
      double L2  = myparams["LambdaQCD"]*myparams["LambdaQCD"];
      double EoN = myparams["EoverN"];
      double CG  = myparams["CaggQCD"];

      // Values taken from 1511.02867, eq. 42
      //parentparams.setValue("gagg", 2.227E-12/FA);
      parentparams.setValue("gagg", 1E-9*alpha_red*std::fabs(EoN-CG)/fa);
      // Defined with Ce = cos(beta')^2/3
      parentparams.setValue("gaee", m_electron*myparams["Caee"]/fa);
      parentparams.setValue("fa", fa);
      parentparams.setValue("ma0", 1E+3*L2/fa);
      parentparams.setValue("Tcrit", myparams["Tcrit"]);
      parentparams.setValue("beta", myparams["beta"]);
      parentparams.setValue("thetai", myparams["thetai"]);
  }
#undef PARENT
#undef MODEL

// KSVZ axion model
#define MODEL KSVZAxion
#define PARENT QCDAxion
  START_MODEL
  DEFINEPARS(fa,Tcrit,beta,thetai,LambdaQCD,EoverN,CaggQCD)

  INTERPRET_AS_PARENT_FUNCTION(KSVZAxion_IAPfunc)

  void MODEL_NAMESPACE::KSVZAxion_IAPfunc (const ModelParameters &myparams, ModelParameters &parentparams)
  {
      logger()<<"Running interpret_as_parent calculations for KSVZAxion -> QCDAxion ..."<<EOM;

      const double prefactor = 3.0*alpha_EM*alpha_EM/(2.0*pi);
      const double scale     = 1.0;

      double EoN     = myparams["EoverN"];
      double LQCD    = myparams["LambdaQCD"];
      double CaggQCD = myparams["CaggQCD"];
      double fa      = myparams["fa"];

      parentparams.setValue("EoverN", EoN);
      parentparams.setValue("CaggQCD", CaggQCD);
      parentparams.setValue("Caee", prefactor*(EoN*std::log(fa/m_electron) - CaggQCD*std::log(scale/m_electron)));
      parentparams.setValue("fa", fa);
      parentparams.setValue("LambdaQCD", LQCD);
      parentparams.setValue("Tcrit", myparams["Tcrit"]);
      parentparams.setValue("beta", myparams["beta"]);
      parentparams.setValue("thetai", myparams["thetai"]);
  }
#undef PARENT
#undef MODEL

// DFSZ-I axion model
#define MODEL DFSZAxion_I
#define PARENT QCDAxion
  START_MODEL
  DEFINEPARS(fa,Tcrit,beta,thetai,LambdaQCD,EoverN,CaggQCD)
  DEFINEPARS(tanbeta)

  INTERPRET_AS_PARENT_FUNCTION(DFSZAxion_IAPfunc)

  void MODEL_NAMESPACE::DFSZAxion_IAPfunc (const ModelParameters &myparams, ModelParameters &parentparams)
  {
      logger()<<"Running interpret_as_parent calculations for DFSZAxion -> QCDAxion ..."<<EOM;

      double angle = std::atan(myparams["tanbeta"]);
      double S2    = std::sin(angle);
             S2    = S2*S2;

      parentparams.setValue("EoverN", myparams["EoverN"]);
      parentparams.setValue("CaggQCD", myparams["CaggQCD"]);
      parentparams.setValue("Caee", S2/3.0);
      parentparams.setValue("fa", myparams["fa"]);
      parentparams.setValue("LambdaQCD", myparams["LambdaQCD"]);
      parentparams.setValue("Tcrit", myparams["Tcrit"]);
      parentparams.setValue("beta", myparams["beta"]);
      parentparams.setValue("thetai", myparams["thetai"]);
  }
#undef PARENT
#undef MODEL

// DFSZ-II axion model
#define MODEL DFSZAxion_II
#define PARENT QCDAxion
  START_MODEL
  DEFINEPARS(fa,Tcrit,beta,thetai,LambdaQCD,EoverN,CaggQCD)
  DEFINEPARS(tanbeta)

  INTERPRET_AS_PARENT_FUNCTION(DFSZAxion_IAPfunc)

  void MODEL_NAMESPACE::DFSZAxion_IAPfunc (const ModelParameters &myparams, ModelParameters &parentparams)
  {
      logger()<<"Running interpret_as_parent calculations for DFSZAxion -> QCDAxion ..."<<EOM;

      double angle = std::atan(myparams["tanbeta"]);
      double S2    = std::sin(angle);
             S2    = S2*S2;

      parentparams.setValue("EoverN", myparams["EoverN"]);
      parentparams.setValue("CaggQCD", myparams["CaggQCD"]);
      parentparams.setValue("Caee", (1.0-S2)/3.0);
      parentparams.setValue("fa", myparams["fa"]);
      parentparams.setValue("LambdaQCD", myparams["LambdaQCD"]);
      parentparams.setValue("Tcrit", myparams["Tcrit"]);
      parentparams.setValue("beta", myparams["beta"]);
      parentparams.setValue("thetai", myparams["thetai"]);
  }
#undef PARENT
#undef MODEL

// Simple ALP model with not explicitely temperature-dependent mass
#define MODEL SimpleALP
#define PARENT GeneralALP
  START_MODEL
  DEFINEPARS(Cagg,Caee,fa,Lambda,thetai)

  INTERPRET_AS_PARENT_FUNCTION(SimpleALP_IAPfunc)

  void MODEL_NAMESPACE::SimpleALP_IAPfunc (const ModelParameters &myparams, ModelParameters &parentparams)
  {
      logger()<<"Running interpret_as_parent calculations for SimpleALP -> GeneralALP ..."<<EOM;

      const double alpha_red = alpha_EM/sqrt(2.0*pi);

      double L2 = myparams["Lambda"]*myparams["Lambda"];
      double FA = myparams["fa"];

      parentparams.setValue("gagg", 1E-9*alpha_red*myparams["Cagg"]/FA);
      parentparams.setValue("gaee", m_electron*myparams["Caee"]/FA);
      parentparams.setValue("fa", FA);
      parentparams.setValue("ma0", 1E+3*L2/FA);
      parentparams.setValue("Tcrit", 1.0E99);
      parentparams.setValue("beta", 0);
      parentparams.setValue("thetai", myparams["thetai"]);
  }
#undef PARENT
#undef MODEL

#endif /* defined(__GeneralAxion_hpp__) */
