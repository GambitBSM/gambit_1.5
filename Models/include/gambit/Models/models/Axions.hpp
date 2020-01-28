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
///  \date 2016 Oct
///  \date 2017 Feb, May, Jul
///  \date 2018 Feb
///  \date 2019 Feb
///
///  *********************************************

#ifndef __GeneralALP_hpp__
#define __GeneralALP_hpp__

#include "gambit/Models/models/CosmoEnergyInjection.hpp"
#include "gambit/Models/models/CosmoModels.hpp"

// General axion model with parametric temperature-dependent mass and cosmological applications.
#define MODEL GeneralCosmoALP
  START_MODEL
  // Physical units: gagg [GeV^-1], gaee [dimensionless], fa [GeV], ma0 [eV]
  //                 Tchi [MeV], beta [dimensionless], thetai [dimensionless], f0_thermal [dimensionless], T_R [MeV]
  DEFINEPARS(gagg,gaee,fa,ma0,Tchi,beta,thetai)
  DEFINEPARS(f0_thermal, T_R)
  
  // Friendship with "DecayingDM_photon" (Mapping is defined in Axions.cpp)
  // (Energy injection into CMB)
  INTERPRET_AS_X_FUNCTION(DecayingDM_photon,GeneralCosmoALP_to_DecayingDM_photon)
  // The mapping CosmoALP_to_DecayingDM_photon depends on the lifetime and the fraction rho_a/rho_cdm (mapping of the mass is trivial).
  INTERPRET_AS_X_DEPENDENCY(DecayingDM_photon,lifetime,double)
  INTERPRET_AS_X_DEPENDENCY(DecayingDM_photon,DM_fraction,double)
  
  // Friendship with "etaBBN_rBBN_rCMB_dNurBBN_dNurCMB" (Mapping is defined in Axions.cpp)
  INTERPRET_AS_X_FUNCTION(etaBBN_rBBN_rCMB_dNurBBN_dNurCMB,GeneralCosmoALP_to_etaBBN_rBBN_rCMB_dNurBBN_dNurCMB)
  INTERPRET_AS_X_DEPENDENCY(etaBBN_rBBN_rCMB_dNurBBN_dNurCMB,external_dNeff_etaBBN,map_str_dbl)
  INTERPRET_AS_X_DEPENDENCY(etaBBN_rBBN_rCMB_dNurBBN_dNurCMB,eta0,double)
#undef MODEL

// Simplified general axion model with parametric temperature-independent mass and cosmological applications.
#define MODEL CosmoALP
#define PARENT GeneralCosmoALP
  START_MODEL
  // Units for these parameters are the same as for the GeneralCosmoALP.
  DEFINEPARS(Cagg,fa,ma0,thetai)
  DEFINEPARS(f0_thermal,T_R)
  // Translation to parent, all defined in Axions.cpp:
  INTERPRET_AS_PARENT_FUNCTION(CosmoALP_to_GeneralCosmoALP)
#undef PARENT
#undef MODEL

// General axion model with parametric temperature-dependent mass.
#define MODEL GeneralALP
#define PARENT GeneralCosmoALP
  START_MODEL
  // Physical units: gagg [GeV^-1], gaee [dimensionless], fa [GeV], ma0 [eV]
  //                 Tchi [MeV], beta [dimensionless], thetai [dimensionless]
  DEFINEPARS(gagg,gaee,fa,ma0,Tchi,beta,thetai)
  // Translation to parent, all defined in Axions.cpp:
  INTERPRET_AS_PARENT_FUNCTION(GeneralALP_to_GeneralCosmoALP)
#undef PARENT
#undef MODEL

// QCD axion model
#define MODEL QCDAxion
#define PARENT GeneralALP
  START_MODEL
  // Units for these parameters are the same as for the GeneralALP.
  DEFINEPARS(fa,Tchi,beta,thetai)
  // Physical units: LambdaChi [MeV], EoverN [dimensionless], CaggQCD [dimensionless]
  //                 Caee [dimensionless]
  DEFINEPARS(LambdaChi,EoverN,CaggQCD,Caee)
  // Translation to parent, all defined in Axions.cpp:
  INTERPRET_AS_PARENT_FUNCTION(QCDAxion_to_GeneralALP)
#undef PARENT
#undef MODEL

// KSVZ axion model
#define MODEL KSVZAxion
#define PARENT QCDAxion
  START_MODEL
  // Units for these parameters are the same as for the QCDAxion.
  DEFINEPARS(fa,Tchi,beta,thetai,LambdaChi,EoverN,CaggQCD)
  INTERPRET_AS_PARENT_FUNCTION(KSVZAxion_to_QCDAxion)
#undef PARENT
#undef MODEL

// DFSZ-I axion model
#define MODEL DFSZAxion_I
#define PARENT QCDAxion
  START_MODEL
  // Units for these parameters are the same as for the QCDAxion.
  DEFINEPARS(fa,Tchi,beta,thetai,LambdaChi,EoverN,CaggQCD)
  // Physical units: tanbeta [dimensionless]
  DEFINEPARS(tanbeta)
  INTERPRET_AS_PARENT_FUNCTION(DFSZAxion_I_to_QCDAxion)
#undef PARENT
#undef MODEL

// DFSZ-II axion model
#define MODEL DFSZAxion_II
#define PARENT QCDAxion
  START_MODEL
  // Units for these parameters are the same as for the QCDAxion.
  DEFINEPARS(fa,Tchi,beta,thetai,LambdaChi,EoverN,CaggQCD)
  // Physical units: tanbeta [dimensionless]
  DEFINEPARS(tanbeta)
  INTERPRET_AS_PARENT_FUNCTION(DFSZAxion_II_to_QCDAxion)
#undef PARENT
#undef MODEL

// ConstantMassALP model with temperature-independent mass and QCD-axion-inspired couplings
#define MODEL ConstantMassALP
#define PARENT GeneralALP
  START_MODEL
  // Units for these parameters are the same as for the GeneralALP.
  DEFINEPARS(fa,thetai)
  // Physical units: Cagg [dimensionless], Caee [dimensionless], Lambda [MeV]
  DEFINEPARS(Cagg,Caee,Lambda)
  INTERPRET_AS_PARENT_FUNCTION(ConstantMassALP_to_GeneralALP)
#undef PARENT
#undef MODEL

#endif
