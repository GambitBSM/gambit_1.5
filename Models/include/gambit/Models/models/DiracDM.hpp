//  GAMBIT: Global and Modular BSM Inference Tool
//  *********************************************
//
//  Dirac DM
//
//  *********************************************
//
//  Authors
//  =======
//
//  (add name and date if you modify)
//
//  Ankit Beniwal
//  2016 August, 2017 June
//
//  Sebastian Wild
//  2018 January
//
//  Sanjay Bloor
//  2018 August
//
//  *********************************************

#ifndef __DiracSingletDM_Z2_hpp__
#define __DiracSingletDM_Z2_hpp__

#define MODEL DiracSingletDM_Z2
  START_MODEL
  DEFINEPARS(mF, lF, xi)
#undef MODEL

#define MODEL DiracSingletDM_Z2_sps
#define PARENT DiracSingletDM_Z2
  START_MODEL
  DEFINEPARS(mF, lF_s, lF_ps)
  INTERPRET_AS_PARENT_FUNCTION(DiracSingletDM_Z2_sps_to_DiracSingletDM_Z2)
#undef PARENT
#undef MODEL

#endif
