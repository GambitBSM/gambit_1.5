//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  GAMBIT models for nuclear level density
///  functions
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Anders Kvellestad
///          anders.kvellestad@fys.uio.no
///  \date 2020 Mar
///
///  *********************************************


#ifndef __nldmodels_hpp__
#define __nldmodels_hpp__

#define MODEL NLDModelCT_and_discretes
  START_MODEL
  DEFINEPARS(nld_T, nld_Eshift, nld_Ecrit)
#undef MODEL

#define MODEL NLDModel10
  START_MODEL
  DEFINEPARS(nld_p1,nld_p2,nld_p3,nld_p4,nld_p5,nld_p6,nld_p7,nld_p8,nld_p9,nld_p10)
#undef MODEL

#define MODEL NLDModel15
  START_MODEL
  DEFINEPARS(nld_p1,nld_p2,nld_p3,nld_p4,nld_p5,nld_p6,nld_p7,nld_p8,nld_p9,nld_p10,nld_p11,nld_p12,nld_p13,nld_p14,nld_p15)
#undef MODEL

#define MODEL NLDModel20
  START_MODEL
  DEFINEPARS(nld_p1,nld_p2,nld_p3,nld_p4,nld_p5,nld_p6,nld_p7,nld_p8,nld_p9,nld_p10,nld_p11,nld_p12,nld_p13,nld_p14,nld_p15,nld_p16,nld_p17,nld_p18,nld_p19,nld_p20)
#undef MODEL

#endif /* defined(__nldmodels_hpp__) */














