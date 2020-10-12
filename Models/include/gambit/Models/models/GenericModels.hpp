//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  GAMBIT generic models.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///   
///  \author Anders Kvellestad
///          a.kvellestad@imperial.ac.uk
///  \date 2018 Oct
///
///  *********************************************


#ifndef __genericmodels_hpp__
#define __genericmodels_hpp__

// This is a generic model with 5 free parameters
#define MODEL GenericModel5
  START_MODEL
  DEFINEPARS(p1,p2,p3,p4,p5)
#undef MODEL

// This is a generic model with 10 free parameters
#define MODEL GenericModel10
  START_MODEL
  DEFINEPARS(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
#undef MODEL

// This is a generic model with 15 free parameters
#define MODEL GenericModel15
  START_MODEL
  DEFINEPARS(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15)
#undef MODEL

// This is a generic model with 20 free parameters
#define MODEL GenericModel20
  START_MODEL
  DEFINEPARS(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20)
#undef MODEL

#endif /* defined(__genericmodels_hpp__) */














