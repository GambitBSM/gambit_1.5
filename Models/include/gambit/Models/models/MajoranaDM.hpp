//  GAMBIT: Global and Modular BSM Inference Tool
//  *********************************************
//
//  Majorana DM 
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
//  *********************************************

#ifndef __MajoranaDM_hpp__
#define __MajoranaDM_hpp__

#define MODEL MajoranaDM
  START_MODEL
  DEFINEPARS(mX, lX, cosXI)
#undef MODEL

#define MODEL MajoranaDM_sps
#define PARENT MajoranaDM
  START_MODEL
  DEFINEPARS(mX, lX_s, lX_ps)
  INTERPRET_AS_PARENT_FUNCTION(MajoranaDM_sps_to_MajoranaDM)
#undef PARENT
#undef MODEL

#define MODEL MajoranaDM_xi
#define PARENT MajoranaDM
  START_MODEL
  DEFINEPARS(mX, lX, xi)
  INTERPRET_AS_PARENT_FUNCTION(MajoranaDM_xi_to_MajoranaDM)
#undef PARENT
#undef MODEL
  
#endif
