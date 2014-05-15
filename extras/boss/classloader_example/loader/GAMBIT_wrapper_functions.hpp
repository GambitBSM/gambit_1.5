#ifndef __GAMBIT_WRAPPER_FUNCTIONS_HPP__
#define __GAMBIT_WRAPPER_FUNCTIONS_HPP__

#include "GAMBIT_wrapper_classes.hpp"


// Function pointers to be filled by dynamic loading
void (*printT_WRAPPER)(Abstract_T&) = NULL;
void (*doubleT_WRAPPER)(Abstract_T&) = NULL;
Abstract_X* (*doubleX_byVal_WRAPPER)(Abstract_X&) = NULL;



void printT(T_gambit t_in)
{
	printT_WRAPPER(*t_in.BEptr);
}


void doubleT(T_gambit& t_in)
{
	doubleT_WRAPPER(*t_in.BEptr);
}


X_gambit doubleX_byVal(X_gambit x_in)
{
	return X_gambit( doubleX_byVal_WRAPPER(*x_in.BEptr) );
}


#endif /* __GAMBIT_WRAPPER_FUNCTIONS_HPP__ */ 