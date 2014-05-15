#include <iostream>
#include "classes.hpp"

//
// Global functions
//
void printT(T t_in)
{
	std::cout << std::endl;
	std::cout << "This is global function 'printT(T)':" << std::endl;
	std::cout << " member i = " << t_in.i << std::endl;
	std::cout << " member d = " << t_in.d << std::endl;
	std::cout << std::endl;
}


void doubleT(T& t_in)
{
	std::cout << std::endl;
	std::cout << "This is global function 'doubleT(T&)':" << std::endl;
	std::cout << std::endl;
	
	t_in.i = 2*t_in.i;
	t_in.d = 2*t_in.d;
}


X doubleX_byVal(X x_in)
{
	std::cout << std::endl;
	std::cout << "This is global function 'doubleX_byVal(X)':" << std::endl;
	std::cout << std::endl;

	X x_new = x_in;
	x_new.t.i = 2*x_new.t.i;
	x_new.t.d = 2*x_new.t.d;
	return x_new;
}


//
// Generated wrapper functions
//
void printT_GAMBIT(Abstract_T& abs_t_in)
{
	printT( dynamic_cast<T&>(abs_t_in));
}


void doubleT_GAMBIT(Abstract_T& abs_t_in)
{
	doubleT( dynamic_cast<T&>(abs_t_in));
}


Abstract_X* doubleX_byVal_GAMBIT(Abstract_X& abs_x_in)
{
	X temp_x = doubleX_byVal( dynamic_cast<X&>(abs_x_in) );
	Abstract_X* abs_x_new = new X(temp_x);
	return abs_x_new;
}
