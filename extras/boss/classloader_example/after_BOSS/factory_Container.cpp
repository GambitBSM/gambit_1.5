#include "classes.hpp"

Abstract_Container<Abstract_X>* Factory_Container_X()
{
	return new Container<X>();
}

Abstract_Container<int>* Factory_Container_int()
{
	return new Container<int>();
}

Abstract_Container<Abstract_T>* Factory_Container_T()
{
	return new Container<T>();
}
