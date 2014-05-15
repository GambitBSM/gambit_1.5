#include <iostream>
#include "classes.hpp"


//
// Class T member methods 
//

void T::printMe()
{
    std::cout << "This is printMe from class T: " << std::endl;
    std::cout << "  i = " << i << std::endl;
    std::cout << "  d = " << d << std::endl;
}


//
// Class X member methods 
//

T X::getT()
{
    return t;
}


void X::setT(T t_in)
{
    t = t_in;
}


//
// Create some dummy instances to make sure these classes are compiled
//

Container<T> dummy_container_T;
Container<X> dummy_container_X;
Container<int> dummy_container_int;
