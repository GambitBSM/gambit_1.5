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

// Container<T> dummy_container_T;
// Container<X> dummy_container_X;
// Container<int> dummy_container_int;


// ---------- Test main function ----------------
// int main()
// {
//     std::cout << std::endl;

//     DummyNameSpace::U u;

//     // static const int arr[] = {1,2,3,4};
//     // std::vector<int> test_vec (arr, arr + sizeof(arr)/sizeof(arr[0]) );

//     // u.printVec(test_vec);
//     // std::cout << std::endl;

//     // u.changeVec(test_vec, 1, 10);
//     // u.printVec(test_vec);

//     static const MyInt arr[] = {MyInt(1),MyInt(2),MyInt(3),MyInt(4)};
//     std::vector<MyInt> test_vec (arr, arr + sizeof(arr)/sizeof(arr[0]) );

//     u.printVecMyInt(test_vec);
//     std::cout << std::endl;

//     return 0;
// }


