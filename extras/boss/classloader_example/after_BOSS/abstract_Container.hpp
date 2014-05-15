#ifndef __ABSTRACT_CONTAINER_HPP__
#define __ABSTRACT_CONTAINER_HPP__
#include <iostream>

// Presumably, this guy needs to know about the other Abstract_Classes.
// ... Oh, I get it. It should be in...
// MOVED TO: abstract_classes_extra.hpp
// class Abstract_T;
// class Abstract_X;

// FORWARD DECLARED IN: abstract_classes_extra.hpp
template <typename T1>
class Abstract_Container{};

template <>
class Abstract_Container<int>
{
    public:

        virtual int& var_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual void printMsg() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_Container<int>* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
 
        virtual void pointerAssign_GAMBIT(Abstract_Container<int> *in){std::cout << "Called virtual function" << std::endl;}  
};

// MOVED TO: abstract_classes_extra.hpp
/* template <> 
class Abstract_Container<Abstract_T>
{
    public:

        virtual Abstract_T& var_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual void printMsg() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_Container<Abstract_T>* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
 
        virtual void pointerAssign_GAMBIT(Abstract_Container<Abstract_T> *in){std::cout << "Called virtual function" << std::endl;}  
}; */

// MOVED TO: abstract_classes_extra.hpp
/* template <>
class Abstract_Container<Abstract_X>
{
    public:

        virtual Abstract_X& var_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual void printMsg() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_Container<Abstract_X>* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
 
        virtual void pointerAssign_GAMBIT(Abstract_Container<Abstract_X> *in){std::cout << "Called virtual function" << std::endl;}  
}; */

// MOVED TO: abstract_classes_extra.hpp
// template <> // This specialization is only needed to go from <X> to <Abstract_X>
// class Abstract_Container<X> : public Abstract_Container<Abstract_X>  
// {};

// MOVED TO: abstract_classes_extra.hpp
// template <> // This specialization is only needed to go from <T> to <Abstract_T>
// class Abstract_Container<T> : public Abstract_Container<Abstract_T>
// {};

#endif  /* __ABSTRACT_CONTAINER_HPP__ */
