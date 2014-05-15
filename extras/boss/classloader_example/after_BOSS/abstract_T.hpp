#ifndef __ABSTRACT_T_HPP__
#define __ABSTRACT_T_HPP__
#include <iostream>

// MOVED TO: abstract_classes_extra.hpp
// class T;

// FORWARD DECLARED IN: abstract_classes_extra.hpp
class Abstract_T
{
    public:

        // Class methods
        virtual void printMe() {std::cout << "Called virtual function" << std::endl;}

        virtual int& i_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_T* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual void pointerAssign_GAMBIT(Abstract_T *in){std::cout << "Called virtual function" << std::endl;}     

        virtual double& d_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
};

#endif  /* __ABSTRACT_T_HPP__ */
