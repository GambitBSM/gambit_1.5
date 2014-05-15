#ifndef __ABSTRACT_X_HPP__
#define __ABSTRACT_X_HPP__
#include <iostream>

// MOVED TO: abstract_classes_extra.hpp
// class X;

// FORWARD DECLARED IN: abstract_classes_extra.hpp
class Abstract_X
{
    public:

        // Class methods
        virtual Abstract_T& t_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_T* getT_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_X* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual void pointerAssign_GAMBIT(Abstract_X *in){std::cout << "Called virtual function" << std::endl;}  

        Abstract_T* getT()
        {
            return getT_GAMBIT();
        }

        virtual void setT_GAMBIT(Abstract_T& t_in) {std::cout << "Called virtual function" << std::endl;};
        void setT(Abstract_T& t_in)
        {
            setT_GAMBIT(t_in);
        }
};

#endif  /* __ABSTRACT_X_HPP__ */
