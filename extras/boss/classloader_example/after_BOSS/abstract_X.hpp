#ifndef __ABSTRACT_X_HPP__
#define __ABSTRACT_X_HPP__
#include <iostream>

#include "abstract_T.hpp"


class Abstract_X
{
    public:

        // Class methods
        virtual Abstract_T& t_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_X* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual void pointerAssign_GAMBIT(Abstract_X *in) {std::cout << "Called virtual function" << std::endl;}  

        virtual Abstract_T* getT_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
        Abstract_T* getT()
        {
            return getT_GAMBIT();
        }

        virtual void setT_GAMBIT(Abstract_T& t_in) {std::cout << "Called virtual function" << std::endl;};
        void setT(Abstract_T& t_in)
        {
            setT_GAMBIT(t_in);
        }

        virtual void refTest_GAMBIT(Abstract_T& t_in, int& i_in) {std::cout << "Called virtual function" << std::endl;};
        void refTest(Abstract_T& t_in, int& i_in)
        {
            refTest_GAMBIT(t_in, i_in);
        }

        virtual int**& testFunc_GAMBIT(Abstract_T* t1, Abstract_T& t2, int**& ipp, double d) {std::cout << "Called virtual function" << std::endl;};
        int**& testFunc(Abstract_T* t1, Abstract_T& t2, int**& ipp, double d)
        {
            return testFunc_GAMBIT(t1, t2, ipp, d);
        }

        virtual ~Abstract_X()
        {
            std::cout << "(Destructor of Abstract_X)" << std::endl;
        }

};

#endif  /* __ABSTRACT_X_HPP__ */
