#ifndef __ABSTRACT_U_HPP__
#define __ABSTRACT_U_HPP__
#include <iostream>

class Abstract_U
{
    public:

        virtual void memberFunc() {std::cout << "Called virtual function" << std::endl;};

    public:
        virtual void pointerAssign_GAMBIT(Abstract_U* in) {std::cout << "Called virtual function" << std::endl;};
        virtual Abstract_U* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
        
        virtual ~Abstract_U()
        {
        	std::cout << "(Destructor of Abstract_U)" << std::endl;
        }

};


#endif  /* __ABSTRACT_U_HPP__ */
