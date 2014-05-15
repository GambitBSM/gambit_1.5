#ifndef __CLASSES_HPP__
#define __CLASSES_HPP__

#include <iostream>

class T
{
    public:
        // Member variables
        int i;
        double d;

        // Constructor
        T() : i(1), d(3.14) {}

        // Class methods
        void printMe();
};



class X
{
    public:
        // Member variables
        T t;

        // Constructor
        X() {}

        // Class methods
        T getT();

        void setT(T& t_in);

        void setT2(T** t_in)
        {
            t = **t_in;
        }
};



template <typename Type>
class Container
{
    public:
        // Member variables
        Type var;

        // Constructor
        Container() {}

        void printMsg()
        {
            std::cout << std::endl;
            std::cout << "A message from class 'Container'." << std::endl;
            std::cout << std::endl;
        }
};


#endif  /* __CLASSES_HPP__ */
