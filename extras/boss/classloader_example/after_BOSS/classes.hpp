// BOSS modifies this file in place.
#ifndef __CLASSES_HPP__
#define __CLASSES_HPP__

#include <iostream>
#include "abstract_classes.hpp"
#include "abstract_classes_extra.hpp"

class T : public virtual Abstract_T
{
    public:
        // Member variables
        int i;
        double d;

        // Constructor
        T() : i(1), d(3.14) {}

        // Class methods
        void printMe();


        //
        // Generated members:
        //

        Abstract_T* pointerCopy_GAMBIT()
        {
            return new T(*this);
        }

        void pointerAssign_GAMBIT(Abstract_T *in)
        {
            *this = *dynamic_cast<T*>(in);
        }        

        int& i_ref_GAMBIT()
        {
            return i;
        }

        double& d_ref_GAMBIT()
        {
            return d;
        }

};



class X : public virtual Abstract_X
{
    public:
        // Member variables
        T t;

        // Constructor
        X() {}

        // Class methods
        T getT();

        void setT(T& t_in);

        //
        // Generated members:
        //

        Abstract_X* pointerCopy_GAMBIT()
        {
            return new X(*this);
        }

        void pointerAssign_GAMBIT(Abstract_X *in)
        {
            *this = *dynamic_cast<X*>(in);
        }             

        Abstract_T& t_ref_GAMBIT()
        {
            return t;
        }

        Abstract_T* getT_GAMBIT()
        {
            return new T(getT());
        }
    
        void setT_GAMBIT(Abstract_T& t_in)
        {
            setT( dynamic_cast<T&>(t_in) );
        }

};


template <typename Type>
class Container : public virtual Abstract_Container<Type>
{
    public:
        // Member variables
        Type var;

        // Constructor
        Container() {}

        Container(Type in) : var(in) {}

        void printMsg()
        {
            std::cout << std::endl;
            std::cout << "A message from class 'Container'." << std::endl;
            std::cout << std::endl;
        }

        //
        // Generated members:
        //

        Abstract_Container<Type>* pointerCopy_GAMBIT()
        {
            return new Container(*this);
        }

        void pointerAssign_GAMBIT(Abstract_Container<Type> *in)
        {
            *this = *dynamic_cast<Container*>(in);
        }             

        Type& var_ref_GAMBIT()
        {
            return var;
        }

};


#endif  /* __CLASSES_HPP__ */
