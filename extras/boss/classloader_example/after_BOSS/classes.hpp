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

        Abstract_T* ptr_copy()
        {
            return new T(*this);
        }

        void ptr_assign(Abstract_T *in)
        {
            *this = *dynamic_cast<T*>(in);
        }        

        int& iRef();

        double& dRef();
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

        void setT(T t_in);

        void setT2(T** t_in, T** t_in2)
        {
            std::cout << " --- passing in pointers: " << *t_in << "  -  " << *t_in2 << std::endl;
            *t_in = *t_in2;
            std::cout << " --- returning pointers : " << *t_in << "  -  " << *t_in2 << std::endl;
            // std::cout << "This is setT2" << std::endl;
            // std::cout << "got pointer: " << t_in << std::endl;
            // t = **t_in;
        }

        //
        // Generated members:
        //

        Abstract_X* ptr_copy()
        {
            return new X(*this);
        }

        void ptr_assign(Abstract_X *in)
        {
            *this = *dynamic_cast<X*>(in);
        }             

        Abstract_T& tRef()
        {
            return t;
        }
        // T& tRef()
        // {
        //     return t;
        // }
        // Abstract_T& tRef_wrapper()
        // {
        //     return tRef();
        // }

        Abstract_T* getT_wrapper()
        {
            return new T(getT());
        }
    
        void setT_wrapper(Abstract_T& t_in)
        {
            setT( dynamic_cast<T&>(t_in) );
        }

        void setT2_wrapper(Abstract_T** t_in, Abstract_T** t_in2)
        {
            // setT2( dynamic_cast<T**>(t_in) );
            // setT2( dynamic_cast<T**>( &dynamic_cast<T*>(*t_in) );
            // setT2( t_in );
            std::cout << "This is setT2_wrapper" << std::endl;
            
            // T* t_in_tmp = dynamic_cast<T*>(*t_in);
            // T** t_in_tmp2 = &t_in_tmp;

            // T* t_in2_tmp = dynamic_cast<T*>(*t_in2);
            // T** t_in2_tmp2 = &t_in2_tmp;

            // setT2(t_in_tmp2, t_in2_tmp2);

            setT2( reinterpret_cast<T**>(t_in), reinterpret_cast<T**>(t_in2) );
        }

};


template <typename Type>
class Container : public virtual Abstract_Container<Type>
{
    public:

        Type var;

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

        Abstract_Container<Type>* ptr_copy()
        {
            return new Container(*this);
        }

        void ptr_assign(Abstract_Container<Type> *in)
        {
            *this = *dynamic_cast<Container*>(in);
        }             

        Type& varRef()
        {
            return var;
        }
};


#endif  /* __CLASSES_HPP__ */
