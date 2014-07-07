#ifndef __CLASSES_HPP__
#define __CLASSES_HPP__

#include <iostream>
#include <string>
#include <vector>

// namespace MyNamespace
// {
//     template<typename T>
//     class SomeClass {};

//     SomeClass<int> dummy_instance;
// }

class MyInt
{
    public:
        
        int val;
        
        MyInt(int val_in) : val(val_in)
        {}
};


namespace DummyNameSpace
{
    class U
    {
        public:

            U()
            {
                std::cout << "(Constructor of U)" << std::endl;
            }

            ~U()
            {
                std::cout << "(Destructor of U)" << std::endl;
            }

            // void memberFunc(std::string s_in, std::vector<int> a)
            // {
            //     std::cout << "This is memberFunc from class U" << std::endl;
            // }

            void printStr(std::string s_in)
            {
                std::cout << "This is printStr in class U. Got string: " << std::endl; 
                std::cout << s_in << std::endl;
            }

            void printVec(std::vector<int> v_in)
            {
                std::cout << "This is printVec in class U. Got vector: " << std::endl;
                for (std::vector<int>::const_iterator i = v_in.begin(); i != v_in.end(); ++i)
                {
                    std::cout << " " << *i;
                }
                std::cout << std::endl;
            }

            // void printVecMyInt(std::vector<MyInt> v_in)
            // {
            //     std::cout << "This is printVecMyInt in class U. Got vector: " << std::endl;
            //     for (std::vector<MyInt>::const_iterator i = v_in.begin(); i != v_in.end(); ++i)
            //     {
            //         std::cout << " " << (*i).val;
            //     }
            //     std::cout << std::endl;
            // }

            std::vector<int>& changeVec(std::vector<int>& v_in, int el, int new_val)
            {
                std::cout << "This is changeVec in class U. Changing vector content... " << std::endl;

                v_in[el] = new_val;

                std::cout << "...done. Now returning the modified vector." << std::endl;

                return v_in;
            }

    };
}


class T : public DummyNameSpace::U, private std::string
{
    public:
        // Member variables
        int i;
        double d;

        // Constructor
        T() : i(1), d(3.14)
        {
            std::cout << "(Constructor of T)" << std::endl;
        }
        
        T(int i_in, double d_in) : i(i_in), d(d_in) 
        {
            std::cout << "(Constructor of T)" << std::endl;
        }

        // Destructor
        ~T()
        {
            std::cout << "(Destructor of T)" << std::endl;
        }

        // Class methods
        void printMe();
};



class X
{
    public:
        // Member variables
        T t;

        // Constructor
        X()
        {
            std::cout << "(Constructor of X)" << std::endl;
        }

        X(T t_in) : t(t_in)
        {
            std::cout << "(Constructor of X)" << std::endl;
        }

        // Destructor
        ~X()
        {
            std::cout << "(Destructor of X)" << std::endl;
        }

        // Class methods
        T getT();

        void setT(T t_in);

        void refTest(T& t_in, int& i_in)
        {
            T new_t;
            new_t.i = 123;
            new_t.d = 1.23;
            t_in = new_t;

            int new_i = 987;
            i_in = new_i;
        }

        int**& testFunc(T* t1, T t2, int**& ipp, double d)
        {
            **ipp += 1;
            return ipp;
        }
};


// template <typename Type>
// class Container
// {
//     public:
//         // Member variables
//         Type var;

//         // Constructor
//         Container() {}

//         Container(Type in) : var(in) {}

//         void printMsg()
//         {
//             std::cout << std::endl;
//             std::cout << "A message from class 'Container'." << std::endl;
//             std::cout << std::endl;
//         }

// };

#endif  /* __CLASSES_HPP__ */
