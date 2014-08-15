#ifndef __WRAPPERBASE__
#define __WRAPPERBASE__

#include <memory>

template <typename T>
class Deleter
{
    protected:
        bool member_variable;

    public:
        Deleter() : member_variable(false) {}
        Deleter(bool memvar_in) : member_variable(memvar_in) {}

        void setMemberVariable(bool memvar_in)
        {
            member_variable = memvar_in;
        }

        void operator()(T* ptr)
        {
            if (member_variable==false) { delete ptr; }
        }
};

template <typename T>
class WrapperBase
{
    public:
        std::shared_ptr<T> BEptr;

        // Constructor
        WrapperBase(T* BEptr_in, bool memvar_in)
        {
            BEptr.reset(BEptr_in, Deleter<T>(memvar_in));
        }

        // Special member function to set member_variable in Deleter: 
        void _setMemberVariable(bool memvar_in)
        {
            std::get_deleter<Deleter<T> >(BEptr)->setMemberVariable(memvar_in);
        }
};

#endif /* __WRAPPERBASE__ */
