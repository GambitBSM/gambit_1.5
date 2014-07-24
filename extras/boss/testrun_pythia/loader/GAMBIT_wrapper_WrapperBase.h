#ifndef __GAMBIT_WRAPPER_WRAPPERBASE_H__
#define __GAMBIT_WRAPPER_WRAPPERBASE_H__

template <typename T>
class WrapperBase
{
    protected:
        bool member_variable;
    public:
        T* BEptr;

        // Constructor
        WrapperBase(T* BEptr_in, bool memvar_in) : BEptr(BEptr_in), member_variable(memvar_in) {}

        // Special member function to set member_variable: 
        void _set_member_variable(bool memvar_in) { member_variable = memvar_in; }

        // Destructor: 
        ~WrapperBase()
        {
            if(member_variable==false) { delete BEptr; }
        }
};

#endif /* __GAMBIT_WRAPPER_WRAPPERBASE_H__ */
