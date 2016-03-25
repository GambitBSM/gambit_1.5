#ifndef __abstract_AlphaStrong_Pythia_8_212_h__
#define __abstract_AlphaStrong_Pythia_8_212_h__

#include "gambit/Backends/abstractbase.hpp"
#include "forward_decls_abstract_classes.h"
#include "forward_decls_wrapper_classes.h"
#include <cstddef>

#include "identification.hpp"

// Forward declaration needed by the destructor pattern.
void set_delete_BEptr(CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::AlphaStrong*, bool);


// Forward declaration needed by the destructor pattern.
void wrapper_deleter(CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::AlphaStrong*);


// Forward declaration for wrapper_creator.
void wrapper_creator(CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Abstract_AlphaStrong*);


namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    
    namespace Pythia8
    {
        class Abstract_AlphaStrong : public virtual AbstractBase
        {
            public:
    
                virtual void init(double, int, int, bool) =0;
    
                virtual void init__BOSS(double, int, int) =0;
    
                virtual void init__BOSS(double, int) =0;
    
                virtual void init__BOSS(double) =0;
    
                virtual void init__BOSS() =0;
    
                virtual void setThresholds(double, double, double) =0;
    
                virtual double alphaS(double) =0;
    
                virtual double alphaS1Ord(double) =0;
    
                virtual double alphaS2OrdCorr(double) =0;
    
                virtual double Lambda3() const =0;
    
                virtual double Lambda4() const =0;
    
                virtual double Lambda5() const =0;
    
                virtual double Lambda6() const =0;
    
                virtual double muThres(int) =0;
    
                virtual double muThres2(int) =0;
    
                virtual double facCMW(int) =0;
    
            public:
                virtual void pointer_assign__BOSS(Abstract_AlphaStrong*) =0;
                virtual Abstract_AlphaStrong* pointer_copy__BOSS() =0;
    
            private:
                AlphaStrong* wptr;
                bool delete_wrapper;
            public:
                AlphaStrong* get_wptr() { return wptr; }
                void set_wptr(AlphaStrong* wptr_in) { wptr = wptr_in; }
                bool get_delete_wrapper() { return delete_wrapper; }
                void set_delete_wrapper(bool del_wrp_in) { delete_wrapper = del_wrp_in; }
    
            public:
                Abstract_AlphaStrong()
                {
                    wptr = 0;
                    delete_wrapper = false;
                }
    
                Abstract_AlphaStrong(const Abstract_AlphaStrong&)
                {
                    wptr = 0;
                    delete_wrapper = false;
                }
    
                Abstract_AlphaStrong& operator=(const Abstract_AlphaStrong&) { return *this; }
    
                virtual void init_wrapper()
                {
                    if (wptr == 0)
                    {
                        wrapper_creator(this);
                        delete_wrapper = true;
                    }
                }
    
                AlphaStrong* get_init_wptr()
                {
                    init_wrapper();
                    return wptr;
                }
    
                AlphaStrong& get_init_wref()
                {
                    init_wrapper();
                    return *wptr;
                }
    
                virtual ~Abstract_AlphaStrong()
                {
                    if (wptr != 0)
                    {
                        set_delete_BEptr(wptr, false);
                        if (delete_wrapper == true)
                        {
                            wrapper_deleter(wptr);
                            wptr = 0;
                            delete_wrapper = false;
                        }
                    }
                }
        };
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"


#endif /* __abstract_AlphaStrong_Pythia_8_212_h__ */
