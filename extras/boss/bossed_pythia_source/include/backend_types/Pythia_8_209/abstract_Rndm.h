#ifndef __abstract_Rndm_Pythia_8_209_h__
#define __abstract_Rndm_Pythia_8_209_h__

#include "gambit/Backends/abstractbase.hpp"
#include "forward_decls_abstract_classes.h"
#include "forward_decls_wrapper_classes.h"
#include <vector>
#include <string>
#include <cstddef>

#ifndef ENUMS_DECLARED
#define ENUMS_DECLARED
#include "enum_decl_copies.h"
#endif

#include "identification.hpp"

// Forward declaration needed by the destructor pattern.
void wrapper_deleter(CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Rndm*);


namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    
    namespace Pythia8
    {
        class Abstract_Rndm : virtual public AbstractBase
        {
            private:
                // IGNORED: Variable  -- Name: DEFAULTSEED  -- XML id: _26383
                // IGNORED: Field  -- Name: initRndm  -- XML id: _26384
                // IGNORED: Field  -- Name: i97  -- XML id: _26385
                // IGNORED: Field  -- Name: j97  -- XML id: _26386
                // IGNORED: Field  -- Name: seedSave  -- XML id: _26387
                // IGNORED: Field  -- Name: sequence  -- XML id: _26388
                // IGNORED: Field  -- Name: u  -- XML id: _26389
                // IGNORED: Field  -- Name: c  -- XML id: _26390
                // IGNORED: Field  -- Name: cd  -- XML id: _26391
                // IGNORED: Field  -- Name: cm  -- XML id: _26392
                // IGNORED: Field  -- Name: useExternalRndm  -- XML id: _26393
                // IGNORED: Field  -- Name: rndmEngPtr  -- XML id: _26394
            public:
    
                virtual void init(int) =0;
    
                virtual void init__BOSS() =0;
    
                virtual double flat() =0;
    
                virtual double exp() =0;
    
                virtual double xexp() =0;
    
                virtual double gauss() =0;
    
                virtual std::pair<double,double> gauss2() =0;
    
                virtual int pick(const std::vector<double,std::allocator<double> >&) =0;
    
                virtual bool dumpState(std::basic_string<char,std::char_traits<char>,std::allocator<char> >) =0;
    
                virtual bool readState(std::basic_string<char,std::char_traits<char>,std::allocator<char> >) =0;
    
            public:
                virtual void pointerAssign__BOSS(Abstract_Rndm*) =0;
                virtual Abstract_Rndm* pointerCopy__BOSS() =0;
    
            private:
                mutable Rndm* wptr;
    
            public:
                Abstract_Rndm()
                {
                }
    
                void wrapper__BOSS(Rndm* wptr_in)
                {
                    wptr = wptr_in;
                    is_wrapped(true);
                    can_delete_wrapper(true);
                }
    
                Rndm* wrapper__BOSS()
                {
                    return wptr;
                }
    
                virtual ~Abstract_Rndm()
                {
                    if (can_delete_wrapper())
                    {
                        can_delete_me(false);
                        wrapper_deleter(wptr);
                    }
                }
        };
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"


#endif /* __abstract_Rndm_Pythia_8_209_h__ */
