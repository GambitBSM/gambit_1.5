/*
 * Master class(es) for class loading with libfirst
 * To be generated automatically by script, and
 * included in both GAMBIT and the backend code.
 * 
 * \author Pat Scott
 * \date 2013-06-23
 * 
 */

namespace GAMBIT
{ namespace Classes
  { namespace LibFirst
    { namespace v1_0
      { class Base
        { virtual double weighMe() {return 0;} ;
          public:
            virtual double returnvalue() {return 0;};
            virtual ~Base() {};
            virtual double valueMe(double) {return 0;};
        };
        class cattleClass : public Base  // the inheritance mode here needs to match the mode with which the true class inherits from the true parent class
        { public:   // the access modes and parameter signature for all functions need to match those of the true class
            virtual ~cattleClass() {};
            virtual double weighMe() {return 0;};  // and for some reason all the methods need to have implementations, not just declarations
            virtual double valueMe(double) {return 0;};
        };
      }
    }
  }
}

