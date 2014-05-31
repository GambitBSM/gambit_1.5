/*
 * A dummy C++ library for testing GAMBIT backend setup 
 * 
 * \author Anders Kvellestad
 * \date 2013-03-26
 *
 * Modified: Pat Scott 2013-06-23
 */

#include <iostream>
#include <string>
#include "libfirst_classloader.hpp"

//
// Some global variables:
//
int someInt;
double someDouble;
bool isInitialized = false;
std::string prefix = "libfirst: ";

class Base : public GAMBIT::Classes::LibFirst::v1_0::Base  // GAMBIT interface classes should always be inherited from publicly
{
  virtual double weighMe() {return 30.0;};
  public:
    double returnvalue() { return 3; }
    double valueMe(double pricePerKilo) {return 100.0;};
    Base();
    ~Base();
};

// Factory function for cattleClass objects (to be added by script)
Base* BaseFactory()
{ 
  return new Base();
}

// Some class that leaves its crap everywhere
class cattleClass : public Base, public GAMBIT::Classes::LibFirst::v1_0::cattleClass  // Script-like mod here: make cattleClass derive from GAMBIT::Classes::LibFirst::v1_0::cattleClass
{
  public: 
    cattleClass(double weight);
    ~cattleClass();
    double weighMe();
    double valueMe(double pricePerKilo);
  private:
    double weight;
};

// Factory function for cattleClass objects (to be added by script)
cattleClass* cattleClassFactory(double weight)
{ 
  return new cattleClass(weight);
}

// Constructor and destructor for base class
Base::Base() { std::cout << prefix << "Made a Base." << std::endl; };
Base::~Base() { std::cout << prefix << "Deleted a Base." << std::endl; };


//
// Some functions:
//

  
// 'initialization'
void initialize(int a)
{
  std::cout << std::endl;
  std::cout << prefix << "This is function 'initialize'." << std::endl;
  someInt = a;
  isInitialized = true;
  std::cout << prefix << "Initialization done. Variable 'someInt' set to: " << someInt << std::endl;
  std::cout << prefix << "Making some instances of the Base class " << std::endl;
  Base myBase;
  Base* Base2 = new Base;
  Base* Base3 = BaseFactory();
  std::cout << prefix << "Now I'm going to make some stack cattle..." << std::endl;
  cattleClass Daisy(500.0);
  cattleClass Buttercup(600.0);
  std::cout << prefix << "OK, that seemed to work." << std::endl;
  std::cout << prefix << "Now I'm going to make some heap cattle..." << std::endl;
  cattleClass* Bluebell = new cattleClass(500.0);
  cattleClass* Lacey = cattleClassFactory(400.0);
  std::cout << prefix << "OK, that seemed to work too." << std::endl;  
} // end initialize

// Constructor for the cattleClass
cattleClass::cattleClass(double weight_in) : weight(weight_in)
{
  std::cout << prefix << "Made cattleClass object." << std::endl;
}

// Destructor for the cattleClass
cattleClass::~cattleClass()
{
  std::cout << prefix << "Slaughter." << std::endl;
}


// 'calculation'
void someFunction()
{
  std::cout << std::endl;
  std::cout << prefix << "This is function 'someFunction'." << std::endl;

  if (isInitialized)
  {
    std::cout << prefix << "Will now perform a calculation..." << std::endl;
    someDouble = 3.1415*someInt;
    std::cout << prefix << "Result stored in variable 'someDouble' is: " << someDouble << std::endl;
  }
  else
  {
    std::cout << prefix << "Not initialized. Cannot perform calculation." << std::endl;
  }
} // end someFunction


// Weigh function for the cattleClass
double cattleClass::weighMe()
{
  std::cout << prefix << "Just got lotsa marrow." << std::endl;
  return weight;
}

// Pricing function for the cattleClass
double cattleClass::valueMe(double pricePerKilo)
{
  std::cout << prefix << "Pricecheck on heifer?" << std::endl;
  return weight*pricePerKilo;
}

// return 'result'
double returnResult()
{
  std::cout << "I'm returnResult() from libfirst.so, and I'm feeling well." << std::endl;
  return someDouble;
}


