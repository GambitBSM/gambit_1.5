
/// Main class 
class outer 
{

  private:

    /// Abstract version of first inner class 
    class inner1
    {
      public:
        virtual double get_f_of_mA(double) = 0;
        virtual void set_mA(double) = 0;
    };

    /// Abstract version of second inner class 
    class inner2
    {
      public:
        virtual double get_f_of_mB(double, double) = 0;
        virtual void set_mB(double) = 0;
    };


  protected:

    /// Constructor
    outer(inner1& inner1instance, inner2& inner2instance) : Insides1 (inner1instance), Insides2 (inner2instance) {}


  public:

    /// Internal references to instances of the inner classes.
    inner1& Insides1;
    inner2& Insides2;

};


/// Specialised, derived version of outer class 
template <typename T>
class specialised_outer : public outer
{

  private:

    /// Specialised, derived version of first inner class 
    template <typename U>
    class specialised_inner1 : public outer::inner1
    { 
      private:
        U mA;
      public:
        void set_mA(double input_mA)     { mA = U(input_mA); }
        double get_f_of_mA(double par1) { return double(mA) * par1; }
    };   

    /// Specialised, derived version of second inner class 
    template <typename U>
    class specialised_inner2 : public outer::inner2
    { 
      private:
        U mB;
      public:
        void set_mB(double input_mB)                 { mB = U(input_mB); }
        double get_f_of_mB(double par1, double par2) { return double(mB) * par1 + par2; }
    };   

    /// Internal instances of the derived inner classes
    specialised_inner1<T> myInner1;
    specialised_inner2<T> myInner2;

  public:

    /// Only public method or variable is the constructor.
    specialised_outer() : outer(myInner1, myInner2) {}

};

