#ifndef __CalculatedState__
#define __CalculatedState__

/* These CalculatedState classes are used by the modules, backends, etc. 
//     They provide an easy way to check if the module/backend/whatever
//     has already calculated the value of some observable. ---Abram  
// NOTE: Once these classes are fully ready, neither the users nor anyone in
//     the collaboration needs to change them. Thus, to see how useful these
//     classes are, please read only the demo files. */

#include <iostream>
#include <boost/type_traits/is_same.hpp>
#include "Exceptions.hpp"

using namespace std;

////// begin: CalculatedState1
//    the calculated state for a single ObjectOfInterest with a tag 
template<typename Tag1>
struct CalculatedState1 {
    bool calculated;
    CalculatedState1() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag1, TagT>::value)
            return calculated;
        else
            // there is no previous ObjectOfInterest... Throw an exception or something.
            NoTagSupportException<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag1, TagT>::value)
            calculated = true;
        else
            // there is no previous ObjectOfInterest... Throw an exception or something.
            NoTagSupportException<TagT>();
    }
};
////// end: CalculatedState1


////// begin: CalculatedState2
//    the calculated state for 2 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2>
struct CalculatedState2 {
    // refer back to the 1 parameter CalculatedState
    CalculatedState1<Tag1> prevCalculatedState;
    bool calculated;
    CalculatedState2() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag2, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag2, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState2


////// begin: CalculatedState3
//    the calculated state for 3 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3>
struct CalculatedState3 {
    // refer back to the 2 parameter CalculatedState
    CalculatedState2<Tag1, Tag2> prevCalculatedState;
    bool calculated;
    CalculatedState3() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag3, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag3, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState3


////// begin: CalculatedState4
//    the calculated state for 4 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4>
struct CalculatedState4 {
    // refer back to the 3 parameter CalculatedState
    CalculatedState3<Tag1, Tag2, Tag3> prevCalculatedState;
    bool calculated;
    CalculatedState4() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag4, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag4, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState4


////// begin: CalculatedState5
//    the calculated state for 5 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5>
struct CalculatedState5 {
    // refer back to the 4 parameter CalculatedState
    CalculatedState4<Tag1, Tag2, Tag3, Tag4> prevCalculatedState;
    bool calculated;
    CalculatedState5() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag5, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag5, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState5


////// begin: CalculatedState6
//    the calculated state for 6 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6>
struct CalculatedState6 {
    // refer back to the 5 parameter CalculatedState
    CalculatedState5<Tag1, Tag2, Tag3, Tag4, Tag5> prevCalculatedState;
    bool calculated;
    CalculatedState6() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag6, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag6, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState6


////// begin: CalculatedState7
//    the calculated state for 7 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7>
struct CalculatedState7 {
    // refer back to the 6 parameter CalculatedState
    CalculatedState6<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6> prevCalculatedState;
    bool calculated;
    CalculatedState7() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag7, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag7, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState7


////// begin: CalculatedState8
//    the calculated state for 8 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8>
struct CalculatedState8 {
    // refer back to the 7 parameter CalculatedState
    CalculatedState7<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7> prevCalculatedState;
    bool calculated;
    CalculatedState8() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag8, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag8, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState8


////// begin: CalculatedState9
//    the calculated state for 9 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9>
struct CalculatedState9 {
    // refer back to the 8 parameter CalculatedState
    CalculatedState8<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8> prevCalculatedState;
    bool calculated;
    CalculatedState9() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag9, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag9, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState9


////// begin: CalculatedState10
//    the calculated state for 10 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10>
struct CalculatedState10 {
    // refer back to the 9 parameter CalculatedState
    CalculatedState9<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9> prevCalculatedState;
    bool calculated;
    CalculatedState10() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag10, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag10, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState10


////// begin: CalculatedState11
//    the calculated state for 11 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11>
struct CalculatedState11 {
    // refer back to the 10 parameter CalculatedState
    CalculatedState10<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10> prevCalculatedState;
    bool calculated;
    CalculatedState11() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag11, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag11, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState11


////// begin: CalculatedState12
//    the calculated state for 12 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12>
struct CalculatedState12 {
    // refer back to the 11 parameter CalculatedState
    CalculatedState11<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11> prevCalculatedState;
    bool calculated;
    CalculatedState12() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag12, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag12, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState12


////// begin: CalculatedState13
//    the calculated state for 13 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13>
struct CalculatedState13 {
    // refer back to the 12 parameter CalculatedState
    CalculatedState12<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12> prevCalculatedState;
    bool calculated;
    CalculatedState13() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag13, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag13, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState13


////// begin: CalculatedState14
//    the calculated state for 14 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14>
struct CalculatedState14 {
    // refer back to the 13 parameter CalculatedState
    CalculatedState13<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13> prevCalculatedState;
    bool calculated;
    CalculatedState14() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag14, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag14, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState14


////// begin: CalculatedState15
//    the calculated state for 15 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15>
struct CalculatedState15 {
    // refer back to the 14 parameter CalculatedState
    CalculatedState14<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14> prevCalculatedState;
    bool calculated;
    CalculatedState15() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag15, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag15, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState15


////// begin: CalculatedState16
//    the calculated state for 16 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16>
struct CalculatedState16 {
    // refer back to the 15 parameter CalculatedState
    CalculatedState15<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15> prevCalculatedState;
    bool calculated;
    CalculatedState16() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag16, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag16, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState16


////// begin: CalculatedState17
//    the calculated state for 17 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17>
struct CalculatedState17 {
    // refer back to the 16 parameter CalculatedState
    CalculatedState16<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16> prevCalculatedState;
    bool calculated;
    CalculatedState17() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag17, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag17, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState17


////// begin: CalculatedState18
//    the calculated state for 18 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18>
struct CalculatedState18 {
    // refer back to the 17 parameter CalculatedState
    CalculatedState17<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17> prevCalculatedState;
    bool calculated;
    CalculatedState18() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag18, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag18, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState18


////// begin: CalculatedState19
//    the calculated state for 19 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19>
struct CalculatedState19 {
    // refer back to the 18 parameter CalculatedState
    CalculatedState18<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18> prevCalculatedState;
    bool calculated;
    CalculatedState19() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag19, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag19, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState19


////// begin: CalculatedState20
//    the calculated state for 20 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20>
struct CalculatedState20 {
    // refer back to the 19 parameter CalculatedState
    CalculatedState19<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19> prevCalculatedState;
    bool calculated;
    CalculatedState20() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag20, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag20, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState20


////// begin: CalculatedState21
//    the calculated state for 21 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21>
struct CalculatedState21 {
    // refer back to the 20 parameter CalculatedState
    CalculatedState20<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20> prevCalculatedState;
    bool calculated;
    CalculatedState21() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag21, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag21, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState21


////// begin: CalculatedState22
//    the calculated state for 22 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22>
struct CalculatedState22 {
    // refer back to the 21 parameter CalculatedState
    CalculatedState21<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21> prevCalculatedState;
    bool calculated;
    CalculatedState22() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag22, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag22, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState22


////// begin: CalculatedState23
//    the calculated state for 23 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23>
struct CalculatedState23 {
    // refer back to the 22 parameter CalculatedState
    CalculatedState22<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22> prevCalculatedState;
    bool calculated;
    CalculatedState23() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag23, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag23, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState23


////// begin: CalculatedState24
//    the calculated state for 24 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24>
struct CalculatedState24 {
    // refer back to the 23 parameter CalculatedState
    CalculatedState23<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23> prevCalculatedState;
    bool calculated;
    CalculatedState24() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag24, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag24, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState24


////// begin: CalculatedState25
//    the calculated state for 25 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25>
struct CalculatedState25 {
    // refer back to the 24 parameter CalculatedState
    CalculatedState24<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24> prevCalculatedState;
    bool calculated;
    CalculatedState25() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag25, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag25, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState25


////// begin: CalculatedState26
//    the calculated state for 26 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26>
struct CalculatedState26 {
    // refer back to the 25 parameter CalculatedState
    CalculatedState25<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25> prevCalculatedState;
    bool calculated;
    CalculatedState26() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag26, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag26, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState26


////// begin: CalculatedState27
//    the calculated state for 27 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27>
struct CalculatedState27 {
    // refer back to the 26 parameter CalculatedState
    CalculatedState26<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26> prevCalculatedState;
    bool calculated;
    CalculatedState27() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag27, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag27, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState27


////// begin: CalculatedState28
//    the calculated state for 28 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28>
struct CalculatedState28 {
    // refer back to the 27 parameter CalculatedState
    CalculatedState27<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27> prevCalculatedState;
    bool calculated;
    CalculatedState28() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag28, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag28, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState28


////// begin: CalculatedState29
//    the calculated state for 29 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29>
struct CalculatedState29 {
    // refer back to the 28 parameter CalculatedState
    CalculatedState28<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28> prevCalculatedState;
    bool calculated;
    CalculatedState29() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag29, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag29, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState29


////// begin: CalculatedState30
//    the calculated state for 30 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30>
struct CalculatedState30 {
    // refer back to the 29 parameter CalculatedState
    CalculatedState29<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29> prevCalculatedState;
    bool calculated;
    CalculatedState30() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag30, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag30, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState30


////// begin: CalculatedState31
//    the calculated state for 31 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31>
struct CalculatedState31 {
    // refer back to the 30 parameter CalculatedState
    CalculatedState30<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30> prevCalculatedState;
    bool calculated;
    CalculatedState31() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag31, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag31, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState31


////// begin: CalculatedState32
//    the calculated state for 32 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32>
struct CalculatedState32 {
    // refer back to the 31 parameter CalculatedState
    CalculatedState31<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31> prevCalculatedState;
    bool calculated;
    CalculatedState32() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag32, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag32, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState32


////// begin: CalculatedState33
//    the calculated state for 33 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33>
struct CalculatedState33 {
    // refer back to the 32 parameter CalculatedState
    CalculatedState32<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32> prevCalculatedState;
    bool calculated;
    CalculatedState33() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag33, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag33, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState33


////// begin: CalculatedState34
//    the calculated state for 34 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34>
struct CalculatedState34 {
    // refer back to the 33 parameter CalculatedState
    CalculatedState33<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33> prevCalculatedState;
    bool calculated;
    CalculatedState34() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag34, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag34, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState34


////// begin: CalculatedState35
//    the calculated state for 35 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35>
struct CalculatedState35 {
    // refer back to the 34 parameter CalculatedState
    CalculatedState34<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34> prevCalculatedState;
    bool calculated;
    CalculatedState35() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag35, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag35, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState35


////// begin: CalculatedState36
//    the calculated state for 36 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36>
struct CalculatedState36 {
    // refer back to the 35 parameter CalculatedState
    CalculatedState35<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35> prevCalculatedState;
    bool calculated;
    CalculatedState36() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag36, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag36, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState36


////// begin: CalculatedState37
//    the calculated state for 37 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37>
struct CalculatedState37 {
    // refer back to the 36 parameter CalculatedState
    CalculatedState36<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36> prevCalculatedState;
    bool calculated;
    CalculatedState37() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag37, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag37, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState37


////// begin: CalculatedState38
//    the calculated state for 38 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38>
struct CalculatedState38 {
    // refer back to the 37 parameter CalculatedState
    CalculatedState37<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37> prevCalculatedState;
    bool calculated;
    CalculatedState38() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag38, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag38, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState38


////// begin: CalculatedState39
//    the calculated state for 39 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38, typename Tag39>
struct CalculatedState39 {
    // refer back to the 38 parameter CalculatedState
    CalculatedState38<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38> prevCalculatedState;
    bool calculated;
    CalculatedState39() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag39, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag39, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState39


////// begin: CalculatedState40
//    the calculated state for 40 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38, typename Tag39, typename Tag40>
struct CalculatedState40 {
    // refer back to the 39 parameter CalculatedState
    CalculatedState39<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39> prevCalculatedState;
    bool calculated;
    CalculatedState40() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag40, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag40, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState40


////// begin: CalculatedState41
//    the calculated state for 41 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38, typename Tag39, typename Tag40, 
         typename Tag41>
struct CalculatedState41 {
    // refer back to the 40 parameter CalculatedState
    CalculatedState40<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40> prevCalculatedState;
    bool calculated;
    CalculatedState41() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag41, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag41, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState41


////// begin: CalculatedState42
//    the calculated state for 42 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38, typename Tag39, typename Tag40, 
         typename Tag41, typename Tag42>
struct CalculatedState42 {
    // refer back to the 41 parameter CalculatedState
    CalculatedState41<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41> prevCalculatedState;
    bool calculated;
    CalculatedState42() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag42, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag42, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState42


////// begin: CalculatedState43
//    the calculated state for 43 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38, typename Tag39, typename Tag40, 
         typename Tag41, typename Tag42, typename Tag43>
struct CalculatedState43 {
    // refer back to the 42 parameter CalculatedState
    CalculatedState42<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42> prevCalculatedState;
    bool calculated;
    CalculatedState43() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag43, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag43, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState43


////// begin: CalculatedState44
//    the calculated state for 44 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38, typename Tag39, typename Tag40, 
         typename Tag41, typename Tag42, typename Tag43, typename Tag44>
struct CalculatedState44 {
    // refer back to the 43 parameter CalculatedState
    CalculatedState43<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43> prevCalculatedState;
    bool calculated;
    CalculatedState44() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag44, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag44, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState44


////// begin: CalculatedState45
//    the calculated state for 45 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38, typename Tag39, typename Tag40, 
         typename Tag41, typename Tag42, typename Tag43, typename Tag44, 
         typename Tag45>
struct CalculatedState45 {
    // refer back to the 44 parameter CalculatedState
    CalculatedState44<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43, Tag44> prevCalculatedState;
    bool calculated;
    CalculatedState45() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag45, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag45, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState45


////// begin: CalculatedState46
//    the calculated state for 46 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38, typename Tag39, typename Tag40, 
         typename Tag41, typename Tag42, typename Tag43, typename Tag44, 
         typename Tag45, typename Tag46>
struct CalculatedState46 {
    // refer back to the 45 parameter CalculatedState
    CalculatedState45<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43, Tag44, Tag45> prevCalculatedState;
    bool calculated;
    CalculatedState46() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag46, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag46, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState46


////// begin: CalculatedState47
//    the calculated state for 47 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38, typename Tag39, typename Tag40, 
         typename Tag41, typename Tag42, typename Tag43, typename Tag44, 
         typename Tag45, typename Tag46, typename Tag47>
struct CalculatedState47 {
    // refer back to the 46 parameter CalculatedState
    CalculatedState46<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43, Tag44, Tag45, Tag46> prevCalculatedState;
    bool calculated;
    CalculatedState47() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag47, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag47, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState47


////// begin: CalculatedState48
//    the calculated state for 48 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38, typename Tag39, typename Tag40, 
         typename Tag41, typename Tag42, typename Tag43, typename Tag44, 
         typename Tag45, typename Tag46, typename Tag47, typename Tag48>
struct CalculatedState48 {
    // refer back to the 47 parameter CalculatedState
    CalculatedState47<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43, Tag44, Tag45, Tag46, Tag47> prevCalculatedState;
    bool calculated;
    CalculatedState48() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag48, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag48, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState48


////// begin: CalculatedState49
//    the calculated state for 49 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38, typename Tag39, typename Tag40, 
         typename Tag41, typename Tag42, typename Tag43, typename Tag44, 
         typename Tag45, typename Tag46, typename Tag47, typename Tag48, 
         typename Tag49>
struct CalculatedState49 {
    // refer back to the 48 parameter CalculatedState
    CalculatedState48<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43, Tag44, Tag45, Tag46, Tag47, Tag48> prevCalculatedState;
    bool calculated;
    CalculatedState49() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag49, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag49, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState49


////// begin: CalculatedState50
//    the calculated state for 50 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36, 
         typename Tag37, typename Tag38, typename Tag39, typename Tag40, 
         typename Tag41, typename Tag42, typename Tag43, typename Tag44, 
         typename Tag45, typename Tag46, typename Tag47, typename Tag48, 
         typename Tag49, typename Tag50>
struct CalculatedState50 {
    // refer back to the 49 parameter CalculatedState
    CalculatedState49<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43, Tag44, Tag45, Tag46, Tag47, Tag48, 
                Tag49> prevCalculatedState;
    bool calculated;
    CalculatedState50() : calculated(false) {}  // constructor
    void reset() {  // reset to false
        calculated = false;
        prevCalculatedState.reset();
    }
    
    template<typename TagT> 
    bool isCalculated() {
        // return the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag50, TagT>::value)
            return calculated;
        else
            return prevCalculatedState.isCalculated<TagT>();
    }
    
    template<typename TagT> 
    void setCalculated() {
        // set the CalculatedState of some ObjectOfInterest with a tag
        if(boost::is_same<Tag50, TagT>::value)
            calculated = true;
        else
            prevCalculatedState.setCalculated<TagT>();
    }
};
////// end: CalculatedState50


#endif /* defined(__CalculatedState__) */
