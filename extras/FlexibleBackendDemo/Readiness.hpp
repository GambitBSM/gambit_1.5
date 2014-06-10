#ifndef __Readiness__
#define __Readiness__

/* These Readiness classes are used by the modules, backends, etc. 
//     They provide an easy way to check if the module/backend/whatever
//     has everything it needs to execute. ---Abram  
// NOTE: Once these classes are fully ready, neither the users nor anyone in
//     the collaboration needs to change them. Thus, to see how useful these
//     classes are, please read only the demo files. */

#include <iostream>
#include <boost/type_traits/is_same.hpp>
#include "Exceptions.hpp"
#include "Objects.hpp"

using namespace std;

////// begin: Readiness1
//    the state of readiness for a single ObjectOfInterest with a tag 
template<typename Tag1>
struct Readiness1 {
    bool ready;
    Readiness1() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
    }
    
    bool isReady() {  // return the value
        return ready;
    }
    
    void whatIsNotReady() {  // return the value
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag1>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag1, TagT>::value)
            ready = true;
        else
            // there is no previous ObjectOfInterest... Throw an exception or something.
            throw NoTagSupportException<TagT>();
    }
};
////// end: Readiness1


////// begin: Readiness2
//    the state of readiness for 2 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2>
struct Readiness2 {
    // refer back to the 1 parameter Readiness
    Readiness1<Tag1> prevReadiness;
    bool ready;
    Readiness2() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag2>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag2, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness2


////// begin: Readiness3
//    the state of readiness for 3 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3>
struct Readiness3 {
    // refer back to the 2 parameter Readiness
    Readiness2<Tag1, Tag2> prevReadiness;
    bool ready;
    Readiness3() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag3>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag3, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness3


////// begin: Readiness4
//    the state of readiness for 4 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4>
struct Readiness4 {
    // refer back to the 3 parameter Readiness
    Readiness3<Tag1, Tag2, Tag3> prevReadiness;
    bool ready;
    Readiness4() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag4>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag4, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness4


////// begin: Readiness5
//    the state of readiness for 5 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5>
struct Readiness5 {
    // refer back to the 4 parameter Readiness
    Readiness4<Tag1, Tag2, Tag3, Tag4> prevReadiness;
    bool ready;
    Readiness5() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag5>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag5, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness5


////// begin: Readiness6
//    the state of readiness for 6 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6>
struct Readiness6 {
    // refer back to the 5 parameter Readiness
    Readiness5<Tag1, Tag2, Tag3, Tag4, Tag5> prevReadiness;
    bool ready;
    Readiness6() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag6>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag6, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness6


////// begin: Readiness7
//    the state of readiness for 7 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7>
struct Readiness7 {
    // refer back to the 6 parameter Readiness
    Readiness6<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6> prevReadiness;
    bool ready;
    Readiness7() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag7>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag7, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness7


////// begin: Readiness8
//    the state of readiness for 8 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8>
struct Readiness8 {
    // refer back to the 7 parameter Readiness
    Readiness7<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7> prevReadiness;
    bool ready;
    Readiness8() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag8>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag8, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness8


////// begin: Readiness9
//    the state of readiness for 9 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9>
struct Readiness9 {
    // refer back to the 8 parameter Readiness
    Readiness8<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8> prevReadiness;
    bool ready;
    Readiness9() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag9>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag9, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness9


////// begin: Readiness10
//    the state of readiness for 10 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10>
struct Readiness10 {
    // refer back to the 9 parameter Readiness
    Readiness9<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9> prevReadiness;
    bool ready;
    Readiness10() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag10>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag10, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness10


////// begin: Readiness11
//    the state of readiness for 11 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11>
struct Readiness11 {
    // refer back to the 10 parameter Readiness
    Readiness10<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10> prevReadiness;
    bool ready;
    Readiness11() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag11>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag11, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness11


////// begin: Readiness12
//    the state of readiness for 12 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12>
struct Readiness12 {
    // refer back to the 11 parameter Readiness
    Readiness11<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11> prevReadiness;
    bool ready;
    Readiness12() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag12>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag12, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness12


////// begin: Readiness13
//    the state of readiness for 13 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13>
struct Readiness13 {
    // refer back to the 12 parameter Readiness
    Readiness12<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12> prevReadiness;
    bool ready;
    Readiness13() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag13>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag13, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness13


////// begin: Readiness14
//    the state of readiness for 14 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14>
struct Readiness14 {
    // refer back to the 13 parameter Readiness
    Readiness13<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13> prevReadiness;
    bool ready;
    Readiness14() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag14>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag14, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness14


////// begin: Readiness15
//    the state of readiness for 15 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15>
struct Readiness15 {
    // refer back to the 14 parameter Readiness
    Readiness14<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14> prevReadiness;
    bool ready;
    Readiness15() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag15>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag15, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness15


////// begin: Readiness16
//    the state of readiness for 16 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16>
struct Readiness16 {
    // refer back to the 15 parameter Readiness
    Readiness15<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15> prevReadiness;
    bool ready;
    Readiness16() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag16>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag16, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness16


////// begin: Readiness17
//    the state of readiness for 17 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17>
struct Readiness17 {
    // refer back to the 16 parameter Readiness
    Readiness16<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16> prevReadiness;
    bool ready;
    Readiness17() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag17>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag17, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness17


////// begin: Readiness18
//    the state of readiness for 18 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18>
struct Readiness18 {
    // refer back to the 17 parameter Readiness
    Readiness17<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17> prevReadiness;
    bool ready;
    Readiness18() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag18>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag18, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness18


////// begin: Readiness19
//    the state of readiness for 19 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19>
struct Readiness19 {
    // refer back to the 18 parameter Readiness
    Readiness18<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18> prevReadiness;
    bool ready;
    Readiness19() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag19>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag19, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness19


////// begin: Readiness20
//    the state of readiness for 20 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20>
struct Readiness20 {
    // refer back to the 19 parameter Readiness
    Readiness19<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19> prevReadiness;
    bool ready;
    Readiness20() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag20>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag20, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness20


////// begin: Readiness21
//    the state of readiness for 21 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21>
struct Readiness21 {
    // refer back to the 20 parameter Readiness
    Readiness20<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20> prevReadiness;
    bool ready;
    Readiness21() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag21>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag21, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness21


////// begin: Readiness22
//    the state of readiness for 22 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22>
struct Readiness22 {
    // refer back to the 21 parameter Readiness
    Readiness21<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21> prevReadiness;
    bool ready;
    Readiness22() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag22>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag22, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness22


////// begin: Readiness23
//    the state of readiness for 23 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23>
struct Readiness23 {
    // refer back to the 22 parameter Readiness
    Readiness22<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22> prevReadiness;
    bool ready;
    Readiness23() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag23>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag23, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness23


////// begin: Readiness24
//    the state of readiness for 24 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24>
struct Readiness24 {
    // refer back to the 23 parameter Readiness
    Readiness23<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23> prevReadiness;
    bool ready;
    Readiness24() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag24>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag24, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness24


////// begin: Readiness25
//    the state of readiness for 25 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25>
struct Readiness25 {
    // refer back to the 24 parameter Readiness
    Readiness24<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24> prevReadiness;
    bool ready;
    Readiness25() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag25>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag25, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness25


////// begin: Readiness26
//    the state of readiness for 26 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26>
struct Readiness26 {
    // refer back to the 25 parameter Readiness
    Readiness25<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25> prevReadiness;
    bool ready;
    Readiness26() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag26>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag26, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness26


////// begin: Readiness27
//    the state of readiness for 27 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27>
struct Readiness27 {
    // refer back to the 26 parameter Readiness
    Readiness26<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26> prevReadiness;
    bool ready;
    Readiness27() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag27>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag27, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness27


////// begin: Readiness28
//    the state of readiness for 28 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28>
struct Readiness28 {
    // refer back to the 27 parameter Readiness
    Readiness27<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27> prevReadiness;
    bool ready;
    Readiness28() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag28>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag28, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness28


////// begin: Readiness29
//    the state of readiness for 29 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29>
struct Readiness29 {
    // refer back to the 28 parameter Readiness
    Readiness28<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28> prevReadiness;
    bool ready;
    Readiness29() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag29>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag29, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness29


////// begin: Readiness30
//    the state of readiness for 30 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30>
struct Readiness30 {
    // refer back to the 29 parameter Readiness
    Readiness29<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29> prevReadiness;
    bool ready;
    Readiness30() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag30>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag30, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness30


////// begin: Readiness31
//    the state of readiness for 31 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31>
struct Readiness31 {
    // refer back to the 30 parameter Readiness
    Readiness30<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30> prevReadiness;
    bool ready;
    Readiness31() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag31>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag31, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness31


////// begin: Readiness32
//    the state of readiness for 32 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32>
struct Readiness32 {
    // refer back to the 31 parameter Readiness
    Readiness31<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31> prevReadiness;
    bool ready;
    Readiness32() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag32>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag32, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness32


////// begin: Readiness33
//    the state of readiness for 33 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33>
struct Readiness33 {
    // refer back to the 32 parameter Readiness
    Readiness32<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32> prevReadiness;
    bool ready;
    Readiness33() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag33>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag33, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness33


////// begin: Readiness34
//    the state of readiness for 34 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34>
struct Readiness34 {
    // refer back to the 33 parameter Readiness
    Readiness33<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33> prevReadiness;
    bool ready;
    Readiness34() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag34>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag34, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness34


////// begin: Readiness35
//    the state of readiness for 35 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35>
struct Readiness35 {
    // refer back to the 34 parameter Readiness
    Readiness34<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34> prevReadiness;
    bool ready;
    Readiness35() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag35>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag35, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness35


////// begin: Readiness36
//    the state of readiness for 36 ObjectOfInterests with tags 
template<typename Tag1, typename Tag2, typename Tag3, typename Tag4, 
         typename Tag5, typename Tag6, typename Tag7, typename Tag8, 
         typename Tag9, typename Tag10, typename Tag11, typename Tag12, 
         typename Tag13, typename Tag14, typename Tag15, typename Tag16, 
         typename Tag17, typename Tag18, typename Tag19, typename Tag20, 
         typename Tag21, typename Tag22, typename Tag23, typename Tag24, 
         typename Tag25, typename Tag26, typename Tag27, typename Tag28, 
         typename Tag29, typename Tag30, typename Tag31, typename Tag32, 
         typename Tag33, typename Tag34, typename Tag35, typename Tag36>
struct Readiness36 {
    // refer back to the 35 parameter Readiness
    Readiness35<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35> prevReadiness;
    bool ready;
    Readiness36() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag36>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag36, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness36


////// begin: Readiness37
//    the state of readiness for 37 ObjectOfInterests with tags 
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
struct Readiness37 {
    // refer back to the 36 parameter Readiness
    Readiness36<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36> prevReadiness;
    bool ready;
    Readiness37() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag37>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag37, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness37


////// begin: Readiness38
//    the state of readiness for 38 ObjectOfInterests with tags 
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
struct Readiness38 {
    // refer back to the 37 parameter Readiness
    Readiness37<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37> prevReadiness;
    bool ready;
    Readiness38() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag38>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag38, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness38


////// begin: Readiness39
//    the state of readiness for 39 ObjectOfInterests with tags 
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
struct Readiness39 {
    // refer back to the 38 parameter Readiness
    Readiness38<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38> prevReadiness;
    bool ready;
    Readiness39() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag39>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag39, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness39


////// begin: Readiness40
//    the state of readiness for 40 ObjectOfInterests with tags 
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
struct Readiness40 {
    // refer back to the 39 parameter Readiness
    Readiness39<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39> prevReadiness;
    bool ready;
    Readiness40() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag40>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag40, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness40


////// begin: Readiness41
//    the state of readiness for 41 ObjectOfInterests with tags 
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
struct Readiness41 {
    // refer back to the 40 parameter Readiness
    Readiness40<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40> prevReadiness;
    bool ready;
    Readiness41() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag41>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag41, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness41


////// begin: Readiness42
//    the state of readiness for 42 ObjectOfInterests with tags 
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
struct Readiness42 {
    // refer back to the 41 parameter Readiness
    Readiness41<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41> prevReadiness;
    bool ready;
    Readiness42() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag42>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag42, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness42


////// begin: Readiness43
//    the state of readiness for 43 ObjectOfInterests with tags 
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
struct Readiness43 {
    // refer back to the 42 parameter Readiness
    Readiness42<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42> prevReadiness;
    bool ready;
    Readiness43() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag43>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag43, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness43


////// begin: Readiness44
//    the state of readiness for 44 ObjectOfInterests with tags 
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
struct Readiness44 {
    // refer back to the 43 parameter Readiness
    Readiness43<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43> prevReadiness;
    bool ready;
    Readiness44() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag44>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag44, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness44


////// begin: Readiness45
//    the state of readiness for 45 ObjectOfInterests with tags 
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
struct Readiness45 {
    // refer back to the 44 parameter Readiness
    Readiness44<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43, Tag44> prevReadiness;
    bool ready;
    Readiness45() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag45>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag45, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness45


////// begin: Readiness46
//    the state of readiness for 46 ObjectOfInterests with tags 
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
struct Readiness46 {
    // refer back to the 45 parameter Readiness
    Readiness45<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43, Tag44, Tag45> prevReadiness;
    bool ready;
    Readiness46() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag46>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag46, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness46


////// begin: Readiness47
//    the state of readiness for 47 ObjectOfInterests with tags 
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
struct Readiness47 {
    // refer back to the 46 parameter Readiness
    Readiness46<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43, Tag44, Tag45, Tag46> prevReadiness;
    bool ready;
    Readiness47() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag47>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag47, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness47


////// begin: Readiness48
//    the state of readiness for 48 ObjectOfInterests with tags 
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
struct Readiness48 {
    // refer back to the 47 parameter Readiness
    Readiness47<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43, Tag44, Tag45, Tag46, Tag47> prevReadiness;
    bool ready;
    Readiness48() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag48>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag48, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness48


////// begin: Readiness49
//    the state of readiness for 49 ObjectOfInterests with tags 
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
struct Readiness49 {
    // refer back to the 48 parameter Readiness
    Readiness48<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43, Tag44, Tag45, Tag46, Tag47, Tag48> prevReadiness;
    bool ready;
    Readiness49() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag49>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag49, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness49


////// begin: Readiness50
//    the state of readiness for 50 ObjectOfInterests with tags 
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
struct Readiness50 {
    // refer back to the 49 parameter Readiness
    Readiness49<Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, 
                Tag7, Tag8, Tag9, Tag10, Tag11, Tag12, 
                Tag13, Tag14, Tag15, Tag16, Tag17, Tag18, 
                Tag19, Tag20, Tag21, Tag22, Tag23, Tag24, 
                Tag25, Tag26, Tag27, Tag28, Tag29, Tag30, 
                Tag31, Tag32, Tag33, Tag34, Tag35, Tag36, 
                Tag37, Tag38, Tag39, Tag40, Tag41, Tag42, 
                Tag43, Tag44, Tag45, Tag46, Tag47, Tag48, 
                Tag49> prevReadiness;
    bool ready;
    Readiness50() : ready(false) {}  // constructor
    void reset() {  // reset to false
        ready = false;
        prevReadiness.reset();
    }
    
    bool isReady() {  // return the value
        return (ready && prevReadiness.isReady());
    }
    
    void whatIsNotReady() {  // return the value
        prevReadiness.whatIsNotReady();
        if (!ready)
            std::cout<<"  "<<Tags::tagName<Tag50>()<<"\n";
    }
    
    template<typename TagT> 
    void setReady() {
        // set the readiness of some ObjectOfInterest with a tag
        if(boost::is_same<Tag50, TagT>::value)
            ready = true;
        else
            prevReadiness.setReady<TagT>();
    }
};
////// end: Readiness50


#endif /* defined(__Readiness__) */
