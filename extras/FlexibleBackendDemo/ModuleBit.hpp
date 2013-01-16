#ifndef __ModuleBit__
#define __ModuleBit__

#define __DEBUG__

/* These macros create Bit (or Backend, or whatever) classes, connect those
//     Bits with "Required" and "Provided" Objects Of Interest, and set up 
//     an interface to the Readiness instance which determines whether 
//     or not the Bit is ready to provide. ---Abram  
// NOTE: Once these macros are fully ready, neither the users nor anyone in
//     the collaboration needs to change them. Thus, to see how useful these
//     macros are, please read only the demo files. */

#include <iostream>
#include <typeinfo>
#include "Exceptions.hpp"
#include "Objects.hpp"
#include "Readiness.hpp"

#define CREATE_BIT(BIT) \
/* Create a generic "Bit" which can "Provide" and "Require" Objects Of Interest
//     (OOI) which are specified by "TagType". For now I'll just run with this...
//     I'll put in some placeholder functions which "Report to the core".
//     If the rest of GAMBIT likes this, I can work with the core guys on
//     how to implement those functions later. ---Abram */ \
struct BIT {  \
    /* All of the template functions below need to be specialized for actual OOIs. 
    //     The other macros accomplish this task. */ \
    template <typename TagType> \
    inline void setValue(typename ObjectType<TagType>::type & newValue) { \
        /* Set the value of this OOI (one which is required by this Bit)*/ \
        throw NoTagSupportException<TagType>(); \
    } \
    \
    template <typename TagType> \
    inline typename ObjectType<TagType>::type & getValue() { \
        /* Calculate/Fetch the value of this OOI (one which is provided by this Bit). */ \
        throw NoTagSupportException<TagType>(); \
        return * ObjectOfInterest<TagType>::value; \
    } \
    \
    template <typename TagType> \
    inline bool isRequired() { \
        /* Returns whether or not this Bit requires this OOI<TagType> */ \
        return false; \
    } \
    \
    template <typename TagType> \
    inline bool canBeProvided() { \
        /* Returns whether or not this Bit provides this OOI<TagType> */ \
        return false; \
    } \
    \
    template <typename TagType> \
    void iRequireReport() { \
        /* Report to the core that this Bit requires this OOI<TagType> */ \
    } \
    \
    template <typename TagType> \
    void iProvideReport() { \
        /* Report to the core that this Bit provides this OOI<TagType> */ \
    } \
    \
    inline bool isReady(); /* Declared now, implemented in READY_CALCULATED_BIT */ \
    inline void whatIsNotReady(); /* Declared now, implemented in READY_CALCULATED_BIT */ \
    template <typename TagType> \
    inline void setReady(); /* Declared now, implemented in READY_CALCULATED_BIT */ \
    template <typename TagType> \
    inline bool isCalculated(); /* Declared now, implemented in READY_CALCULATED_BIT */ \
    template <typename TagType> \
    inline void setCalculated(); /* Declared now, implemented in READY_CALCULATED_BIT */ \
    inline void reset(); /* Declared now, implemented in READY_CALCULATED_BIT */ \
/*    inline void calculateAll(); /* Declared now, implemented by the user */ \
};


#define REQUIRE_OBJECT(BIT, TAG_TYPE, READINESS_INSTANCE) \
/* Register that the Bit requires a certain OOI with tag TAG_TYPE. */ \
template <> \
inline void BIT::setValue<Tags::TAG_TYPE>(ObjectType<Tags::TAG_TYPE>::type & newValue) { \
    /* Set the value of this OOI (one which is required by this Bit) */ \
    *ObjectOfInterest<Tags::TAG_TYPE>::value = newValue; \
    if(BIT::isReady()) \
        /* If the bit is ready, setting a value will mean a "new event".
        //     Thus, reset the readiness.*/ \
        BIT::reset(); \
    READINESS_INSTANCE.setReady<Tags::TAG_TYPE>(); \
} \
\
template <> \
inline ObjectType<Tags::TAG_TYPE>::type & BIT::getValue<Tags::TAG_TYPE>() { \
    /* CANNOT BOTH REQUIRE AND PROVIDE AN OOI... */ \
    throw NotProvidedException<Tags::TAG_TYPE>(); \
    return * ObjectOfInterest<Tags::TAG_TYPE>::value; \
} \
\
template <> \
inline bool BIT::isRequired<Tags::TAG_TYPE>() { \
    return true; \
} \
\
template <> \
inline bool BIT::canBeProvided<Tags::TAG_TYPE>() { \
    /* CANNOT BOTH REQUIRE AND PROVIDE AN OOI... */ \
    return false; \
} \
\
template <> \
void BIT::iRequireReport<Tags::TAG_TYPE>() { \
    /* Report to the core that this Bit requires this OOI<TagType> ...
    //     Fake functionality below. Implement later. */ \
    std::cout<<"Dear Core, I require OOI with tag: \n"; \
    std::cout<<Tags::tagName<Tags::TAG_TYPE>()<<"\n"; \
} \
\
template <> \
void BIT::iProvideReport<Tags::TAG_TYPE>() { \
    /* CANNOT BOTH REQUIRE AND PROVIDE AN OOI... */ \
    throw NotProvidedException<Tags::TAG_TYPE>(); \
}



#define PROVIDE_OBJECT(BIT, TAG_TYPE) \
/* Register that the Bit provides a certain OOI with tag TAG_TYPE. */ \
template <> \
inline void BIT::setValue<Tags::TAG_TYPE>(ObjectType<Tags::TAG_TYPE>::type & value) { \
    /* CANNOT BOTH REQUIRE AND PROVIDE AN OOI... */ \
    throw NotRequiredException<Tags::TAG_TYPE>(); \
} \
\
template <> \
inline bool BIT::isRequired<Tags::TAG_TYPE>() { \
    /* CANNOT BOTH REQUIRE AND PROVIDE AN OOI... */ \
    return false; \
} \
\
template <> \
inline bool BIT::canBeProvided<Tags::TAG_TYPE>() { \
    return true; \
} \
\
template <> \
void BIT::iRequireReport<Tags::TAG_TYPE>() { \
    /* CANNOT BOTH REQUIRE AND PROVIDE AN OOI... */ \
    throw NotRequiredException<Tags::TAG_TYPE>(); \
} \
\
template <> \
void BIT::iProvideReport<Tags::TAG_TYPE>() { \
    /* Report to the core that this Bit provides this OOI<TagType> ...
    //     Fake functionality below. Implement later. */ \
    std::cout<<"Dear Core, I provide OOI with tag: \n"; \
    std::cout<<Tags::tagName<Tags::TAG_TYPE>()<<"\n"; \
}



#define READY_CALCULATED_BIT(BIT, READINESS_INSTANCE, CALCULATED_STATE_INSTANCE) \
inline bool BIT::isReady() { \
    /* Connect the Bit with its Readiness class. */ \
    return READINESS_INSTANCE.isReady(); \
} \
\
inline void BIT::whatIsNotReady() { \
    /* Connect the Bit with its Readiness class. */ \
    READINESS_INSTANCE.whatIsNotReady(); \
} \
\
template <typename TagType> \
inline void BIT::setReady() { \
    /* Connect the Bit with its Readiness class. */ \
    READINESS_INSTANCE.setReady<TagType>(); \
} \
\
template <typename TagType> \
inline bool BIT::isCalculated() { \
    /* Connect the Bit with its CalculatedState class. */ \
    return CALCULATED_STATE_INSTANCE.isCalculated<TagType>(); \
} \
\
template <typename TagType> \
inline void BIT::setCalculated() { \
    /* Connect the Bit with its CalculatedState class. */ \
    CALCULATED_STATE_INSTANCE.setCalculated<TagType>(); \
} \
\
inline void BIT::reset() { \
    READINESS_INSTANCE.reset(); \
    CALCULATED_STATE_INSTANCE.reset(); \
} \


#ifdef __DEBUG__
#define SET_CALCULATION(BIT, TAG_TYPE, CALCULATION_BODY) \
/* Set up the getValue method, which calculates the OOI. */ \
template <> \
inline ObjectType<Tags::TAG_TYPE>::type & BIT::getValue<Tags::TAG_TYPE>() { \
    /* Calculate/Fetch the value of this OOI (one which is provided by this Bit). */ \
    if(!BIT::isReady()) { \
        std::cout<<#BIT<<" is not ready. Give it the following requirements with setValue<Tag>(value):\n"; \
        BIT::whatIsNotReady(); \
        throw NotReadyException(); \
    } \
    if(!BIT::isCalculated<Tags::TAG_TYPE>()) { \
        /* The observable needs to be calculated first */ \
        *ObjectOfInterest<Tags::TAG_TYPE>::value = CALCULATION_BODY; \
        BIT::setCalculated<Tags::TAG_TYPE>(); \
    } else \
        std::cout<<"DEBUG: "<<#BIT<<" has already calculated "<<Tags::tagName<Tags::TAG_TYPE>()<<"\n"; \
    return * ObjectOfInterest<Tags::TAG_TYPE>::value; \
} \

#else
#define SET_CALCULATION(BIT, TAG_TYPE, CALCULATION_BODY) \
/* Set up the getValue method, which calculates the OOI. */ \
template <> \
inline ObjectType<Tags::TAG_TYPE>::type & BIT::getValue<Tags::TAG_TYPE>() { \
    /* Calculate/Fetch the value of this OOI (one which is provided by this Bit). */ \
    if(!BIT::isReady()) { \
        std::cout<<#BIT<<" is not ready. Give it the following requirements with setValue<Tag>(value):\n"; \
        BIT::whatIsNotReady(); \
        throw NotReadyException(); \
    } \
    if(!BIT::isCalculated<Tags::TAG_TYPE>()) { \
        /* The observable needs to be calculated first */ \
        *ObjectOfInterest<Tags::TAG_TYPE>::value = CALCULATION_BODY; \
        BIT::setCalculated<Tags::TAG_TYPE>(); \
    } \
    \
    return * ObjectOfInterest<Tags::TAG_TYPE>::value; \
} \

#endif /* defined(__DEBUG__) */


#endif /* defined(__ModuleBit__) */



