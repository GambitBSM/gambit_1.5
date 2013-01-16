#include <iostream>
#include <math.h>
#include "Objects.hpp"
#include "ModuleBit.hpp"
#include "Readiness.hpp"
#include "CalculatedState.hpp"

/* Not everyone was convinced with my last demo, so I decided to improve it. 
//     Please focus on these "demo" files. They show the true power of all 
//     those macros I wrote. If you are worried about how complicated those
//     macros are, don't be worried; neither the collaboration members nor the 
//     users will ever have to modify them once we have them finalized.
//     Also, here are a lot more comments for your viewing pleasure. ---Abram */

using namespace std;


/* **** OBJECT OF INTEREST (OOI) SET UP ****
// These macros set up the Objects Of Interest (OOIs). They associate tags with 
//     variable types. Here, all the variable types are doubles. However, note
//     that these macros would be *extremely* useful for handling data of 
//     much different types (to be demonstrated in a future example code). */

// One fake OOI, for a debug test:
ASSOCIATE_TAG_OBJECT(HeadsWillRoll, int)
// One Tag which is not even associated with an OOI, for a debug test:
namespace Tags {
    struct MarblesRollToo{};
}

// OOIs which will be required by the Backend, set up by macros:
ASSOCIATE_TAG_OBJECT(Px, double)
ASSOCIATE_TAG_OBJECT(Py, double)
ASSOCIATE_TAG_OBJECT(Pz, double)
ASSOCIATE_TAG_OBJECT(P0, double)
// Instantiate a Readiness class to keep track of when the Backend has all 
//     of its requirements met:
Readiness4<Tags::Px, Tags::Py, Tags::Pz, Tags::P0> readiness;

// OOIs which will be provided by the Backend, set up by macros:
ASSOCIATE_TAG_OBJECT(PT, double)
ASSOCIATE_TAG_OBJECT(Mass, double)
ASSOCIATE_TAG_OBJECT(Theta, double)
ASSOCIATE_TAG_OBJECT(Phi, double)
// Instantiate a CalculatedState class to keep track of whether
//     or not the Backend has already calculated a certain observable:
CalculatedState4<Tags::PT, Tags::Mass, Tags::Theta, Tags::Phi> calculatedState;



// **** BIT/BACKEND SET UP ****

// This macro easily creates the Bit/Backend.
CREATE_BIT(KinematicsCalculator)

// These macros state that the Backend we created requires some of these
//     OOIs in order to provide anything in return. The readiness class is 
//     interfaced for each OOI as well.
REQUIRE_OBJECT(KinematicsCalculator, Px, readiness)
REQUIRE_OBJECT(KinematicsCalculator, Py, readiness)
REQUIRE_OBJECT(KinematicsCalculator, Pz, readiness)
REQUIRE_OBJECT(KinematicsCalculator, P0, readiness)

// This macro finalizes the interface between the readiness class, the
//     calćulated state class, and the Backend.
READY_CALCULATED_BIT(KinematicsCalculator, readiness, calculatedState)

// These macros state that the Backend can provide some other OOIs.
PROVIDE_OBJECT(KinematicsCalculator, PT)
PROVIDE_OBJECT(KinematicsCalculator, Mass)
PROVIDE_OBJECT(KinematicsCalculator, Theta)
PROVIDE_OBJECT(KinematicsCalculator, Phi)



// IMPORTANT: Up until this point, the user has not had to hard code anything.
//     The user has simply used all of our convenient macros and classes.
//     In fact, there is even a macro for supplying a simple calculation:
SET_CALCULATION(KinematicsCalculator, PT, sqrt(OOI(Px) * OOI(Px) + OOI(Py) * OOI(Py)) )
// NOTE: The OOI macro gives the value of the ObjectOfInterest.


// The user can also specify calculations with some custom methods:
inline double calculateMass(double inPx, double inPy, double inPz, double inP0) {
    return sqrt(inP0*inP0 - inPx*inPx - inPy*inPy - inPz*inPz);
}
SET_CALCULATION(KinematicsCalculator, Mass, calculateMass(OOI(Px), OOI(Py), OOI(Pz), OOI(P0)) )


// ... Which is very handy if the calculations are complicated:
inline double calculateTheta(double inPt, double inPz) {
    double tempTheta;
    double PI = acos(-1.);
    
    tempTheta = atan(fabs(inPt / inPz));
    if (inPz < 0.)
        tempTheta = PI - tempTheta;
    
    return tempTheta;
}
inline double calculateTheta(double inPx, double inPy, double inPz) {
    double tempPt = sqrt(inPx*inPx + inPy*inPy);
    return calculateTheta(tempPt, inPz);
}
SET_CALCULATION(KinematicsCalculator, Theta, calculateTheta(OOI(Px), OOI(Py), OOI(Pz)) )


inline double calculatePhi(double inPx, double inPy) {
    double tempPhi;
    double PI = acos(-1.);
    
    tempPhi = atan(fabs(inPy / inPx));
    if(inPx < 0. && inPy > 0.)
        tempPhi = PI - tempPhi;
    if(inPx < 0. && inPy < 0.)
        tempPhi = PI + tempPhi;
    if(inPx > 0. && inPy < 0.)
        tempPhi = 2*PI - tempPhi;

  return tempPhi;
}
SET_CALCULATION(KinematicsCalculator, Phi, calculatePhi(OOI(Px), OOI(Py)) )


int main() {
    // For the purposes of this demonstration, think of the main() program
    //     as a toy Core.
    
    // Let's say for now that the Core controls all the actual variables.
    double px, py, pz, p0, pT, mass, theta, phi;
    double inputPx = 6.;
    double inputPy = 8.;
    double inputPz = 30.;
    double inputP0 = 100.;
    int inputCrap = 45;
    double temp;
    
    // These macros pass the addresses of those variables to the OOIs:
    PASS_ADDRESS(Px, &px)
    PASS_ADDRESS(Py, &py)
    PASS_ADDRESS(Pz, &pz)
    PASS_ADDRESS(P0, &p0)
    PASS_ADDRESS(PT, &pT)
    PASS_ADDRESS(Mass, &mass)
    PASS_ADDRESS(Theta, &theta)
    PASS_ADDRESS(Phi, &phi)
    
    // Instantiate that Backend we created:
    KinematicsCalculator kinematicsBackend;
    
    cout<<"\nMost initializations are performed with easy to use macros.\n";
    cout<<"    Please look at the source: 'demoKinematicsCalculator.cpp'\n";
    
    cout<<"\nNow testing the functionality and Exceptions for this Backend...\n\n";
    
    cout<<"\nFirst, let's try to retrieve some value,\n";
    cout<<"    although KinematicsCalculator is not ready.\n";
    cout<<"  TRY: temp = kinematicsBackend.getValue<Tags::PT>();\n";
    try {
        temp = kinematicsBackend.getValue<Tags::PT>();
    }
    catch (exception &e)
    {
        cout<<"****Exception caught\n";
        cout<<e.what()<<"\n";
    }
    
    cout<<"\nLet's try to provide a value for some OOI\n";
    cout<<"    which KinematicsCalculator cannot handle.\n";
    cout<<"  TRY: kinematicsBackend.setValue<Tags::HeadsWillRoll>(inputCrap);\n";
    try {
        kinematicsBackend.setValue<Tags::HeadsWillRoll>(inputCrap);
    }
    catch (exception &e)
    {
        cout<<"****Exception caught\n";
        cout<<e.what()<<"\n";
    }
    
    cout<<"\nLet's try to provide a value for some OOI\n";
    cout<<"    which is not even initialized at all.\n";
    cout<<"  TRY: kinematicsBackend.setValue<Tags::MarblesRollToo>(inputCrap);\n";
    cout<<"  OOPS: That gives a compile-time error.\n";
    cout<<"  TRY instead: kinematicsBackend.getValue<Tags::MarblesRollToo>();\n";
     try {
        kinematicsBackend.getValue<Tags::MarblesRollToo>();
    }
    catch (exception &e)
    {
        cout<<"****Exception caught\n";
        cout<<e.what()<<"\n";
    }
    
    cout<<"\nLet's provide some required OOIs\n";
    cout<<"  TRY: kinematicsBackend.setValue<Tags::Px>(inputPx);\n";
    cout<<"  TRY: kinematicsBackend.setValue<Tags::Pz>(inputPz);\n";
    try {
        kinematicsBackend.setValue<Tags::Px>(inputPx);
        kinematicsBackend.setValue<Tags::Pz>(inputPz);
    }
    catch (exception &e)
    {
        cout<<"****Exception caught\n";
        cout<<e.what()<<"\n";
    }
    
    cout<<"\nAgain, let's try to retrieve some value,\n";
    cout<<"    although KinematicsCalculator is not ready.\n";
    cout<<"  TRY: temp = kinematicsBackend.getValue<Tags::PT>();\n";
    try {
        temp = kinematicsBackend.getValue<Tags::PT>();
    }
    catch (exception &e)
    {
        cout<<"****Exception caught\n";
        cout<<e.what()<<"\n";
    }
    
    cout<<"\nLet's provide the rest of the required OOIs\n";
    cout<<"  TRY: kinematicsBackend.setValue<Tags::Py>(inputPy);\n";
    cout<<"  TRY: kinematicsBackend.setValue<Tags::P0>(inputP0);\n";
    try {
        kinematicsBackend.setValue<Tags::Py>(inputPy);
        kinematicsBackend.setValue<Tags::P0>(inputP0);
    }
    catch (exception &e)
    {
        cout<<"****Exception caught\n";
        cout<<e.what()<<"\n";
    }
    
    cout<<"\nFinally, let's try to retrieve some value,\n";
    cout<<"    now that KinematicsCalculator is ready.\n";
    cout<<"  TRY: temp = kinematicsBackend.getValue<Tags::PT>();\n";
    try {
        temp = kinematicsBackend.getValue<Tags::PT>();
    }
    catch (exception &e)
    {
        cout<<"\nException caught\n";
        cout<<e.what()<<"\n";
    }
    cout<<"temp is "<<temp<<"\n";
    cout<<"temp should be "<<sqrt(inputPx*inputPx + inputPy*inputPy)<<"\n";
    
    cout<<"\nKinematicsCalculator can detect if it has already calculated,\n";
    cout<<"    some value, to not waste time by calculating it again.\n";
    cout<<"  TRY: temp = kinematicsBackend.getValue<Tags::PT>();\n";
    try {
        temp = kinematicsBackend.getValue<Tags::PT>();
    }
    catch (exception &e)
    {
        cout<<"\nException caught\n";
        cout<<e.what()<<"\n";
    }
    cout<<"temp is "<<temp<<"\n";
    cout<<"temp should be "<<sqrt(inputPx*inputPx + inputPy*inputPy)<<"\n";
    
    cout<<"\nIf we provide one new OOI as input, KinematicsCalculator\n";
    cout<<"    considers it a 'new event'. \n";
    cout<<"  TRY: inputPy = 0.6;\n";
    cout<<"  AND: kinematicsBackend.setValue<Tags::Py>(inputPy);\n";
    cout<<"  AND: temp = kinematicsBackend.getValue<Tags::PT>();\n";
    try {
        inputPy = 0.6;
        kinematicsBackend.setValue<Tags::Py>(inputPy);
        temp = kinematicsBackend.getValue<Tags::PT>();
    }
    catch (exception &e)
    {
        cout<<"****Exception caught\n";
        cout<<e.what()<<"\n";
    }
    
    cout<<"\nThus, the remaining requirements must also be set. \n";
    cout<<"  TRY: inputPx = 0.8; inputPz = 1.2; inputP0 = 4.6;\n";
    cout<<"  AND: kinematicsBackend.setValue<Tags::Px>(inputPx);\n";
    cout<<"  AND: kinematicsBackend.setValue<Tags::Pz>(inputPz);\n";
    cout<<"  AND: kinematicsBackend.setValue<Tags::P0>(inputP0);\n";
    cout<<"  AND: temp = kinematicsBackend.getValue<Tags::PT>();\n";
    try {
        inputPx = 0.8; inputPz = 1.2; inputP0 = 4.6;
        kinematicsBackend.setValue<Tags::Px>(inputPx);
        kinematicsBackend.setValue<Tags::Pz>(inputPz);
        kinematicsBackend.setValue<Tags::P0>(inputP0);
        temp = kinematicsBackend.getValue<Tags::PT>();
    }
    catch (exception &e)
    {
        cout<<"****Exception caught\n";
        cout<<e.what()<<"\n";
    }
    cout<<"temp is "<<temp<<"\n";
    cout<<"temp should be "<<sqrt(inputPx*inputPx + inputPy*inputPy)<<"\n";
    cout<<"NOTE: Work in progress: Some requirements may not need to be reset\n";
    cout<<"      for each new event. This is easy enough to implement.\n";
    
    cout<<"\nFinally, let's test the remaining calculable observables: \n";
    cout<<"  TRY: temp = kinematicsBackend.getValue<Tags::PT>();\n";
    cout<<"  AND: temp = kinematicsBackend.getValue<Tags::Mass>();\n";
    cout<<"  AND: temp = kinematicsBackend.getValue<Tags::Theta>();\n";
    cout<<"  AND: temp = kinematicsBackend.getValue<Tags::Phi>();\n";
    try {
        temp = kinematicsBackend.getValue<Tags::PT>();
        cout<<"PT is "<<temp<<", and should be "<<sqrt(inputPx*inputPx + inputPy*inputPy)<<"\n";
        temp = kinematicsBackend.getValue<Tags::Mass>();
        cout<<"Mass is "<<temp<<", and should be "<<calculateMass(inputPx, inputPy, inputPz, inputP0)<<"\n";
        temp = kinematicsBackend.getValue<Tags::Theta>();
        cout<<"Theta is "<<temp<<", and should be "<<calculateTheta(inputPx, inputPy, inputPz)<<"\n";
        temp = kinematicsBackend.getValue<Tags::Phi>();
        cout<<"Phi is "<<temp<<", and should be "<<calculatePhi(inputPx, inputPy)<<"\n";
    }
    catch (exception &e)
    {
        cout<<"****Exception caught\n";
        cout<<e.what()<<"\n";
    }
    
    cout<<"\n\n   End of demo. Have a nice day!\n\n";
}



