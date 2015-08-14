//==========================================================================
// This file has been automatically generated for Pythia 8
// MadGraph5_aMC@NLO v. 2.3.0, 2015-07-01
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Pythia8_Sigma_MC4BSM_2012_UFO_qq_p1p1_H
#define Pythia8_Sigma_MC4BSM_2012_UFO_qq_p1p1_H

#include <complex> 

#include "Pythia8/SigmaProcess.h"
#include "Parameters_MC4BSM_2012_UFO.h"

using namespace std; 

namespace Pythia8 
{
//==========================================================================
// A class for calculating the matrix elements for
// Process: u u~ > p1 p1 WEIGHTED=2
//--------------------------------------------------------------------------

class Sigma_MC4BSM_2012_UFO_qq_p1p1 : public Sigma2Process 
{
  public:

    // Constructor.
    Sigma_MC4BSM_2012_UFO_qq_p1p1() {}

    // Initialize process.
    virtual void initProc(); 

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin(); 

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat(); 

    // Select flavour, colour and anticolour.
    virtual void setIdColAcol(); 

    // Evaluate weight for decay angles.
    virtual double weightDecay(Event& process, int iResBeg, int iResEnd); 

    // Info on the subprocess.
    virtual string name() const {return "qq_p1p1 (MC4BSM_2012_UFO)";}

    virtual int code() const {return 10005;}

    virtual string inFlux() const {return "qqbarSame";}
    int id3Mass() const {return 9000006;}
    int id4Mass() const {return 9000006;}

    // Tell Pythia that sigmaHat returns the ME^2
    virtual bool convertM2() const {return true;}

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    static const int nwavefuncs = 6; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 2; 
    std::complex<double> amp[namplitudes]; 
    double matrix_uux_p1p1(); 

    // Constants for array limits
    static const int nexternal = 4; 
    static const int nprocesses = 2; 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // Pointer to the model parameters
    Parameters_MC4BSM_2012_UFO * pars; 

}; 

}  // end namespace Pythia8

#endif  // Pythia8_Sigma_MC4BSM_2012_UFO_qq_p1p1_H

