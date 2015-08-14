//==========================================================================
// This file has been automatically generated for Pythia 8
// MadGraph5_aMC@NLO v. 2.3.0, 2015-07-01
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Pythia8_Sigma_MC4BSM_2012_UFO_qq_evevx_H
#define Pythia8_Sigma_MC4BSM_2012_UFO_qq_evevx_H

#include <complex> 

#include "Pythia8/SigmaProcess.h"
#include "Parameters_MC4BSM_2012_UFO.h"

using namespace std; 

namespace Pythia8 
{
//==========================================================================
// A class for calculating the matrix elements for
// Process: u u~ > ev ev~ WEIGHTED=4
// Process: c c~ > ev ev~ WEIGHTED=4
// Process: d d~ > ev ev~ WEIGHTED=4
// Process: s s~ > ev ev~ WEIGHTED=4
//--------------------------------------------------------------------------

class Sigma_MC4BSM_2012_UFO_qq_evevx : public Sigma2Process 
{
  public:

    // Constructor.
    Sigma_MC4BSM_2012_UFO_qq_evevx() {}

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
    virtual string name() const {return "qq_evevx (MC4BSM_2012_UFO)";}

    virtual int code() const {return 10002;}

    virtual string inFlux() const {return "qqbarSame";}
    int id3Mass() const {return 9000009;}
    int id4Mass() const {return 9000009;}
    virtual int resonanceA() const {return 23;}
    virtual bool isSChannel() const {return true;}
    // Tell Pythia that sigmaHat returns the ME^2
    virtual bool convertM2() const {return true;}

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    static const int nwavefuncs = 8; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 4; 
    std::complex<double> amp[namplitudes]; 
    double matrix_uux_evevx(); 
    double matrix_ddx_evevx(); 

    // Constants for array limits
    static const int nexternal = 4; 
    static const int nprocesses = 4; 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // Pointer to the model parameters
    Parameters_MC4BSM_2012_UFO * pars; 

}; 

}  // end namespace Pythia8

#endif  // Pythia8_Sigma_MC4BSM_2012_UFO_qq_evevx_H

