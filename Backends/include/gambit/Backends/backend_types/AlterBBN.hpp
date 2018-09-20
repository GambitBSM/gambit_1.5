//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Definitions of the container structure
///  for the AlterBBN backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 Jun
///  
///
///  *********************************************


#ifndef __AlterBBN_types_hpp__
#define __AlterBBN_types_hpp__

#define NNUCREAC 100
#define NTABMAX 1000

namespace Gambit
{

  struct relicparam
/* structure containing the cosmological model parameters */
  {
    int entropy_model,energy_model;
    double dd0,ndd,Tdend,Tddeq; // dark density
    double sd0,nsd,Tsend; // dark entropy
    double Sigmad0,nSigmad,TSigmadend; // dark entropy injection
    double Sigmarad0,nSigmarad,TSigmaradend; // standard entropy injection
    double nt0,nnt,Tnend; // non-thermal production of relics
    
    double quintn2,quintn3,quintn4,quintT12,quintT23,quintT34; // effective quintessence model

    int phi_model; // decaying scalar field model switch
    double eta_phi,Gamma_phi,rhot_phi_Tmax,rho_phi; // eta_phi = b / m_phi
    double rhot_phi0,Tphi0;
    double T_RH;
    double Sigmatildestar;
    double Sigmatildestar_max;
    double Tstdstar_max;

    double mgravitino; // gravitino mass
            
    double relicmass;
    int scalar;
    
    int solver; // switch for linear or logarithmic differential equation solver
    
    double T; // Temperature in GeV
    double Y; // Y=n/s
    double Tfo,Tmax; // Freeze out and maximal temperature
    
    int full_comput; // Switch to deactivate the fast freeze out temperature determination
  
    double table_eff[276][3];   // Reads values from the SgStar files
   
    int use_table_rhoPD;
    double table_rhoPD[2][NTABMAX];
    int size_table_rhoPD;

    /*---------------------*/
    /* AlterBBN parameters */
    /*---------------------*/
    
    int err;
    int failsafe;
    double eta0;                // Initial Baryon to photon ratio
    double Nnu;                 // Number of Neutrinos (e+- included)
    double dNnu;                // Number of extra neutrinos (delta N_nu)
    double life_neutron,life_neutron_error;     // neutron lifetime
    double xinu1,xinu2,xinu3;   // [e-,neutrino], [muon,neutrino], [tau,neutrino] respectively (degeneracy parameters)
    double m_chi;               // Mass of WIMP
    double g_chi;
    double Tinit;
    int wimp;                   // Switch to enable (1) / disable (0) wimps
    int SMC_wimp;               // wimp coupling to SM particles. 1 for EM, 2 for neutrino, 3 for neut. and eq. neut.
    int selfConjugate;          // 1/0 for self-conjugate/non-self-conjugate WIMP
    int fermion;
    int EM_coupled, neut_coupled, neuteq_coupled;
    double chi2;
    int nobs;    
  };


struct errorparam
/* structure containing the cosmological model parameters */
 {
   int failsafe;                // failsafe mode
    int errnumber;              // process number for error calculation
    double random[NNUCREAC+2];  // random numbers for Monte Carlo
    double life_neutron;
 };
 
}
#endif /* defined __AlterBBN_types_hpp__ */
