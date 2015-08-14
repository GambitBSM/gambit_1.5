//==========================================================================
// This file has been automatically generated for Pythia 8 by
// MadGraph5_aMC@NLO v. 2.3.0, 2015-07-01
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "Sigma_MC4BSM_2012_UFO_qq_evevx.h"
#include "HelAmps_MC4BSM_2012_UFO.h"

using namespace Pythia8_MC4BSM_2012_UFO; 

namespace Pythia8 
{

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: u u~ > ev ev~ WEIGHTED=4
// Process: c c~ > ev ev~ WEIGHTED=4
// Process: d d~ > ev ev~ WEIGHTED=4
// Process: s s~ > ev ev~ WEIGHTED=4

//--------------------------------------------------------------------------
// Initialize process.

void Sigma_MC4BSM_2012_UFO_qq_evevx::initProc() 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_MC4BSM_2012_UFO::getInstance(); 
  pars->setIndependentParameters(particleDataPtr, couplingsPtr, slhaPtr); 
  pars->setIndependentCouplings(); 
  // Set massive/massless matrix elements for c/b/mu/tau
  mcME = 0.; 
  mbME = particleDataPtr->m0(5); 
  mmuME = 0.; 
  mtauME = particleDataPtr->m0(15); 
  jamp2[0] = new double[1]; 
  jamp2[1] = new double[1]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void Sigma_MC4BSM_2012_UFO_qq_evevx::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(particleDataPtr, couplingsPtr, slhaPtr, alpS); 
  pars->setDependentCouplings(); 
  // Reset color flows
  for(int i = 0; i < 1; i++ )
    jamp2[0][i] = 0.; 
  for(int i = 0; i < 1; i++ )
    jamp2[1][i] = 0.; 

  // Local variables and constants
  const int ncomb = 16; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1}, {-1, -1,
      -1, 1}, {-1, -1, 1, -1}, {-1, -1, 1, 1}, {-1, 1, -1, -1}, {-1, 1, -1, 1},
      {-1, 1, 1, -1}, {-1, 1, 1, 1}, {1, -1, -1, -1}, {1, -1, -1, 1}, {1, -1,
      1, -1}, {1, -1, 1, 1}, {1, 1, -1, -1}, {1, 1, -1, 1}, {1, 1, 1, -1}, {1,
      1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {36, 36, 36, 36}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
    t[i] = 0.; 
  }

  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  // For now, call setupForME() here
  id1 = 2; 
  id2 = -2; 
  if( !setupForME())
  {
    return; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_uux_evevx(); 
        t[1] = matrix_ddx_evevx(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[2] = matrix_uux_evevx(); 
        t[3] = matrix_ddx_evevx(); 
        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc]; 
          tsum += t[iproc]; 
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true; 
          ngood++; 
          igood[ngood] = ihel; 
        }
      }
    }
    jhel = 0; 
    sum_hel = min(sum_hel, ngood); 
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_uux_evevx(); 
      t[1] = matrix_ddx_evevx(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[2] = matrix_uux_evevx(); 
      t[3] = matrix_ddx_evevx(); 
      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double Sigma_MC4BSM_2012_UFO_qq_evevx::sigmaHat() 
{
  // Select between the different processes
  if(id1 == -3 && id2 == 3)
  {
    // Add matrix elements for processes with beams (-3, 3)
    return matrix_element[3]; 
  }
  else if(id1 == -1 && id2 == 1)
  {
    // Add matrix elements for processes with beams (-1, 1)
    return matrix_element[3]; 
  }
  else if(id1 == 4 && id2 == -4)
  {
    // Add matrix elements for processes with beams (4, -4)
    return matrix_element[0]; 
  }
  else if(id1 == 2 && id2 == -2)
  {
    // Add matrix elements for processes with beams (2, -2)
    return matrix_element[0]; 
  }
  else if(id1 == -4 && id2 == 4)
  {
    // Add matrix elements for processes with beams (-4, 4)
    return matrix_element[2]; 
  }
  else if(id1 == -2 && id2 == 2)
  {
    // Add matrix elements for processes with beams (-2, 2)
    return matrix_element[2]; 
  }
  else if(id1 == 3 && id2 == -3)
  {
    // Add matrix elements for processes with beams (3, -3)
    return matrix_element[1]; 
  }
  else if(id1 == 1 && id2 == -1)
  {
    // Add matrix elements for processes with beams (1, -1)
    return matrix_element[1]; 
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.; 
  }
}

//--------------------------------------------------------------------------
// Select identity, colour and anticolour.

void Sigma_MC4BSM_2012_UFO_qq_evevx::setIdColAcol() 
{
  if(id1 == -3 && id2 == 3)
  {
    // Pick one of the flavor combinations
    int flavors[1][2] = {{9000009, -9000009}}; 
    vector<double> probs; 
    double sum = matrix_element[3]; 
    probs.push_back(matrix_element[3]/sum); 
    int choice = rndmPtr->pick(probs); 
    id3 = flavors[choice][0]; 
    id4 = flavors[choice][1]; 
  }
  else if(id1 == -1 && id2 == 1)
  {
    // Pick one of the flavor combinations
    int flavors[1][2] = {{9000009, -9000009}}; 
    vector<double> probs; 
    double sum = matrix_element[3]; 
    probs.push_back(matrix_element[3]/sum); 
    int choice = rndmPtr->pick(probs); 
    id3 = flavors[choice][0]; 
    id4 = flavors[choice][1]; 
  }
  else if(id1 == 4 && id2 == -4)
  {
    // Pick one of the flavor combinations (9000009, -9000009)
    int flavors[1][2] = {{9000009, -9000009}}; 
    vector<double> probs; 
    double sum = matrix_element[0]; 
    probs.push_back(matrix_element[0]/sum); 
    int choice = rndmPtr->pick(probs); 
    id3 = flavors[choice][0]; 
    id4 = flavors[choice][1]; 
  }
  else if(id1 == 2 && id2 == -2)
  {
    // Pick one of the flavor combinations (9000009, -9000009)
    int flavors[1][2] = {{9000009, -9000009}}; 
    vector<double> probs; 
    double sum = matrix_element[0]; 
    probs.push_back(matrix_element[0]/sum); 
    int choice = rndmPtr->pick(probs); 
    id3 = flavors[choice][0]; 
    id4 = flavors[choice][1]; 
  }
  else if(id1 == -4 && id2 == 4)
  {
    // Pick one of the flavor combinations
    int flavors[1][2] = {{9000009, -9000009}}; 
    vector<double> probs; 
    double sum = matrix_element[2]; 
    probs.push_back(matrix_element[2]/sum); 
    int choice = rndmPtr->pick(probs); 
    id3 = flavors[choice][0]; 
    id4 = flavors[choice][1]; 
  }
  else if(id1 == -2 && id2 == 2)
  {
    // Pick one of the flavor combinations
    int flavors[1][2] = {{9000009, -9000009}}; 
    vector<double> probs; 
    double sum = matrix_element[2]; 
    probs.push_back(matrix_element[2]/sum); 
    int choice = rndmPtr->pick(probs); 
    id3 = flavors[choice][0]; 
    id4 = flavors[choice][1]; 
  }
  else if(id1 == 3 && id2 == -3)
  {
    // Pick one of the flavor combinations (9000009, -9000009)
    int flavors[1][2] = {{9000009, -9000009}}; 
    vector<double> probs; 
    double sum = matrix_element[1]; 
    probs.push_back(matrix_element[1]/sum); 
    int choice = rndmPtr->pick(probs); 
    id3 = flavors[choice][0]; 
    id4 = flavors[choice][1]; 
  }
  else if(id1 == 1 && id2 == -1)
  {
    // Pick one of the flavor combinations (9000009, -9000009)
    int flavors[1][2] = {{9000009, -9000009}}; 
    vector<double> probs; 
    double sum = matrix_element[1]; 
    probs.push_back(matrix_element[1]/sum); 
    int choice = rndmPtr->pick(probs); 
    id3 = flavors[choice][0]; 
    id4 = flavors[choice][1]; 
  }
  setId(id1, id2, id3, id4); 
  // Pick color flow
  int ncolor[2] = {1, 1}; 
  if((id1 == 2 && id2 == -2 && id3 == 9000009 && id4 == -9000009) || (id1 == 4
      && id2 == -4 && id3 == 9000009 && id4 == -9000009))
  {
    vector<double> probs; 
    double sum = jamp2[0][0]; 
    for(int i = 0; i < ncolor[0]; i++ )
      probs.push_back(jamp2[0][i]/sum); 
    int ic = rndmPtr->pick(probs); 
    static int colors[1][8] = {{1, 0, 0, 1, 0, 0, 0, 0}}; 
    setColAcol(colors[ic][0], colors[ic][1], colors[ic][2], colors[ic][3],
        colors[ic][4], colors[ic][5], colors[ic][6], colors[ic][7]);
  }
  else if((id1 == 1 && id2 == -1 && id3 == 9000009 && id4 == -9000009) || (id1
      == 3 && id2 == -3 && id3 == 9000009 && id4 == -9000009))
  {
    vector<double> probs; 
    double sum = jamp2[1][0]; 
    for(int i = 0; i < ncolor[1]; i++ )
      probs.push_back(jamp2[1][i]/sum); 
    int ic = rndmPtr->pick(probs); 
    static int colors[1][8] = {{1, 0, 0, 1, 0, 0, 0, 0}}; 
    setColAcol(colors[ic][0], colors[ic][1], colors[ic][2], colors[ic][3],
        colors[ic][4], colors[ic][5], colors[ic][6], colors[ic][7]);
  }
  else if((id1 == -2 && id2 == 2 && id3 == 9000009 && id4 == -9000009) || (id1
      == -4 && id2 == 4 && id3 == 9000009 && id4 == -9000009))
  {
    vector<double> probs; 
    double sum = jamp2[0][0]; 
    for(int i = 0; i < ncolor[0]; i++ )
      probs.push_back(jamp2[0][i]/sum); 
    int ic = rndmPtr->pick(probs); 
    static int colors[1][8] = {{0, 1, 1, 0, 0, 0, 0, 0}}; 
    setColAcol(colors[ic][0], colors[ic][1], colors[ic][2], colors[ic][3],
        colors[ic][4], colors[ic][5], colors[ic][6], colors[ic][7]);
  }
  else if((id1 == -1 && id2 == 1 && id3 == 9000009 && id4 == -9000009) || (id1
      == -3 && id2 == 3 && id3 == 9000009 && id4 == -9000009))
  {
    vector<double> probs; 
    double sum = jamp2[1][0]; 
    for(int i = 0; i < ncolor[1]; i++ )
      probs.push_back(jamp2[1][i]/sum); 
    int ic = rndmPtr->pick(probs); 
    static int colors[1][8] = {{0, 1, 1, 0, 0, 0, 0, 0}}; 
    setColAcol(colors[ic][0], colors[ic][1], colors[ic][2], colors[ic][3],
        colors[ic][4], colors[ic][5], colors[ic][6], colors[ic][7]);
  }
}

//--------------------------------------------------------------------------
// Evaluate weight for angles of decay products in process

double Sigma_MC4BSM_2012_UFO_qq_evevx::weightDecay(Event& process, int iResBeg,
    int iResEnd)
{
  // Just use isotropic decay (default)
  return 1.; 
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void Sigma_MC4BSM_2012_UFO_qq_evevx::calculate_wavefunctions(const int perm[],
    const int hel[])
{
  // Calculate wavefunctions for all processes
  double p[nexternal][4]; 
  int i; 

  // Convert Pythia 4-vectors to double[]
  for(i = 0; i < nexternal; i++ )
  {
    p[i][0] = pME[i].e(); 
    p[i][1] = pME[i].px(); 
    p[i][2] = pME[i].py(); 
    p[i][3] = pME[i].pz(); 
  }

  // Calculate all wavefunctions
  ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  FFV1P0_3(w[0], w[1], pars->GC_2, pars->ZERO, pars->ZERO, w[4]); 
  FFV2_5_3(w[0], w[1], pars->GC_28, pars->GC_40, pars->mdl_MZ, pars->mdl_WZ,
      w[5]);
  FFV1P0_3(w[0], w[1], pars->GC_1, pars->ZERO, pars->ZERO, w[6]); 
  FFV2_3_3(w[0], w[1], pars->GC_27, pars->GC_40, pars->mdl_MZ, pars->mdl_WZ,
      w[7]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[3], w[2], w[4], pars->GC_13, amp[0]); 
  FFV1_0(w[3], w[2], w[5], pars->GC_43, amp[1]); 
  FFV1_0(w[3], w[2], w[6], pars->GC_13, amp[2]); 
  FFV1_0(w[3], w[2], w[7], pars->GC_43, amp[3]); 


}
double Sigma_MC4BSM_2012_UFO_qq_evevx::matrix_uux_evevx() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 2; 
  const int ncolor = 1; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1}; 
  static const double cf[ncolor][ncolor] = {{3}}; 

  // Calculate color flows
  jamp[0] = -amp[0] - amp[1]; 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}

double Sigma_MC4BSM_2012_UFO_qq_evevx::matrix_ddx_evevx() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 2; 
  const int ncolor = 1; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1}; 
  static const double cf[ncolor][ncolor] = {{3}}; 

  // Calculate color flows
  jamp[0] = -amp[2] - amp[3]; 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[1][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}


}  // end namespace Pythia8
