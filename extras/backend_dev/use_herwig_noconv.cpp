/*
 * Test Herwig backend _without_ convenience functions
 * 
 * \author Anders Kvellestad
 * \date 2013-04-07
 * 
 */

#include <iostream>
#include "backend_rollcall.hpp"

int main()
{
  std::cout << std::endl << std::endl;
  std::cout << "BACKEND: Herwig" << std::endl;
  std::cout << "---------------" << std::endl;
  std::cout << std::endl;
  
  
  /* ==================================
   *  General setup and initialization
   * ================================== */

  std::string outputPrefix = "Herwig_backend: ";

  /* Choose process: */
  std::cout << std::endl << outputPrefix << "Set process specifications" << std::endl;

  hwproc_type tempHwproc = GAMBIT::Backend::Herwig::getHwproc();
  tempHwproc.pbeam1 = 3500.;
  tempHwproc.pbeam2 = 3500.;
  tempHwproc.iproc  = 1500;
  tempHwproc.maxev  = 999999999;	// Should be ridiculously large number so that we never get to it
  GAMBIT::Backend::Herwig::setHwproc(tempHwproc);

  hwbmch_type tempHwbmch = GAMBIT::Backend::Herwig::getHwbmch();
  strcpy(tempHwbmch.part1, "P       "); // Must be 8 characters long
  strcpy(tempHwbmch.part2, "P       ");
  GAMBIT::Backend::Herwig::setHwbmch(tempHwbmch);

  /* Initialize common blocks */ 
  std::cout << std::endl << outputPrefix << "Initialize common blocks" << std::endl;
  GAMBIT::Backend::Herwig::hwigin();

  /* Overwrite any default HERWIG common block values */
  hwevnt_type myHwevnt = GAMBIT::Backend::Herwig::getHwevnt(); 
  myHwevnt.maxpr = 0; // Number of events printed to screen
  GAMBIT::Backend::Herwig::setHwevnt(myHwevnt); 


  /* =========================================
   *  Initialization for each parameter point
   * ========================================= */
  
  /* Read parameter point provided from Core/ScannerBit (not implemented) */
  std::cout << std::endl << outputPrefix << "Read parameter point (not implemented yet)" << std::endl;

  /* Compute parameter dependent quantities */
  std::cout << std::endl << outputPrefix << "Compute parameter dependent constants" << std::endl;
  GAMBIT::Backend::Herwig::hwuinc();
  
  /* Set seed */
  int seed1 = 1234;
  int seed2 = 5678;
  std::cout << std::endl << outputPrefix << "Set seeds to: " << seed1 << ", " << seed2 << std::endl;
  hwevnt_type tempHwevnt = GAMBIT::Backend::Herwig::getHwevnt();
  tempHwevnt.nrn[0] = seed1;
  tempHwevnt.nrn[1] = seed2;
  GAMBIT::Backend::Herwig::setHwevnt(tempHwevnt);
  
  /* Initialize process */
  std::cout << std::endl << outputPrefix << "Initialize elementary process"<< std::endl;
  GAMBIT::Backend::Herwig::hweini();
  
  
  /* ==================
   *  Event generation
   * ================== */

  bool hadronize = false;
  bool soft = false;

  /* Generate a number of events */
  for(int i = 0; i<10 ;i++)
  {
    std::cout << std::endl << outputPrefix << "Generate event"<< std::endl;

    /* Initialize event */
    GAMBIT::Backend::Herwig::hwuine();
    /* Generate hard process */
    GAMBIT::Backend::Herwig::hwepro();
    /* Parton shower */
    GAMBIT::Backend::Herwig::hwbgen();
    /* Heavy particle decay */
    GAMBIT::Backend::Herwig::hwdhob();

    /* Hadronization: */
    if (hadronize)
    {
      /* Cluster formation */
      GAMBIT::Backend::Herwig::hwcfor();
      /* Cluster decay */
      GAMBIT::Backend::Herwig::hwcdec();
      /* Unstable particle decay */
      GAMBIT::Backend::Herwig::hwdhad();
      /* Heavy flavour decay */
      GAMBIT::Backend::Herwig::hwdhvy();
      /* Soft underlying event */
      if (soft) { GAMBIT::Backend::Herwig::hwmevt(); }
    }

    /* Finishing touches to event */
    GAMBIT::Backend::Herwig::hwufne();

    /* Read event data */
    std::cout << std::endl << outputPrefix << "Read event data and return it in a suitable format (not implemented yet)" << std::endl;

  } // End event loop


  /* ==========
   *  Finalize
   * ========== */
 
  /* Finalize elementary process */
  std::cout << std::endl << outputPrefix << "Finalize elementary process" << std::endl;
  GAMBIT::Backend::Herwig::hwefin();
 
 

  std::cout << std::endl;
  return 0;
}
