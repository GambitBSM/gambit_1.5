/*
 * Test Herwig backend _with_ convenience functions
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

  /* Choose process: 
   * 
   * Arguments:
   * - energy beam1        (double)
   * - energy beam2        (double)
   * - Herwig process id   (int)
   * - particle type beam1 (char[8])  
   * - particle type beam2 (char[8])
   */ 
  GAMBIT::Backend::Herwig::chooseProcess(3500., 3500., 1500, "P       ", "P       ");

  /* Initialize common blocks */ 
  GAMBIT::Backend::Herwig::initCommonBlocks();

  /* Overwrite any default HERWIG common block values */
  hwevnt_type myHwevnt = GAMBIT::Backend::Herwig::getHwevnt(); 
  myHwevnt.maxpr = 0; // Number of events printed to screen
  GAMBIT::Backend::Herwig::setHwevnt(myHwevnt); 


  /* =========================================
   *  Initialization for each parameter point
   * ========================================= */
  
  /* Read parameter point provided from Core/ScannerBit (not implemented) */
  GAMBIT::Backend::Herwig::readParamPoint();
  
  /* Set seeds and initialize elementary process */
  GAMBIT::Backend::Herwig::initElemProcess(1234, 5678);
  
  
  /* ==================
   *  Event generation
   * ================== */

  /* Generate a number of events */
  for(int i = 0; i<10 ;i++)
  {
    /* Generate event:
     * 
     * Arguments:
     * - hadronize  (bool)
     * - soft event (bool)
     */
    GAMBIT::Backend::Herwig::generateEvent(false, false);

    /* Read event data */
    GAMBIT::Backend::Herwig::readEvent();
  }


  /* ==========
   *  Finalize
   * ========== */
 
  /* Finalize elementary process */
  GAMBIT::Backend::Herwig::finalizeElemProcess();
 
 

  std::cout << std::endl;
  return 0;
}
