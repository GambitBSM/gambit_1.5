/*
 * Sketch of backend to HERWIG 6.51 (Fortran).
 * 
 * \author Anders Kvellestad
 * \date 2013-04-03
 * 
 * Modified: 2013-04-07
 */

#include <string.h>
#include "backend_general.hpp"


/* Start stuff inside include braces; this may need to go in another separate header. */

#ifndef __BACKEND_LIBFORTRANCODE_HPP__
#define __BACKEND_LIBFORTRANCODE_HPP__

  /* We must give the correct size for a bunch of HERWIG arrays in order to
   * successfully manipulate the common blocks. This size must be the
   * same as the parameter 'NMXHEP' in the Herwig source used to compile the
   * shared library. 
   * (In HERWIG 6.51 this parameter is defined in the file 'herwig6510.inc'.) */

  static const int NMXHEP = 4000; 


  /* Structs to be used as types for specific common blocks
   * in the HERWIG library. */

  struct hwproc_type
  {
    double ebeam1, ebeam2, pbeam1, pbeam2;
    int iproc, maxev;
  };

  struct hwbmch_type
  {
    char part1[8], part2[8];
  };

  struct hepevt_type
  {
    int nevhep, nhep, isthep[NMXHEP], idhep[NMXHEP], jmohep[NMXHEP][2], jdahep[NMXHEP][2];
    double phep[NMXHEP][5], vhep[NMXHEP][4];
  };

  struct hwevnt_type
  {
    double avwgt, evwgt, gamwt, tlout, wbigst, wgtmax, wgtsum, wsqsum;
    int idhw[NMXHEP], ierror, istat, lwevt, maxer, maxpr, nowgt, nrn[2], numer, numeru, nwgts, gensof;
  };

#endif // end stuff in include braces


/* Specify the library path along with a backend name. */

#define LIBPATH "./libherwig.so"
#ifdef BACKENDRENAME
  #define BACKENDNAME BACKENDRENAME
#else
  #define BACKENDNAME Herwig
#endif


/* Load library, obtain pointers to the various library symbols
 * and set up a minimal interface consisting of get/set functions
 * for the variables and function pointers for the functions. */

LOAD_LIBRARY

BE_VARIABLE(Hwproc, hwproc_type, "hwproc_", pHwproc)
BE_VARIABLE(Hwbmch, hwbmch_type, "hwbmch_", pHwbmch)
BE_VARIABLE(Hepevt, hepevt_type, "hepevt_", pHepevt)
BE_VARIABLE(Hwevnt, hwevnt_type, "hwevnt_", pHwevnt)

BE_FUNCTION(hwigin, void, (), "hwigin_")
BE_FUNCTION(inimas, void, (), "inimas_")
BE_FUNCTION(hwuinc, void, (), "hwuinc_")
//BE_FUNCTION(hwabeg, void, (), "hwabeg_")
BE_FUNCTION(hweini, void, (), "hweini_")
BE_FUNCTION(hwuine, void, (), "hwuine_")
BE_FUNCTION(hwefin, void, (), "hwefin_")
BE_FUNCTION(hwepro, void, (), "hwepro_")
BE_FUNCTION(hwbgen, void, (), "hwbgen_")
BE_FUNCTION(hwdhob, void, (), "hwdhob_")
BE_FUNCTION(hwcfor, void, (), "hwcfor_")
BE_FUNCTION(hwcdec, void, (), "hwcdec_")
BE_FUNCTION(hwdhad, void, (), "hwdhad_")
BE_FUNCTION(hwdhvy, void, (), "hwdhvy_")
BE_FUNCTION(hwmevt, void, (), "hwmevt_")
BE_FUNCTION(hwufne, void, (), "hwufne_")


/* End of minimal backend setup. 
 * Any additional convenience functions can be defined below using
 * the variables and function pointers. */

namespace GAMBIT                                                           
{                                                                            
  namespace Backend                                                        
  {                                                                          
    namespace BACKENDNAME                                                  
    {                                                                        

      /*
       * Convenience functions: 
       */

      /* A prefix used for all cout statements */
      std::string outputPrefix = std::string(STRINGIFY(BACKENDNAME)) + std::string("_backend: ");


      /* function: chooseProcess */
      void chooseProcess(double pbeam1, double pbeam2, int iproc, const char part1[8], const char part2[8], int maxev = 999999999)
      {
        std::cout << std::endl << outputPrefix << "Set process specifications" << std::endl;
        hwproc_type tempHwproc = getHwproc();
        tempHwproc.pbeam1 = pbeam1;
        tempHwproc.pbeam2 = pbeam2;
        tempHwproc.iproc  = iproc;
        tempHwproc.maxev  = maxev;	// Should be ridiculously large number so that we never get to it
        setHwproc(tempHwproc);

        hwbmch_type tempHwbmch = getHwbmch();
        strcpy(tempHwbmch.part1, part1); // Must be 8 characters long
        strcpy(tempHwbmch.part2, part2);
        setHwbmch(tempHwbmch);
      }


      /* function: initCommonBlocks */
      void initCommonBlocks()
      {
        /* Set common block values to HERWIG default values. */ 
        std::cout << std::endl << outputPrefix << "Initialize common blocks" << std::endl;
        hwigin();
      }


      /* function: readParamPoint */
      void readParamPoint()
      {
        /* Read parameter point */
        std::cout << std::endl << outputPrefix << "Read parameter point (not implemented yet)" << std::endl;
        
        /* Compute parameter dependent quantities */
        std::cout << std::endl << outputPrefix << "Compute parameter dependent constants" << std::endl;
        hwuinc();
      }


      /* function: initElemProcess */
      void initElemProcess(int seed1, int seed2)
      {
        /* Set seed */
        std::cout << std::endl << outputPrefix << "Set seeds to: " << seed1 << ", " << seed2 << std::endl;
        hwevnt_type tempHwevnt = getHwevnt();
        tempHwevnt.nrn[0] = seed1;
        tempHwevnt.nrn[1] = seed2;
        setHwevnt(tempHwevnt);
        
        /* Initialize process */
        std::cout << std::endl << outputPrefix << "Initialize elementary process"<< std::endl;
        hweini();
      }


      /* function: generateEvent */
      void generateEvent(bool hadronize = true, bool soft = true)
      {
        std::cout << std::endl << outputPrefix << "Generate event"<< std::endl;

        /* Initialize event */
        hwuine();
        /* Generate hard process */
        hwepro();
        /* Parton shower */
        hwbgen();
        /* Heavy particle decay */
        hwdhob();

        /* Hadronization: */
        if (hadronize)
        {
          /* Cluster formation */
          hwcfor();
          /* Cluster decay */
          hwcdec();
          /* Unstable particle decay */
          hwdhad();
          /* Heavy flavour decay */
          hwdhvy();
          /* Soft underlying event */
          if (soft) { hwmevt(); }
        }

        /* Finishing touches to event */
        hwufne();
      }


      /* function: readEvent */
      void readEvent()
      {
        std::cout << std::endl << outputPrefix << "Read event data and return it in a suitable format (not implemented yet)" << std::endl;
      }
      

      /* function: finalizeElemProcess */
      void finalizeElemProcess()
      {
        std::cout << std::endl << outputPrefix << "Finalize elementary process" << std::endl;
        hwefin();
      }


    } /* end namespace BACKENDNAME */                                          
  } /* end namespace Backend */                                                
} /* end namespace GAMBIT */                                                   


// Undefine macros to avoid conflict with other backends
#undef LIBPATH 
#undef BACKENDNAME


