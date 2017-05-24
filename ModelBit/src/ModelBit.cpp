//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Function definitions of ModelBit.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2017 Apr
///
///  *********************************************


#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/ModelBit/ModelBit_rollcall.hpp"
#include "gambit/Utils/util_functions.hpp"



namespace Gambit
{

  namespace ModelBit
  {
    using namespace LogTags;

    // Module functions

    // Create the model with SARAH and generate CHep
    void generate_CHep_SARAH(bool &results)
    {
      namespace myPipe = Pipes::generate_CHep_SARAH;
      
      cout << "Running SARAH, please be patient" << endl;

      MString Model = "SSDM";
 
      myPipe::BEreq::SARAH_Start(Model);

      cout << "Making CalcHep files, seat back and wait" << endl;

      myPipe::BEreq::SARAH_MakeCHep();

    }

    // Create the model with SARAH and generate SPheno
    void generate_SPheno_SARAH(bool &results)
    {
      namespace myPipe = Pipes::generate_SPheno_SARAH;

      cout << "Running SARAH, please be patient" << endl;

      MString Model = "SSDM";

      myPipe::BEreq::SARAH_Start(Model);

      cout << "Making SPheno files, this will take a while" << endl;

      myPipe::BEreq::SARAH_MakeSPheno();
    }

  }
}
