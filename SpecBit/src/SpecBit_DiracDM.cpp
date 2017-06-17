//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions of module SpecBit
///
///  SpecBit module functions related to the
///  DiracDM model
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ankit Beniwal
///          (ankit.beniwal@adelaide.edu.au)
///  \date 2016 Oct, Nov
///  \date 2017 Jun
///
///  *********************************************

#include <string>
#include <sstream>

#include "gambit/Elements/gambit_module_headers.hpp"

#include "gambit/Elements/spectrum.hpp"
#include "gambit/Utils/stream_overloads.hpp" 
#include "gambit/Utils/util_macros.hpp"

#include "gambit/SpecBit/SpecBit_rollcall.hpp"
#include "gambit/SpecBit/SpecBit_helpers.hpp"
#include "gambit/SpecBit/QedQcdWrapper.hpp"
#include "gambit/Models/SimpleSpectra/SMHiggsSimpleSpec.hpp"
#include "gambit/Models/SimpleSpectra/DiracDMSimpleSpec.hpp"
#include "gambit/SpecBit/model_files_and_boxes.hpp"

#include "gambit/SpecBit/DiracDMSpec.hpp"

// Switch for debug mode
//#define SPECBIT_DEBUG

namespace Gambit
{

  namespace SpecBit
  {
    using namespace LogTags;
    using namespace flexiblesusy;

    /// Get a (simple) Spectrum object wrapper for the DiracDM model
    void get_DiracDM_spectrum(Spectrum& result)
    {
      namespace myPipe = Pipes::get_DiracDM_spectrum;
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Initialise an object to carry the Dirac plus Higgs sector information
      Models::DiracDMModel diracmodel;
      
      // quantities needed to fill container spectrum, intermediate calculations
      double alpha_em = 1.0 / sminputs.alphainv;
      double C = alpha_em * Pi / (sminputs.GF * pow(2,0.5));
      double sinW2 = 0.5 - pow( 0.25 - C/pow(sminputs.mZ,2) , 0.5);
      double cosW2 = 0.5 + pow( 0.25 - C/pow(sminputs.mZ,2) , 0.5);
      double e = pow( 4*Pi*( alpha_em ),0.5) ;
      
      // Higgs sector
      double mh   = *myPipe::Param.at("mH");
      diracmodel.HiggsPoleMass   = mh;
      
      double vev        = 1. / sqrt(sqrt(2.)*sminputs.GF);
      diracmodel.HiggsVEV        = vev;
      // diracmodel.LambdaH   = GF*pow(mh,2)/pow(2,0.5) ;
      
      // DiracDM sector    
      diracmodel.DiracPoleMass = *myPipe::Param.at("mF");
      diracmodel.DiracLambda   = *myPipe::Param.at("lF");    
      diracmodel.DiraccosXI    = *myPipe::Param.at("cosXI");    
      
       // Check if the EFT is valid (i.e., lF < 4*pi/2*mF)
      if (diracmodel.DiracLambda >= (4*Pi)/(2*diracmodel.DiracPoleMass))
      {
       std::ostringstream msg;
       msg << "Parameter point [mF, lF] = [" << diracmodel.DiracPoleMass
           << " GeV, " << diracmodel.DiracLambda
           << " GeV] is invalid under the EFT validity constraint!";
       invalid_point().raise(msg.str());
      }

      // Standard model
      diracmodel.sinW2 = sinW2;
      
      // gauge couplings
      diracmodel.g1 = e / sinW2;
      diracmodel.g2 = e / cosW2;
      diracmodel.g3   = pow( 4*Pi*( sminputs.alphaS ),0.5) ;
      
      // Yukawas
      double sqrt2v = pow(2.0,0.5)/vev;     
      diracmodel.Yu[0] = sqrt2v * sminputs.mU;
      diracmodel.Yu[1] = sqrt2v * sminputs.mCmC;
      diracmodel.Yu[2] = sqrt2v * sminputs.mT;      
      diracmodel.Ye[0] = sqrt2v * sminputs.mE;
      diracmodel.Ye[1] = sqrt2v * sminputs.mMu;
      diracmodel.Ye[2] = sqrt2v * sminputs.mTau;   
      diracmodel.Yd[0] = sqrt2v * sminputs.mD;
      diracmodel.Yd[1] = sqrt2v * sminputs.mS;
      diracmodel.Yd[2] = sqrt2v * sminputs.mBmB;
            
      // Create a SubSpectrum object to wrap the EW sector information
      Models::DiracDMSimpleSpec diracspec(diracmodel);

      // Retrieve any mass cuts
      static const Spectrum::mc_info mass_cut = myPipe::runOptions->getValueOrDef<Spectrum::mc_info>(Spectrum::mc_info(), "mass_cut");
      static const Spectrum::mr_info mass_ratio_cut = myPipe::runOptions->getValueOrDef<Spectrum::mr_info>(Spectrum::mr_info(), "mass_ratio_cut");

      // We don't supply a LE subspectrum here; an SMSimpleSpec will therefore be automatically created from 'sminputs'
      result = Spectrum(diracspec,sminputs,&myPipe::Param,mass_cut,mass_ratio_cut);

    }    
    
    // print spectrum out, stripped down copy from MSSM version with variable names changed
    void fill_map_from_DiracDMspectrum(std::map<std::string,double>&, const Spectrum&);
   
    void get_DiracDM_spectrum_as_map (std::map<std::string,double>& specmap)
    {
      namespace myPipe = Pipes::get_DiracDM_spectrum_as_map;
      const Spectrum& diracdmspec(*myPipe::Dep::DiracDM_spectrum);
      fill_map_from_DiracDMspectrum(specmap, diracdmspec);
    }
    
    void fill_map_from_DiracDMspectrum(std::map<std::string,double>& specmap, const Spectrum& diracdmspec)
    {
      /// Add everything... use spectrum contents routines to automate task
      static const SpectrumContents::DiracDM contents;
      static const std::vector<SpectrumParameter> required_parameters = contents.all_parameters();
      
      for(std::vector<SpectrumParameter>::const_iterator it = required_parameters.begin();
           it != required_parameters.end(); ++it)
      {
         const Par::Tags        tag   = it->tag();
         const std::string      name  = it->name();
         const std::vector<int> shape = it->shape();

         /// Verification routine should have taken care of invalid shapes etc, so won't check for that here.

         // Check scalar case
         if(shape.size()==1 and shape[0]==1)
         {
           std::ostringstream label;
           label << name <<" "<< Par::toString.at(tag);
           specmap[label.str()] = diracdmspec.get_HE().get(tag,name);
         }
         // Check vector case
         else if(shape.size()==1 and shape[0]>1)
         {
           for(int i = 1; i<=shape[0]; ++i) {
             std::ostringstream label;
             label << name <<"_"<<i<<" "<< Par::toString.at(tag);
             specmap[label.str()] = diracdmspec.get_HE().get(tag,name,i);
           }
         }
         // Check matrix case
         else if(shape.size()==2)
         {
           for(int i = 1; i<=shape[0]; ++i) {
             for(int j = 1; j<=shape[0]; ++j) {
               std::ostringstream label;
               label << name <<"_("<<i<<","<<j<<") "<<Par::toString.at(tag);
               specmap[label.str()] = diracdmspec.get_HE().get(tag,name,i,j);
             }  
           }
         }
         // Deal with all other cases
         else
         {
           // ERROR
           std::ostringstream errmsg;           
           errmsg << "Error, invalid parameter received while converting DiracDMspectrum to map of strings! This should no be possible if the spectrum content verification routines were working correctly; they must be buggy, please report this.";
           errmsg << "Problematic parameter was: "<< tag <<", " << name << ", shape="<< shape; 
           utils_error().forced_throw(LOCAL_INFO,errmsg.str());
         }
      }

    }
    
  } // end namespace SpecBit
} // end namespace Gambit

