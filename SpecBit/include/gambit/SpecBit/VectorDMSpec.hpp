//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  VectorDM derived version of SubSpectrum class.
///
///  *********************************************
///
///  Authors: 
///  <!-- add name and date if you modify -->
///   
///  \author Ankit Beniwal  
///          (ankit.beniwal@adelaide.edu.au)
///  \date 2016 Sep, 2017 Jun
///
///  *********************************************

#ifndef VECTORDMSPEC_H
#define VECTORDMSPEC_H

#include <memory>

#include "gambit/cmake/cmake_variables.hpp"
#include "gambit/Elements/subspectrum.hpp"
#include "gambit/Elements/slhaea_helpers.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "gambit/SpecBit/VectorDMSpec_head.hpp"   // "Header" declarations for VectorDMSpec class

// Flexible SUSY stuff (should not be needed by the rest of gambit)
#include "flexiblesusy/config/config.h"


namespace Gambit
{

   namespace SpecBit
   {
      template <class MI>
      const int VectorDMSpec<MI>::_index_offset = MI::index_offset;

      template <class MI>
      VectorDMSpec<MI>::VectorDMSpec(MI mi, str be_name, str be_version)
         : backend_name(be_name)
         , backend_version(be_version)
         , model_interface(mi)
      {}
  
  
      template <class MI>
      VectorDMSpec<MI>::VectorDMSpec()
      {}
     
      template <class MI>
      VectorDMSpec<MI>::~VectorDMSpec()
      {}
     
     
     
      template <class MI>
      void VectorDMSpec<MI>::RunToScaleOverride(double scale)
      {
        model_interface.model.run_to(scale);
      }
      template <class MI>
      double VectorDMSpec<MI>::GetScale() const
      {
        return model_interface.model.get_scale();
      }
      template <class MI>
      void VectorDMSpec<MI>::SetScale(double scale)
      {
        model_interface.model.set_scale(scale);
      }
      
      template <class MI>
      std::string VectorDMSpec<MI>::AccessError(std::string state) const
      {
        std::string errormsg;
        errormsg = "Error accessing "+ state + " element is out of bounds";
        return errormsg;
      }
     
    
      template <class Model>
      double get_sinthW2_MSbar(const Model& model)
      {
       double sthW2 = Utils::sqr(model.get_g1()) * 0.6 / 
                      (0.6 * Utils::sqr(model.get_g1()) + 
                      Utils::sqr(model.get_g2()));
       return sthW2;
      }
      

      template <class MI>
      typename VectorDMSpec<MI>::GetterMaps VectorDMSpec<MI>::fill_getter_maps()
      {
         typename VectorDMSpec<MI>::GetterMaps map_collection;
         typedef typename MI::Model Model;

         typedef typename MTget::FInfo1 FInfo1;
         typedef typename MTget::FInfo2 FInfo2;

         // Can't use c++11 initialise lists, se have to initialise the index sets like this.
         static const int i01v[] = {0,1};
         static const std::set<int> i01(i01v, Utils::endA(i01v));

         static const int i012v[] = {0,1,2};
         static const std::set<int> i012(i012v, Utils::endA(i012v));

         static const int i0123v[] = {0,1,2,3};
         static const std::set<int> i0123(i0123v, Utils::endA(i0123v));

         static const int i012345v[] = {0,1,2,3,4,5};
         static const std::set<int> i012345(i012345v, Utils::endA(i012345v));
        

        {
            typename MTget::fmap0 tmp_map;
            tmp_map["mu2"] = &Model::get_muH;
            map_collection[Par::mass2].map0 = tmp_map;
         }
        


         //dimensionless - mass dimension 0 parameters

         {
            typename MTget::fmap0 tmp_map;
            tmp_map["g1"]= &Model::get_g1;
            tmp_map["g2"]= &Model::get_g2;
            tmp_map["g3"]= &Model::get_g3;
            tmp_map["lambda_h"]= &Model::get_LamH;
            tmp_map["lambda_hV"]= &Model::get_LamVH;

            map_collection[Par::dimensionless].map0 = tmp_map;
         }
        
        
          {
            typename MTget::fmap0_extraM tmp_map;
            tmp_map["sinW2"] = &get_sinthW2_MSbar<Model>;
            map_collection[Par::dimensionless].map0_extraM = tmp_map;
         }
        


         {
            typename MTget::fmap0 tmp_map;
            tmp_map["vev"] = &Model::get_v;
            tmp_map["lambda_hV"]= &Model::get_LamVH;  // ??? not sure why it wants it here
            map_collection[Par::mass1].map0 = tmp_map;
         }

         // Functions utilising the "extraM" function signature
         // (Zero index, model object as argument)

   
         // Functions utilising the two-index "plain-vanilla" function signature
         // (Two-index member functions of model object)
         {
            typename MTget::fmap2 tmp_map;

            tmp_map["Yd"]= FInfo2( &Model::get_Yd, i012, i012);
            tmp_map["Yu"]= FInfo2( &Model::get_Yu, i012, i012);
            tmp_map["Ye"]= FInfo2( &Model::get_Ye, i012, i012);

            map_collection[Par::dimensionless].map2 = tmp_map;
         }
 
         // Functions utilising the one-index "plain-vanilla" function signature
         // (One-index member functions of model object)
         {
            typename MTget::fmap1 tmp_map;

            //tmp_map["V"] = FInfo1( &Model::get_Mss, i012345 );
            //tmp_map["h0"] = FInfo1( &Model::get_Mhh, i01 );

            //Here we may access the goldstone boson
            //and higgs. maybe too dangerous to keep?
            tmp_map["d"] =    FInfo1( &Model::get_MFd, i012 );
            tmp_map["u"] =    FInfo1( &Model::get_MFu, i012 );
            tmp_map["e-"] =   FInfo1( &Model::get_MFe, i012 );
            tmp_map["e"] =    FInfo1( &Model::get_MFe, i012 );
            tmp_map["dbar"] = FInfo1( &Model::get_MFd, i012 );
            tmp_map["ubar"] = FInfo1( &Model::get_MFu, i012 );
            tmp_map["e+"] =   FInfo1( &Model::get_MFe, i012 );

            map_collection[Par::mass_eigenstate].map1 = tmp_map;
         }
        
        
         {  
            typename MTget::fmap0 tmp_map;
                   
            // ***REMOVED THESE! Leave them to the QedQcdWrapper.***
            // reinstating the Z and W getters as otherwise there is no
            // point in having the setters!
            tmp_map["Z0"] = &Model::get_MVZ_pole_slha;
            //// //tmp_map["g"] = &Model::get_MGluon_pole_slha;
             tmp_map["g"] = &Model::get_MVG_pole_slha;


            map_collection[Par::Pole_Mass].map0 = tmp_map;
         } 

		{
            typename MTget::fmap0 tmp_map;

            tmp_map["V"] =  &Model::get_Mss_pole_slha;
            tmp_map["Vector"] =  &Model::get_Mss_pole_slha;
            tmp_map["h0_1"] = &Model::get_Mhh_pole_slha; //added to match SM Higgs container naming
            tmp_map["A0"] = &Model::get_MAh_pole_slha;

            map_collection[Par::Pole_Mass].map0 = tmp_map;
         }


         return map_collection;
      } 

      // Filler function for setter function pointer maps extractable from "runningpars" container
      template <class MI>
      typename VectorDMSpec<MI>::SetterMaps VectorDMSpec<MI>::fill_setter_maps()
      {
         typename VectorDMSpec<MI>::SetterMaps map_collection;
         typedef typename MI::Model Model;

         typedef typename MTset::FInfo2 FInfo2;

         // Can't use c++11 initialise lists, se have to initialise the index sets like this.
         static const int i01v[] = {0,1};
         static const std::set<int> i01(i01v, Utils::endA(i01v));

         static const int i012v[] = {0,1,2};
         static const std::set<int> i012(i012v, Utils::endA(i012v));

         static const int i0123v[] = {0,1,2,3};
         static const std::set<int> i0123(i0123v, Utils::endA(i0123v));

         static const int i012345v[] = {0,1,2,3,4,5};
         static const std::set<int> i012345(i012345v, Utils::endA(i012345v));
         
         /// @{ mass2 - mass dimension 2 parameters
         //
         // Functions utilising the "plain-vanilla" function signature
         // (Zero index member functions of model object)
         {
            typename MTset::fmap0 tmp_map;
         
            tmp_map["mu2"] = &Model::set_muH;

            map_collection[Par::mass2].map0 = tmp_map;
         }
 
         //
         // Functions utilising the "plain-vanilla" function signature
         // (Zero index member functions of model object)
         {
            typename MTset::fmap0 tmp_map;
            tmp_map["vev"] = &Model::set_v;

            map_collection[Par::mass1].map0 = tmp_map;
         }

         //
         // Functions utilising the "plain-vanilla" function signature
         // (Zero index member functions of model object)
         {
            typename MTset::fmap0 tmp_map;
            tmp_map["g1"]= &Model::set_g1;
            tmp_map["g2"]= &Model::set_g2;
            tmp_map["g3"]= &Model::set_g3;
            tmp_map["lambda_h"]= &Model::set_LamH;
            tmp_map["lambda_hV"]= &Model::set_LamVH;

            map_collection[Par::dimensionless].map0 = tmp_map;
         }

         // Functions utilising the two-index "plain-vanilla" function signature
         // (Two-index member functions of model object)
         {
            typename MTset::fmap2 tmp_map;

            tmp_map["Yd"]= FInfo2( &Model::set_Yd, i012, i012);
            tmp_map["Yu"]= FInfo2( &Model::set_Yu, i012, i012);
            tmp_map["Ye"]= FInfo2( &Model::set_Ye, i012, i012);

            map_collection[Par::dimensionless].map2 = tmp_map;
         }

         return map_collection;
      } 

      /// @}
  
   } // end SpecBit namespace 
} // end Gambit namespace

#endif
