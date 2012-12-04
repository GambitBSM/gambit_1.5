//
//  FortranInterface.h
//  GlobalFits
//
//  Created by Hugh Dickinson on 12/3/12.
//  Copyright (c) 2012 Hugh Dickinson. All rights reserved.
//

#ifndef GlobalFits_FortranInterface_h
#define GlobalFits_FortranInterface_h

#include "FortranFunctor.h"
#include <string>

/*********************************************************************************/

/** Class which provided a simpler interface to use the FortranFunctor. One can always revert to
 * using the bare FortranFunctor if neccessary.
 */
class FortranInterface {
	
	/// Handle for the FORTRAN library (compiled from FORTRAN code)
	void * fFortranLibraryHandle;
	
	/// Check dtatus of dynamic library interface
	void CheckDLStatus() const;
	
	public :
	
	/// Constructor.
	FortranInterface(std::string const & libFileName);
	
	/// Destructor.
	~FortranInterface();
	
	/// Set a user variable
	template<typename TagT>
    void SetVariable(typename VarType<TagT>::type const & value){
		std::string symbolName = FortranFunctor<TagT>::FortranSymbolName(VarName<TagT>::name(), VarName<TagT>::module());
		typename VarType<TagT>::type * varPtr = static_cast<typename VarType<TagT>::type *>(dlsym(fFortranLibraryHandle, symbolName.c_str()));
		CheckDLStatus();
		*varPtr = value;
	}
	
	/// Get a user variable
	template<typename TagT>
    typename VarType<TagT>::type const & GetVariable() const {
		std::string symbolName = FortranFunctor<TagT>::FortranSymbolName(VarName<TagT>::name(), VarName<TagT>::module());
		typename VarType<TagT>::type * varPtr = static_cast<typename VarType<TagT>::type *>(dlsym(fFortranLibraryHandle, symbolName.c_str()));
		CheckDLStatus();
		return *varPtr;
	}
	
	/** Load a function pointer from the ELMAG library, wrap it in a FortranFunctor, and return
	 * that functor.
	 */
	template <typename TagT>
    FortranFunctor<TagT> & Function(){
		return FortranFunctor<TagT>::instance(fFortranLibraryHandle);
	}
	
};


#endif
