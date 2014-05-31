//
//  FortranInterface.C
//  GlobalFits
//
//  Created by Hugh Dickinson on 12/3/12.
//  Copyright (c) 2012 Hugh Dickinson. All rights reserved.
//

#include "FortranInterface.h"

/** Constructor.
 */
FortranInterface::FortranInterface(std::string const & libFileName):
fFortranLibraryHandle(dlopen(libFileName.c_str(), RTLD_LAZY))
{
	CheckDLStatus();
}

/** Destructor
 */
FortranInterface::~FortranInterface(){
	dlclose(fFortranLibraryHandle);
}

/** Check for dynamic library errors
 */
void FortranInterface::CheckDLStatus() const {
	const char * error = dlerror();
	if(error){
		throw(FortranFunctorException(error));
	}
}
