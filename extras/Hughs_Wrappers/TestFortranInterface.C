//
//  TestFortranInterface.C
//  GlobalFits
//
//  Created by Hugh Dickinson on 12/3/12.
//  Copyright (c) 2012 Hugh Dickinson. All rights reserved.
//

#include "FortranInterface.h"
#include <iostream>

// use macros to define our interface
SUBROUTINE(TestSubroutine, testsubroutine, 2)
FUNCTION(TestFunction, testfunction, double, 2)

// main program
int main(int argc, char * argv[]){
	
	try{
		
		FortranInterface fi("testFortranLib.dylib"); // change this string as required
		
		// declare some working variables to pass to the fortran functions (it is important
		// that these have the same type as the variables expected by the FORTRAN function/subroutine
		// to which they are passed.
		float x(1.0), y(2.0);
		
		// Call a subroutine (returns void)
		fi.Function<tags::TestSubroutine>()(x, y);
		
		// call a function (returns double precision in this case)
		double z = fi.Function<tags::TestFunction>()(x, y);
		
		// write result to terminal
		std::cout << "TestFunction returned " << z << std::endl;
		
	}
	catch(FortranFunctorException & e){
		std::cout << e.what() << std::endl;
	}
	
	return 0;
}
