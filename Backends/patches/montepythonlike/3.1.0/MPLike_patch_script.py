r'''   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///  script to patch MontePython to work with gambit
///
///   - harvest all likelihoods included in montepython/likelihoods folder
///   - creapte list of availible likelihoods
///   - replace import statement in each __init__ file to import MontePhytonLike instead of montepython.likelihood_class
///
///   Note: for simplicity the Data object is also defined in MontePythonLike (extra file in standard MontePython)
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2019 June
///
///  *********************************************

'''
import os
import sys
import re

def replace(file, pattern, subst, append_to_beginning=""):
    # Read contents from file as a single string
    file_handle = open(file, 'r')
    file_string = file_handle.read()
    file_handle.close()

    # Use RE package to allow for replacement (also allowing for (multi-line) REGEX)
    file_string = (re.sub(pattern, subst, file_string))

    # Write contents to file; Using mode 'w' truncates the file.
    file_handle = open(file, 'w')
    file_handle.write(append_to_beginning+file_string)
    file_handle.close()

# patch MontePythonLike to work with GAMBIT in a few lines -- just classy!
if __name__ == '__main__':
    
    # create list with all likelihood names contained in montepython/likelihoods/ folder
    output = [dI for dI in os.listdir("likelihoods/") if os.path.isdir(os.path.join('likelihoods/',dI))]
    # (JR) todo: add this line to CosmoBit to cross-check requested Likelihoods with requested ones in yaml files
    
    for like in output:
        # replace importin of montepython.likelihood_class with import of MontePythonLike
        # Note that the system path has to be inserted in the begining of the file as well
        replace("likelihoods/"+like+"/__init__.py", "from montepython.likelihood_class import", 
            "from MontePythonLike import", append_to_beginning="import sys \nsys.path.append('../../')\n" )
        # also replace importing of io_mp module (only contains input/output stream so safe to use with GAMBIT)
        replace("likelihoods/"+like+"/__init__.py", "import montepython.io_mp as io_mp", "import io_mp")

# Note: I tested to init all likes with a python script like that. Most of them work -- did not test the ones where
# data have to be downloaded separately, yet (Planck, wmap, JLA...) and
# some failed -- will check on what happens at some point (TODO for (JR) )