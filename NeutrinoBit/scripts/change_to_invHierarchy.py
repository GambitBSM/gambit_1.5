#!/bin/python
import sys
import os
import shutil
import fileinput

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)


def main(argv):
    
    directory="/net/archive/groups/plgghbt/gambit/gambit_RHN/yaml/"
    invert="/net/archive/groups/plgghbt/gambit/gambit_RHN/yamlinv/"
    for filename in os.listdir(directory ):
        if filename.endswith(".yaml"):
            filenameINV=filename
            filenameINV=filenameINV.replace('.yaml', '_inv.yaml') 
            print filenameINV
            shutil.copy2(directory+filename, invert+filenameINV)


            replaceAll(invert+filenameINV, '      range: [2e-3, 3e-3]', '      range: [-3e-3, -2e-3]')
                    


if __name__=="__main__":
    sys.exit(main(sys.argv))
