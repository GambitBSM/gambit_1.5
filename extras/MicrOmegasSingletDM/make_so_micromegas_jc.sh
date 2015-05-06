#!/bin/sh
##########################################################
# An alternative script to build a MicrOmegas shared
# library.
#
# Author: Jonathan Cornell
# Date: April 2014, May 2015
#
# Notes: To use this you need to have compiled MicrOmegas 
# and the model you are interested in with -fPIC added 
# to FFLAGS in the FlagsForSh file. Run this in the model
# directory. This script adds all the functions from
# micromegas.a (except Fortran wrappers) to the shared
# library along with all their dependencies. The shared
# object is created in the MODEL/lib directory.
###########################################################

# Variable definitions
CXX=gcc
MODEL=SingletDM
# The below is what your compiler needs to link to X11. For
# OSX it is "-L/usr/X11R6/lib -lX11" and for most Linux
# distributions is it "-lX11" I believe
X11="-L/usr/X11R6/lib -lX11"

mkdir ./so_tmp
cd so_tmp
ar -x ../../sources/micromegas.a

# Removing Fortran Wrappers
rm *Fort.o fortran.o

# And now, build the shared library:
$CXX -shared -o libmicromegas$MODEL.so *.o ../../CalcHEP_src/lib/sqme_aux.so ../work/work_aux.a ../../CalcHEP_src/lib/dynamic_me.a ../lib/aLib.a ../../CalcHEP_src/lib/libSLHAplus.a ../../CalcHEP_src/lib/ntools.a ../../CalcHEP_src/lib/num_c.a ../../CalcHEP_src/lib/serv.a $X11

# Move the shared object to the MODEL/lib directory and clean
# up the mess.
cp libmicromegas$MODEL.so ../lib/
cd ..
rm -rf ./so_tmp
