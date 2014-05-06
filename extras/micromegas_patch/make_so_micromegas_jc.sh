#!/bin/sh
##########################################################
# An alternative script to build a micromegas shared
# library.
#
# Author: Jonathan Cornell
# Date: 30 April 2014
#
# Notes: To use this you need to have compiled micromegas 
# and the model you are interested in with -fPIC added 
# to FFLAGS in the FlagsForSh file. Run this in the model
# directory. This script adds all the functions from
# micromegas.a to the shared library along with all their
# dependencies. Note that the library built here depends
# on sqme_aux.so (which is located in 
# MICROMEGAS/CalcHep_src/lib) so make sure that GAMBIT can
# also find this library (by adding the above path to
# $LD_LIBRARY_PATH or whatever). The shared object is
# created in the MODEL/lib directory.
###########################################################

mkdir ./so_tmp
cd so_tmp
ar -x ../../sources/micromegas.a
cp ../../CalcHEP_src/lib/sqme_aux.so ./

# For the MSSM you need to include the following objects from
# aLib.a (I don't include the entire aLib.a library because 
# it leads symbols being multiply defined). More may be
# needed as DarkBit grows.
ar -x ../lib/aLib.a rdLesH.o fillVal.o bsg_nlo.o lambdas.o

# And now, build the shared library:
gcc -shared -o micromegas.so *.o -l:sqme_aux.so -L../../CalcHEP_src/lib -L../work -L../lib -l:work_aux.a -l:dynamic_me.a -l:aLib.a -l:libSLHAplus.a -l:ntools.a -l:num_c.a -l:serv.a

# Move the shared object to the MODEL/lib directory and clean
# up the mess.
cp micromegas.so ../lib/
cd ..
rm -r ./so_tmp
