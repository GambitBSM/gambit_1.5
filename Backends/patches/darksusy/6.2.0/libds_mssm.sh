#!/bin/sh
# Run this file in the darksusy-6.2.0/lib directory to create a shared object for 
# the DarkSUSY MSSM mode (libds_mssm.so)

rm -rf ./libds_mssm
mkdir ./libds_mssm
cd libds_mssm
cp ../libds_mssm.a ./
cp ../libds_core.a ./
cp ../libds_empty.a ./
cp ../libFH.a ./
cp ../libHB.a ./
cp ../libHS.a ./
cp ../libcfitsio.a ./
cp ../libhealpix.a ./
cp ../libisajet.a ./
cp ../libisospin.a ./
for file in ./*.a; do ar x "$file"; done
cd ../
gfortran libds_mssm/*.o -shared -fopenmp -o libds_mssm.so 
rm -rf ./libds_mssm
