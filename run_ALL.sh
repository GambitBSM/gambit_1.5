#!/bin/bash

for f in /net/archive/groups/plgghbt/gambit/gambit_RHN/yaml/*.yaml
do
    echo 'Preparing'$f
    f1=`echo $f  | rev | cut -d '/' -f 1 | rev | cut -d '.' -f  1 `
    echo $f1
    cp run.sh run_$f1.sh
    echo '#SBATCH -o runs/RHN1'$f1'.out' >> run_$f1.sh
    echo '#SBATCH -e runs/RHN1'$f1'.err' >> run_$f1.sh
    echo 'module load plgrid/tools/python/2.7.14 plgrid/tools/cmake/3.10.2 plgrid/libs/boost/1.58.0 plgrid/libs/gsl/2.4 plgrid/libs/hdf5/1.8.19 plgrid/tools/intel/17.0.5 plgrid/libs/mkl/11.3.3 plgrid/libs/eigen/3.2.7 plgrid/libs/hdf5/1.8.16'  >> run_$f1.sh
    echo 'cd /net/archive/groups/plgghbt/gambit/gambit_RHN'  >> run_$f1.sh

    if [[ ! -d  'runs/'$f1 ]]; then 
	mkdir 'runs/'$f1;
    fi
    echo 'mpiexec -outfile-pattern=runs/RHN_'$f1'.out -errfile-pattern=runs/RHN'$f1'.err ./gambit -f '$f >> run_$f1.sh

    
    if [[ ! -d  'runs/RHN/'$f1'/samples/RHN.hdf5' ]]; then
	echo 'Submitting '$f1
	sbatch run_$f1.sh
    if

done
