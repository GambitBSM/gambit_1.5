#!/bin/bash



for f in /net/archive/groups/plgghbt/gambit/gambit_RHN/runs/RHN_inv/*/samples/RHN.hdf5
do
    echo 'Doing' $f
    #f1=`echo $f  | rev | cut -d '/' -f 1 | rev | cut -d '.' -f  1 `
    #echo $f1
    f1='RHN_noslide.hdf5'
    f2='RHN_4sigma.hdf5'
    
    fout=$f
    fout=${fout/'RHN.hdf5'/$f1}
    echo $fout
    
    fout2=$f
    fout2=${fout2/'RHN.hdf5'/$f2}
    echo $fout2

    name=`echo $f  | rev | cut -d '/' -f 3 | rev | cut -d '.' -f  1 ` 
    echo $name
    
    cp template.sh lunch_inv/$name.sh
    echo 'python /net/archive/groups/plgghbt/gambit/gambit_RHN/NeutrinoBit/scripts/cut_script_noslide.py '$f $fout >> lunch_inv/$name.sh
    echo 'python /net/archive/groups/plgghbt/gambit/gambit_RHN/NeutrinoBit/scripts/cut_script_4sigma.py '$fout $fout2 >> lunch_inv/$name.sh
    
    sbatch lunch_inv/$name.sh

    

done
