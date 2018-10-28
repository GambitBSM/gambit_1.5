#!/bin/bash            

for f in /net/archive/groups/plgghbt/gambit/gambit_RHN/runs/RHN_inv/*/samples/RHN_4sigma.hdf5
do
    echo 'Doing' $f
    name=`echo $f  | rev | cut -d '/' -f 3 | rev | cut -d '.' -f  1 `
    echo $name
    
    f1='RHN_4sigma_thin.hdf5'
    

    fout=$f
    fout=${fout/'RHN_4sigma.hdf5'/$f1}
    
    echo $f $fout

    case=''
    # looking at the case:
    if [[ $f == *runs/RHN_inv/e* ]];
    then
	case='"#Ue1 @NeutrinoBit::Ue1"'
    fi
    if [[ $f == *runs/RHN_inv/m* ]];
    then
	case='"#Um1 @NeutrinoBit::Um1"'
    fi
    if [[ $f == *runs/RHN_inv/t* ]];
    then
	case='"#Ut1 @NeutrinoBit::Ut1"'
    fi

    cp template_long.sh lunch_thin_inv/$name.sh
    echo 'python /net/archive/groups/plgghbt/gambit/gambit_RHN/NeutrinoBit/scripts/thin_datasets.py -v '$f $fout '/RHN "#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::M_1" ' $case 1e-5   >> lunch_thin_inv/$name.sh
    
    sbatch lunch_thin_inv/$name.sh


done
