#!/bin/bash
#
#SBATCH --nodes 1                              # Number of nodes
#SBATCH --ntasks-per-node 1                # Number of MPI processes per node.
#SBATCH --mem 3G
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM
#SBATCH -p plgrid                          # Partition to submit to
#SBATCH -A hbtv1                           # Account to charge
#SBATCH --signal=SIGUSR1@300               # Signal and time before walltime to send
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mchrzasz@cern.ch       # Email to which notifications will be sent



cd /net/archive/groups/plgghbt/gambit/gambit_RHN/NeutrinoBit/scripts/promek_scripts
module load plgrid/tools/python/2.7.14 plgrid/tools/cmake/3.10.2 plgrid/libs/boost/1.58.0 plgrid/libs/gsl/2.4 plgrid/libs/hdf5/1.8.19 plgrid/tools/intel/17.0.5 plgrid/libs/mkl/11.3.3 plgrid/libs/eigen/3.2.7 plgrid/libs/hdf5/1.8.16


