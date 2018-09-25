#!/bin/bash
#
#SBATCH --nodes 20                              # Number of nodes
#SBATCH --ntasks-per-node 12                # Number of MPI processes per node.
#SBATCH --mem 120G
#SBATCH -t 0-23:50                         # Runtime in D-HH:MM
#SBATCH -p plgrid                          # Partition to submit to
#SBATCH -A hbtv0                           # Account to charge
#SBATCH --signal=SIGUSR1@300               # Signal and time before walltime to send
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mchrzasz@cern.ch       # Email to which notifications will be sent

