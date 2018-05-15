#!/bin/bash

for i in 123 132 213 231 312 321; do qsub -v OUTPATH=RHN_diff_NH_${i},RORDER=${i} lisa_RHN; done
