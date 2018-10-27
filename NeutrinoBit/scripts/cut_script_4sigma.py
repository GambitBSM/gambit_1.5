#!/usr/bin/env python

import h5py
import numpy as np
import sys
import os



#
# Usage: python cut_hdf5.py input.hdf5 output.hdf5
#



# Get input file path
in_file_path = sys.argv[1]
out_file_path = sys.argv[2]


#
# Read file
#


f = h5py.File(in_file_path, 'r')

group = f["/RHN"]

Ue1 = np.array(group["#Ue1 @NeutrinoBit::Ue1"])
Um1 = np.array(group["#Um1 @NeutrinoBit::Um1"])
Ut1 = np.array(group["#Ut1 @NeutrinoBit::Ut1"])

LogLike = np.array(group["LogLike"])
#LogLike_isvalid = np.array(group["LogLike_isvalid"],dtype=np.bool)

print group.keys()

# This applies validity conditions to a few parameters (and we assume the others follow suit)
#mask = LogLike_isvalid

# mask = np.array([True]*len(LogLike))
mask = LogLike > 0

print "Number of valid combined entries: ", np.sum(mask)




#
# Apply cuts
#

# Construct list of cuts
cut_list = []

LL_bf = np.max(LogLike[mask])  # Get best-fit likelihood value     

# # 2 sigma cut
#cut_list.append( (LogLike[mask] > LL_bf - 3.09) )

# # 3 sigma cut
# cut_list.append( (LogLike[mask] > LL_bf - 5.915) )

# 4 sigma cut
cut_list.append( (LogLike[mask] > LL_bf - 9.665) )

# # Other cuts
#cut_list.append( ( M1[mask] > 100) )

# Ue1, Um1, Ut1 > 10^-10
#cut_list.append( (Ue1[mask] > 1e-10 ) )
#cut_list.append( (Um1[mask] > 1e-10 ) )
#cut_list.append( (Ut1[mask] > 1e-10 ) )
#cut_list.append( (Ut1[mask] < 1e-4 ) )

# Construct np.array of bools used to apply the cuts
if len(cut_list) == 0:
    cuts = np.array( [True]*len(LogLike[mask]) , dtype=np.bool)
else:
    cuts = np.array([np.all(l) for l in zip(*cut_list)], dtype=np.bool)

# Print output file
with h5py.File(out_file_path, 'w') as out_f:

    out_group = out_f.create_group('RHN') 

    for dataset_name in group.keys():
        print "Masked ", dataset_name
        masked_set = group[dataset_name][mask][cuts]
        out_group.create_dataset(dataset_name, data=masked_set)


LogLike = LogLike[mask][cuts]

print
print "Number of points after cuts: ", len(LogLike)
if len(LogLike) > 0:
    print "L range: ", np.min(LogLike), np.max(LogLike)

f.close()

