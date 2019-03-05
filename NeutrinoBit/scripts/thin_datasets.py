#!/usr/bin/env python

import h5py
import numpy as np
import sys
import os
import math
import copy
from optparse import OptionParser


######################################################################################################
# Script for thinning datasets for easier plotting.
# It will thin the data and only keep the points with hightest LogLike withing certain resolution in
# a multi-dimensional variable space, either for specific variables provided or all scan parameters.
# WARNING! If used with specific variables, do not use the resulting samples for plotting anything 
# other than those variables, since data will be lost in other dimensions
# 
# Usage: python thin_datasets.py <in_file> <out_file> <group> <resolution> [<var1> <var2>...]
#        <in_file>                   Input file
#        <out_file>                  Output file
#        <group>                     Group inside the dataset
#        <resolution>                Resolution of the thinning, relative to the sampling area
#        [<var1> <var2>...]          Variables to use for thinning
#
# Options:
#        -h, --help                    Show this help message and exit
#        -v, --verbose                 Print debug output
#        -a, --all                     Use all the scan parameters for thinning
#        -l LOGVARS, --logvars=LOGVARS Set log scale for input variables
#
#######################################################################################################

# Parse options and arguments
parser = OptionParser()
parser.add_option("-v", "--verbose", help = "Print debug output", action='store_true', dest='verbose', default=False)
parser.add_option("-a", "--all", help = "Use all the scan parameters for thinning", action='store_true', dest='all_params', default=False)
parser.add_option("-l", "--logvars", help = "Set log scale for input variables", metavar='LOGVARS')
(options, args) = parser.parse_args()
if len(args) < 4 or (not options.all_params and len(args) < 6):
  print 'Wrong number of arguments \n\
        \n\
Usage: python thin_datasets.py <in_file> <out_file> <group> <resolution> [<var1> <var2>...]\n\
       <in_file>                     Input file\n\
       <out_file>                    Output file\n\
       <group>                       Group inside the dataset\n\
       <resolution>                  Resolution of the thinning, relative to the sampling area\n\
       [<var1> <var2>...]            Variables to use for thinningi\n\
\n\
Options:\n\
       -h, --help                    Show this help message and exit\n\
       -v, --verbose                 Print debug output\n\
       -a, --all                     Use all the scan parameters for thinning\n\
       -l LOGVARS, --logvars=LOGVARS Set log scale for input variables'

  exit()

in_file_path = args[0]
out_file_path = args[1]
in_group = args[2]
resolution = float(args[3])
if not options.all_params:
  param_names = args[4:]
  nparams = len(args)-4

# Read file
f = h5py.File(in_file_path, 'r')

group = f[in_group]

params = []
if options.all_params:
  nparams = 0
  while "unitCubeParameters["+str(nparams)+"]" in group :
    params.append(np.array(group["unitCubeParameters["+str(nparams)+"]"]))
    nparams+=1
else :
  for param_name in param_names:
    params.append(np.array(group[param_name]))

# Set logscales for the params in needed
if options.logvars is not None :
  logvars = options.logvars.split(',')
  for logvar in logvars:
    params[int(logvar)] = np.log10(params[int(logvar)])

# Identify limits and binning
mins = []
maxs = []
bins = []
r = 0
for param in params :
  mins.append(np.min(param))
  maxs.append(np.max(param))
  bins.append(resolution * (maxs[-1] - mins[-1]))
  r = r + bins[-1]**2
r = np.sqrt(r)

mask =  np.zeros(len(params[0]), dtype=bool)

print "Number of original entries: ",len(mask)
if options.verbose: print "Preparing mask..."

params_full = copy.deepcopy(params)

# Take LogLike dataset for sorting
LogLike = np.array(group['LogLike'])
LLsort = LogLike.argsort()
params[:] = [param[LLsort[::-1]] for param in params]

count = 0
while len(params[0]) :
  rad = np.array( np.sum( [(param - param[0])**2 for param in params], axis=0 ) )
  mask = mask | np.all( [params_full[i] == params[i][0] for i in range(nparams)], axis=0)
  params[:] = [param[rad > r] for param in params]

  count += 1
  if count == 100 :
    if options.verbose: print len(params[0]), " points left ..."
    count = 0
 

if options.verbose: print "Number of thinned entries: ", np.sum(mask)

# Print output file
with h5py.File(out_file_path, 'w') as out_f:

    out_group = out_f.create_group(in_group)

    for i,dataset_name in enumerate(group.keys()):
        if options.verbose: print "Masked ", dataset_name
        old_dataset = group[dataset_name]
        masked_set = old_dataset[mask]
        out_group.create_dataset(dataset_name, data=masked_set)


LogLike = LogLike[mask]

print
print "Number of points after cuts: ", len(LogLike)
if len(LogLike) > 0:
    print "L range: ", np.min(LogLike), np.max(LogLike)

f.close()
