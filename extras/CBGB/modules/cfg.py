#
#   ================================
#   |                              |
#   |   CBGB configuration module  |
#   |                              |
#   ================================
#


# Set path to source file:
#src_file_path = 'example/sdecay.f'
src_file_path = 'example/dsmssm.h'
#src_file_path = '../susyhit_temps/hdecay.f'


# Choose 'fixed' or 'free' format:
format = 'fixed'  


# List the common blocks to be loaded:
load_common_blocks = ['mspctm','pacodes','widths','intdof','vrtxs',
                      'smruseful','smcuseful','couplingconstants','sckm',
                      'mixing','mssmtype','mssmpar','mssmswitch',
                      'sfermionmass','mssmwidths','mssmmixing']


# Convert tabs to how many spaces?
tabs_to_n_spaces = 6


# Settings for constructing a GAMBIT capability corresponding to each common block:
capability_prefix = 'cb_'
capability_suffix = ''
