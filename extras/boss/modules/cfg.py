##########################################
#                                        #
#  Global configuration module for BOSS  #
#                                        #
##########################################

#
# Variables in this module can be accessed and altered by 
# all other modules. The values listet here are only default
# values.
#

from collections import OrderedDict


xml_file_name    = ''
id_dict          = OrderedDict() 
file_dict        = OrderedDict()
std_types_dict   = OrderedDict()
typedef_dict     = OrderedDict()
class_dict       = OrderedDict()
func_dict        = OrderedDict()
accepted_types   = [] 
std_headers_used = []

std_types_used     = []
all_types_in_class = []


accepted_paths     = ['pythia8185']
# accepted_paths     = ['softsusy']

std_include_paths  = ['/usr/include/']

accepted_classes   = ['Pythia8::ParticleData', 'Pythia8::RotBstMatrix', 
                      'Pythia8::Vec4', 'Pythia8::Particle', 
                      'Pythia8::Event', 'Pythia8::Pythia']

# accepted_classes   = ['Pythia8::Flag', 'Pythia8::Mode', 
#                       'Pythia8::Parm', 'Pythia8::Word', 
#                       'Pythia8::FVec', 'Pythia8::MVec',
#                       'Pythia8::PVec', 'Pythia8::Settings']


# accepted_classes   = ['DoubleVector']
 # - ComplexMatrix
 # - OpMultiply<Complex>
 # - AltEwsbMssm
 # - ComplexVector
 # - DoubleMatrix
 # - Indexable<Complex,ComplexVector>
 # - Indexable<double,DoubleVector>
 # - QedQcd
 # - SoftParsMssm
 # - MssmSoftsusy
 # - MatIndexable<Complex,ComplexMatrix>
 # - MatIndexable<double,DoubleMatrix>
 # - Complex
 # - MssmSusy
 # - RGE

accepted_functions = []

extra_output_dir      = 'test_output'
# extra_output_dir      = 'extra_output_softsusy'
code_suffix           = '_GAMBIT'
abstr_header_prefix   = 'abstract_'
factory_file_prefix   = 'factory_'
abstr_class_prefix    = 'Abstract__'
all_headers_fname     = 'all_abstract_headers.hpp'
all_typedefs_fname    = 'all_typedefs.hpp'
header_extension      = '.h'
source_extension      = '.cc'
add_path_to_includes  = 'Pythia8'

# header_extension      = '.h'
# source_extension      = '.cpp'
# add_path_to_includes  = ''

indent = 4