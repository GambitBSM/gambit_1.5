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


accepted_paths     = ['original_src']

std_include_paths  = ['/usr/include/']

#accepted_classes = ['T', 'X']
accepted_classes = ['T']

accepted_functions = ['printT']

extra_output_dir      = 'test_output'
# extra_output_dir      = 'extra_output_softsusy'
code_suffix           = '_GAMBIT'
abstr_header_prefix   = 'abstract_'
factory_file_prefix   = 'factory_'
abstr_class_prefix    = 'Abstract__'
all_headers_fname     = 'all_abstract_headers.hpp'
all_typedefs_fname    = 'all_typedefs.hpp'
header_extension      = '.hpp'
source_extension      = '.cpp'
add_path_to_includes  = ''

# header_extension      = '.h'
# source_extension      = '.cpp'
# add_path_to_includes  = ''

indent = 4