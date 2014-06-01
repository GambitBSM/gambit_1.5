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

std_types_in_class = {}
all_types_in_class = {}


accepted_paths     = ['classloader_example/original']

std_include_paths  = ['/usr/include/']

loaded_classes       = ['T', 'X']
loaded_functions     = []

wrapper_class_tree     = True
load_parent_classes    = True
wrap_inherited_members = False

extra_output_dir      = 'test_output'
code_suffix           = '_GAMBIT'
abstr_header_prefix   = 'abstract_'
factory_file_prefix   = 'factory_'
abstr_class_prefix    = 'Abstract__'
wrapper_header_prefix = 'GAMBIT_wrapper_'
all_headers_fname     = 'all_abstract_headers.hpp'
all_typedefs_fname    = 'all_typedefs.hpp'
header_extension      = '.hpp'
source_extension      = '.cpp'
add_path_to_includes  = ''


indent = 4


# Dictionary of what header to include for various standard types
class_header_dict = {
 	"std::array"             : "<array>", 
 	"std::vector"            : "<vector>", 
 	"std::deque"             : "<deque>", 
 	"std::list"              : "<list>", 
 	"std::forward_list"      : "<forward_list>", 
 	"std::set"               : "<set>",  
 	"std::multiset"          : "<set>", 
 	"std::map"               : "<map>", 
 	"std::multimap"          : "<map>", 
 	"std::unordered_set"     : "<unordered_set>", 
 	"std::unordered_multiset": "<unordered_set>", 
 	"std::unordered_map"     : "<unordered_map>", 
 	"std::unordered_multimap": "<unordered_map>", 
 	"std::stack"             : "<stack>", 
 	"std::queue"             : "<queue>",
 	"std::priority_queue"    : "<queue>"
}