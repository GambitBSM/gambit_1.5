###################################
#                                 #
#  Configuration module for BOSS  #
#                                 #
###################################

from collections import OrderedDict



# GAMBIT specific options:

gambit_backend_name    = 'BOSSedPythia'
gambit_backend_version = '1.0'
gambit_base_namespace  = ''
gambit_backend_basedir = 'backend_types'

shared_lib_file_name = 'libpythia8.so'

# Information about the external code:

accepted_paths     = ['pythia8186_original', 'pythia8186']
std_include_paths  = ['/usr/include/']

loaded_classes     = ['Pythia8::Pythia', 'Pythia8::Hist', 'Pythia8::Event', 'Pythia8::Particle', 'Pythia8::Info', 'Pythia8::Vec4']
loaded_functions   = []

wrapper_class_tree     = True
load_parent_classes    = False
wrap_inherited_members = False

extra_output_dir      = 'pythia_BOSS_output'
abstr_header_prefix   = 'abstract_'
wrapper_header_prefix = 'wrapper_'
factory_file_prefix   = 'factory_'

header_extension = '.h'
source_extension = '.cc'

add_path_to_includes = 'Pythia8'

indent = 4


# Dictionary of what header to include for various standard types

known_class_headers = {
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
    "std::priority_queue"    : "<queue>",
    "std::string"            : "<string>",
    "std::istream"           : "<istream>",
    "std::ostream"           : "<ostream>",
    "std::iostream"          : "<iostream>"
}

