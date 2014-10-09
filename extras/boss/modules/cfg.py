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


# GAMBIT specific options:
# -- NOTE: Setting gambit_mode = True will overwrite some of the other choices below, like 'all_wrapper_fname'.
gambit_mode            = True
gambit_backend_name    = 'Minimal'
gambit_backend_version = '0.1'
gambit_base_namespace  = 'Gambit::Backends'


# Information about the external code:

# accepted_paths     = ['pythia8186_original', 'pythia8186']
accepted_paths     = ['minimal_original', 'minimal']

std_include_paths  = ['/usr/include/']

# loaded_classes       = ['Pythia8::Pythia', 'Pythia8::Event', 'Pythia8::Particle']

# loaded_classes       = ['Pythia8::Pythia', 'Pythia8::Hist', 'Pythia8::Event', 'Pythia8::Particle', 'Pythia8::Info', 'Pythia8::Vec4']


# loaded_classes       = ['Pythia8::ParticleDataEntry', 'Pythia8::ResonanceWidths']

# loaded_classes       = ['Pythia8::Pythia', 'Pythia8::Hist', 'Pythia8::Event', 'Pythia8::Particle', 'Pythia8::Info', 'Pythia8::Vec4',
#                         'Pythia8::Sphericity', 'Pythia8::Thrust', 'Pythia8::ClusterJet', 'Pythia8::ParticleData', 'Pythia8::ParticleDataEntry',
#                         'Pythia8::RotBstMatrix', 'Pythia8::Junction', 'Pythia8::Couplings', 'Pythia8::ResonanceWidths', 'Pythia8::DecayChannel']

# loaded_classes       = ['Pythia8::Pythia', 'Pythia8::Settings', 'Pythia8::Hist', 'Pythia8::Event', 'Pythia8::Particle', 'Pythia8::Info', 'Pythia8::Vec4',
#                         'Pythia8::Sphericity', 'Pythia8::Thrust', 'Pythia8::ClusterJet', 'Pythia8::ParticleData', 'Pythia8::ParticleDataEntry',
#                         'Pythia8::RotBstMatrix', 'Pythia8::Junction', 'Pythia8::Couplings', 'Pythia8::ResonanceWidths', 'Pythia8::DecayChannel']

# loaded_functions     = ['']

loaded_classes       = ['NamespaceForX::X', 'NamespaceForY::Y']
loaded_functions     = []

wrapper_class_tree     = True
load_parent_classes    = False
wrap_inherited_members = False

extra_output_dir        = 'output'
code_suffix             = '_GAMBIT'
abstr_header_prefix     = 'Abstract_'
factory_file_prefix     = 'factory_'
abstr_class_prefix      = 'Abstract_'
wrapper_header_prefix   = 'GAMBIT_wrapper_'

all_wrapper_fname       = 'boss_loaded_classes'
wrapper_deleter_fname   = 'wrapper_deleter'
all_typedefs_fname      = 'all_typedefs'
frwd_decls_abs_fname    = 'forward_decls_abstract_classes'
frwd_decls_wrp_fname    = 'forward_decls_wrapper_classes'
wrapper_typedefs_fname  = 'wrapper_typedefs'

header_extension        = '.hpp'
source_extension        = '.cpp'

add_path_to_includes    = ''  #'Pythia8'


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


# Dictionary of what names to use for various operator symbols

operator_names = {
          "="   : "equal",
          "+"   : "plus",
          "-"   : "minus",
          "*"   : "asterix",
          "/"   : "slash",
          "%"   : "percent",
          "&"   : "ampersand",
          "++"  : "plus_plus", 
          "--"  : "minus_minus",
          "+="  : "plus_equal",
          "-="  : "minus_equal",
          "*="  : "asterix_equal",
          "/="  : "slash_equal",
          "%="  : "percent_equal",
          "&="  : "ampersand_equal",
          "|="  : "bar_equal",
          "^="  : "caret_equal",
         "<<="  : "double_angle_bracket_left_equal",
         ">>="  : "double_angle_bracket_right_equal",
          "[]"  : "square_bracket_pair",
          "()"  : "round_bracket_pair",
          "=="  : "double_equal",
          "!="  : "exclamation_equal",
          ">"   : "angle_bracket_right",
          "<"   : "angle_bracket_left",
          ">="  : "angle_bracket_right_equal",
          "<="  : "angle_bracket_left_equal",
          "!"   : "exclamation",
          "&&"  : "double_ampersand",
          "|"   : "bar",
          "^"   : "caret",
          "<<"  : "double_angle_bracket_left",
          ">>"  : "double_angle_bracket_right",
          "->"  : "arrow",
         "->*"  : "arrow_asterix",
          ","   : "comma",
         "new"  : "new",
       "new[]"  : "new_square_bracket_pair",
      "delete"  : "delete",
    "delete[]"  : "delete_square_bracket_pair",
}


# Some global variables

xml_file_name    = ''
id_dict          = OrderedDict() 
file_dict        = OrderedDict()
std_types_dict   = OrderedDict()
typedef_dict     = OrderedDict()
class_dict       = OrderedDict()
func_dict        = OrderedDict()
new_header_files = OrderedDict()
accepted_types   = [] 
std_headers_used = []

std_types_in_class = {}
all_types_in_class = {}

gambit_full_namespace = ''
gambit_backend_safeversion = ''