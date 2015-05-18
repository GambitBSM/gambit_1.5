###################################
#                                 #
#  Configuration module for BOSS  #
#                                 #
###################################

from collections import OrderedDict



# ~~~~~ GAMBIT-specific options ~~~~~

gambit_backend_name    = 'Pythia'
gambit_backend_version = '8.186'
gambit_base_namespace  = ''


# ~~~~~ Information about the external code ~~~~~

include_path = 'pythia8186/include'
source_path  = 'pythia8186/src'

additional_include_paths = []

accepted_paths     = ['pythia8186']

std_include_paths  = ['/usr/include/']

loaded_classes     = [
                      'Pythia8::Pythia',
                      'Pythia8::Hist',
                      'Pythia8::Event',
                      'Pythia8::Particle',
                      'Pythia8::Info',
                      'Pythia8::Vec4',
                      'Pythia8::Rndm',
                      'Pythia8::SlowJet',
                      'Pythia8::ParticleData',
                      'Pythia8::ParticleDataEntry',
                      'Pythia8::Settings',
                      'Pythia8::SigmaTotal',
                      'Pythia8::SigmaProcess',
                      'Pythia8::PartonLevel',
                      'Pythia8::Couplings',
                      'Pythia8::ResonanceGmZ',
                      'Pythia8::CoupSUSY',
                      'Pythia8::SLHAinterface',
                     ]

loaded_functions   = []

ditch = [
          'Pythia8::Pythia::initSLHA',
        ]

# wrapper_class_tree     = True
load_parent_classes    = True
wrap_inherited_members = False

extra_output_dir      = 'pythia_BOSS_output'
abstr_header_prefix   = 'abstract_'
wrapper_header_prefix = 'wrapper_'
factory_file_prefix   = 'factory_'

header_extension = '.h'
source_extension = '.cc'

add_path_to_includes = 'Pythia8'

indent = 4


# ~~~~~ Dictionary standard type headers ~~~~~

known_class_headers = {
    "std::array"             : "<array>", 
    "std::vector"            : "<vector>", 
    "std::deque"             : "<deque>", 
    "std::complex"           : "<complex>", 
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

