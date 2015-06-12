###############################
#                             #
#  Global variables for BOSS  #
#                             #
###############################

#
# These variables can be changed during program execution.
#

from collections import OrderedDict
import os

import modules.cfg as cfg

boss_temp_dir = 'temp'

new_code              = OrderedDict()
xml_file_name         = ''

id_dict               = OrderedDict()
all_id_dict           = OrderedDict()
name_dict             = OrderedDict()
all_name_dict         = OrderedDict()

file_dict             = OrderedDict()
all_classes_dict      = OrderedDict()
std_types_dict        = OrderedDict()
typedef_dict          = OrderedDict()
class_dict            = OrderedDict()
loaded_classes_in_xml = OrderedDict()
func_dict             = OrderedDict()
new_header_files      = OrderedDict()
accepted_types        = [] 
std_headers_used      = []
original_file_paths   = OrderedDict()
functions_done        = []
n_functions_done      = 0
classes_done          = []
n_classes_done        = 0
factory_info          = OrderedDict()

gambit_backend_namespace    = 'CAT_3(BACKENDNAME,_,SAFE_VERSION)'
gambit_backend_safeversion  = cfg.gambit_backend_version.replace('.','_')
gambit_backend_name_full    = cfg.gambit_backend_name + '_' + gambit_backend_safeversion

gambit_backend_types_basedir = 'backend_types'
gambit_backend_dir_complete  = os.path.join(cfg.extra_output_dir, gambit_backend_types_basedir, gambit_backend_name_full) 

gambit_backend_incl_dir = 'gambit/Backends/'
gambit_utils_incl_dir   = 'gambit/Utils/'

code_suffix         = '__BOSS'
abstr_class_prefix  = 'Abstract_'


func_return_utils_fname = 'function_return_utils'
wrapper_deleter_fname   = 'wrapperdeleter'
frwd_decls_abs_fname    = 'forward_decls_abstract_classes'
frwd_decls_wrp_fname    = 'forward_decls_wrapper_classes'
wrapper_typedefs_fname  = 'wrappertypedefs'
abstract_typedefs_fname = 'abstracttypedefs'
enum_decls_wrp_fname    = 'enum_decl_copies'


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
