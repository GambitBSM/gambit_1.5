####################################
#                                  #
#  Utility functions for handling  # 
#  C++ functions with BOSS         #
#                                  #
####################################

import xml.etree.ElementTree as ET
from collections import OrderedDict
import copy
import warnings

import modules.cfg as cfg
import modules.utils as utils


# ======== getArgs ========

def getArgs(func_el):
    
    #
    # Returns a list with one dict per argument.
    # Each dict contains the following keywords:
    # 
    #   'name', 'type', 'kw', 'id', 'native', 'fundamental', 'namespaces', 'accepted_class'
    #

    args = []
    argc = 1
    for sub_el in func_el.getchildren():
        if sub_el.tag == 'Argument':

            arg_dict = OrderedDict()
            if 'name' in sub_el.keys():
                arg_dict['name'] = sub_el.get('name')
            else:
                arg_dict['name'] = 'arg_' + str(argc)
                argc += 1
            arg_dict['type'], arg_dict['kw'], arg_dict['id'] = utils.findType(sub_el)

            arg_type_el = cfg.id_dict[arg_dict['id']]
            arg_dict['native']      = utils.isNative(arg_type_el)
            arg_dict['fundamental'] = utils.isFundamental(arg_type_el)

            arg_dict['accepted_class'] = False
            if arg_type_el.tag == 'Class':
                demangled_name = arg_type_el.get('demangled').replace(' ','')
                if demangled_name in cfg.accepted_classes:
                    arg_dict['accepted_class'] = True

            type_el = cfg.id_dict[arg_dict['id']]
            arg_dict['type_namespaces'] = utils.getNamespaces(type_el)

            args.append(arg_dict)

    return args

# ======== END: getArgs ========



# ======== constrArgsBracket ========

def constrArgsBracket(args, include_arg_name=True, include_arg_type=True, include_namespace=False): #invent_arg_names=False):

    #
    # Requires a list of dicts as input, as returned by 'getArgs' or 'constrWrapperArgs'.
    #

    # Construct bracket with input arguments
    args_seq = ''
    argc = 1
    for arg_dict in args:

        if include_arg_type:
            args_seq += ''.join([ kw+' ' for kw in arg_dict['kw'] ])
            if include_namespace:
                namespaces = arg_dict['type_namespaces']
                if len(namespaces)>0:
                    args_seq += '::'.join(namespaces) + '::'
            args_seq += arg_dict['type']
        if include_arg_name:
            args_seq += ' ' + arg_dict['name']

        args_seq += ', '
    args_seq = args_seq.rstrip(', ')
    args_bracket = '(' + args_seq + ')'

    return args_bracket

# ======== END: constrArgsBracket ========



# ======== constrWrapperName ========

def constrWrapperName(func_el):

    func_name   = func_el.get('name')
    w_func_name = func_name + cfg.code_suffix

    return w_func_name

# ======== END: constrWrapperName ========



# # ======== chooseWrapperReturnType ========

# def chooseWrapperReturnType(return_el):

#     # return_type, return_kw, return_id = utils.findType( cfg.id_dict[func_el.get('returns')] )
#     # return_el = cfg.id_dict[return_id]

#     if utils.isNative(return_el):
#         w_return_type = cfg.abstr_class_prefix + return_type
#     else:
#         w_return_type = return_type

#     return w_return_type

# # ======== END: chooseWrapperReturnType ========



# ======== constrWrapperArgs ========

def constrWrapperArgs(args):

    #
    # Requires a list of dicts as input, as returned by 'getArgs'.
    #

    import modules.classutils as classutils

    w_args = copy.deepcopy(args) 
    for arg_dict in w_args:
        if arg_dict['native'] == True:
            if len(arg_dict['type_namespaces']) > 0:
                # namespaces, type_name = arg_dict['type'].rsplit('::',1)
                # arg_dict['type'] = namespaces + '::' + cfg.abstr_class_prefix + type_name
                arg_dict['type'] = classutils.getAbstractClassName(arg_dict['type'])
            else:
                arg_dict['type'] = cfg.abstr_class_prefix + arg_dict['type']

    return w_args

# ======== END: constrWrapperArgs ========



# ======== constrDeclLine ========

def constrDeclLine(return_type, func_name, args_bracket, keywords=[]):

    decl_line = ''
    for keyw in keywords: 
        decl_line += keyw + ' '
    decl_line += return_type + ' ' + func_name + args_bracket

    return decl_line

# ======== END: constrDeclLine ========



# ======== constrWrapperBody ========

def constrWrapperBody(return_type, func_name, args, return_is_native, keywords=[]):

    #
    # Input:
    # 
    #   - Return type of original function
    #   - List of dicts for original arguments
    #   - Name of original function
    #   - Boolean stating whether the orignal return type is native

    w_func_body = ''
    call_args = []

    for arg_index in range( len(args) ):
        
        arg_dict   = args[arg_index]
        # w_arg_dict = w_args[arg_index]  # Not used?

        orig_type        = arg_dict['type']
        n_pointers       = orig_type.count('*')
        is_ref           = bool(orig_type.count('&'))
        orig_type_base   = orig_type.replace('*','').replace('&','')
        name             = arg_dict['name']
        call_arg         = name

        # Construct sequence of temp variables needed for conversion of native types
        if arg_dict['native'] == True:

            kw_str  = ' '.join(arg_dict['kw'])        
            kw_str += ' '*bool(len(kw_str))

            if n_pointers == 0:
                temp_var = '_temp_'+name
                # w_func_body += ' '*cfg.indent + orig_type_base + ' ' + temp_var + ' = ' + '*(' + name + '.downcast())' + ';\n'
                w_func_body += ' '*cfg.indent + kw_str + orig_type_base + '&'*is_ref + ' ' + temp_var + ' = *(reinterpret_cast<' + kw_str + orig_type_base + '*>(&' + name + '));\n'
                call_arg = temp_var
            else:
                temp_var = '_temp_'+name
                w_func_body += ' '*cfg.indent + kw_str + orig_type_base + '*'*n_pointers + ' ' + temp_var + ' = reinterpret_cast<' + kw_str + orig_type_base + '*'*n_pointers + '>(' + name + ')' + ';\n'
                call_arg = temp_var
                # for i in range(1, n_pointers+1):
                #     temp_var = '_temp_'+name+'_'+str(i)
                #     call_arg = temp_var
                #     if i == 1:
                #         w_func_body += ' '*cfg.indent + orig_type_base + '* ' + temp_var + ' = ' + '(' + n_pointers*'*' + name + ').downcast()' + ';\n'
                #     else:
                #         prev_temp_var = '_temp_'+name+'_'+str(i-1)
                #         w_func_body += ' '*cfg.indent + orig_type_base + i*'*' + ' ' + temp_var + ' = &' + prev_temp_var + ';\n'
                # w_func_body += '\n'

        # Update list of argument names used for original function call
        call_args.append(call_arg)

    temp_res_name     = '_temp_result'
    orig_func_call    = func_name + '(' + ','.join(call_args) + ')'
    return_type_base  = return_type.replace('*','').replace('&','')
    n_pointers_return = return_type.count('*')   
    is_ref            = bool(return_type.count('&'))
    kw_str            = ' '.join(keywords) + ' '*bool(len(keywords))
    
    # if return_is_native:
    #     w_return_type = cfg.abstr_class_prefix + return_type
    # else:
    #     w_return_type = return_type

    # w_return_type_base  = w_return_type.replace('*','').replace('&','')
    # n_pointers_w_return = w_return_type.count('*')

    if return_type == 'void':
        w_func_body += ' '*cfg.indent + orig_func_call + ';\n'
    else:
        if n_pointers_return == 0:
            if is_ref:
                w_func_body += ' '*cfg.indent + kw_str + return_type_base + '* _temp_p = new ' + return_type_base + '(' + orig_func_call + ');\n'
                w_func_body += ' '*cfg.indent + kw_str + return_type_base + '*& _temp_pr = _temp_p;\n'
                w_func_body += ' '*cfg.indent + 'return _temp_pr;\n'

                # Pythia8::Particle*  temp_p   = new Pythia8::Particle(front());
                # Pythia8::Particle*& temp_pr  = temp_p;
                # return temp_pr;

            else:                
                w_func_body += ' '*cfg.indent + 'return new ' + return_type_base + '(' + orig_func_call + ');\n'
        else:
            w_func_body += ' '*cfg.indent + 'return ' + orig_func_call + ';\n'

            # w_func_body += ' '*cfg.indent + 'return reinterpret_cast<' + w_return_type_base + '*'*n_pointers_w_return + '>(' + orig_func_call + ');\n'

        # w_func_body += ' '*cfg.indent + return_type + ' ' + temp_res_name + ' = ' + func_name + '(' + ','.join(call_args) + ')' ';\n'
        # if return_is_native:
        #     w_func_body += ' '*cfg.indent + 'return ' + temp_res_name + '.abstractify()' + ';\n'
        #     #
        #     # FIXME: Take 'pointerness' of return type into account
        #     #
        # else:
        #     w_func_body += ' '*cfg.indent + 'return ' + temp_res_name + ';\n'
    w_func_body = '{\n' + w_func_body + '}'

    return w_func_body

# ======== END: constrWrapperBody ========



# ======== ignoreFunction ========

def ignoreFunction(func_el):

    return_type_accepted = True
    arg_types_accepted   = True

    # Check function return type
    if not utils.isAcceptedType(func_el):
        return_type, return_kw, return_id = utils.findType(func_el)
        print 'non-accepted return type: ' + return_type
        return_type_accepted = False

    # Check argument types
    args = getArgs(func_el)
    for arg_dict in args:
        arg_type_name = arg_dict['type'].replace(' ','')
        if not utils.isAcceptedType(arg_type_name, byname=True):
            print 'non-accepted argument type: ', arg_type_name
            arg_types_accepted = False
            break

    # Return result
    if (not return_type_accepted) or (not arg_types_accepted):
        return True
    else:
        return False 

# ======== END: ignoreFunction ========



# ======== usesNativeType ========

def usesNativeType(func_el):

    uses_native_type = False

    return_type_name, return_kw, return_id = utils.findType(func_el)
    return_is_native = utils.isNative( cfg.id_dict[return_id] )

    args = getArgs(func_el)
    is_arg_native = [arg_dict['native'] for arg_dict in args]

    if (return_is_native) or (True in is_arg_native):
        uses_native_type = True

    return uses_native_type

# ======== END: usesNativeType ========