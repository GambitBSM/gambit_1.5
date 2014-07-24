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
    #   'name', 'type', 'kw', 'id', 'native', 'fundamental', 'namespaces', 'loaded_class', 'type_namespaces'
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
            arg_dict['loaded_class'] = utils.isLoadedClass(arg_type_el)

            type_el = cfg.id_dict[arg_dict['id']]
            arg_dict['type_namespaces'] = utils.getNamespaces(type_el)

            args.append(arg_dict)

    return args

# ======== END: getArgs ========



# # ======== constrArgsBracket ========

# def constrArgsBracket(args, include_arg_name=True, include_arg_type=True, include_namespace=False, use_wrapper_class=False, wrapper_to_pointer=False): #invent_arg_names=False):

#     #
#     # Requires a list of dicts as input, as returned by 'getArgs' or 'constrWrapperArgs'.
#     #

#     import modules.classutils as classutils

#     # Construct bracket with input arguments
#     args_seq = ''
#     argc = 1
#     for arg_dict in args:

#         if include_arg_type:
#             args_seq += ''.join([ kw+' ' for kw in arg_dict['kw'] ])

#             if use_wrapper_class and arg_dict['loaded_class'] == True:
#                 args_seq += classutils.toWrapperType(arg_dict['type'])

#             else:
#                 if include_namespace:
#                     namespaces = arg_dict['type_namespaces']

#                     if len(namespaces)>0:
#                         args_seq += '::'.join(namespaces) + '::'

#                 args_seq += arg_dict['type']

#         if include_arg_type and include_arg_name:
#             args_seq += ' '

#         if include_arg_name:
#             if utils.isLoadedClass(arg_dict['type'], byname=True) and wrapper_to_pointer:
#                 if arg_dict['type'].count('*') == 0:
#                     args_seq += '*' + arg_dict['name'] + '.BEptr'
#                 elif arg_dict['type'].count('*') == 1:
#                     args_seq += arg_dict['name'] + '.BEptr'
#                 else:
#                     raise Exception('funcutils.constrArgsBracket cannot presently deal with arguments of type pointer-to-pointer for wrapper classes.')
#             else:
#                 args_seq += arg_dict['name']

#         args_seq += ', '

#     args_seq = args_seq.rstrip(', ')
#     args_seq = args_seq.strip()
#     args_bracket = '(' + args_seq + ')'

#     return args_bracket

# # ======== END: constrArgsBracket ========



# ======== constrArgsBracket ========

def constrArgsBracket(args, include_arg_name=True, include_arg_type=True, include_namespace=False,
                      cast_to_original=False, use_wrapper_class=False, wrapper_to_pointer=False, 
                      use_wrapper_base_class=False):

    #
    # Requires a list of dicts as input, as returned by 'getArgs' or 'constrWrapperArgs'.
    #

    import modules.classutils as classutils

    # Construct bracket with input arguments
    args_seq = ''
    argc = 1
    for arg_dict in args:

        if include_arg_name and cast_to_original:

            if arg_dict['loaded_class']:

                # We assume that arg_dict['type'] *is* the original type!
                cast_to_type = arg_dict['type']
        
                if include_namespace:
                    namespaces = arg_dict['type_namespaces']
                    if len(namespaces)>0:
                        cast_to_type = '::'.join(namespaces) + '::' + cast_to_type

                # If argument type is not pointer or reference, add a reference operator '&'
                check_type = cast_to_type.split('<')[0]
                if ('*' not in check_type) and ('&' not in check_type):
                    cast_to_type = cast_to_type + '&'

                # Add qualifiers
                if len(arg_dict['kw']) > 0:
                    qualifiers = ' '.join(arg_dict['kw'])
                    cast_to_type = qualifiers + ' ' + cast_to_type

                # args_seq += 'dynamic_cast< ' + cast_to_type + ' >(' + arg_dict['name'] + ')'
                args_seq += 'reinterpret_cast< ' + cast_to_type + ' >(' + arg_dict['name'] + ')'  # Switched to reinterpret cast to be able to cast forward declared types

            else:

                args_seq += arg_dict['name']


        else:
            if include_arg_type:
                args_seq += ''.join([ kw+' ' for kw in arg_dict['kw'] ])

                if use_wrapper_class and arg_dict['loaded_class'] == True:
                    args_seq += classutils.toWrapperType(arg_dict['type'], use_base_type=use_wrapper_base_class)

                else:
                    if include_namespace:
                        args_seq += arg_dict['type']
                    else:
                        args_seq += utils.removeNamespace(arg_dict['type'])
                        # namespaces = arg_dict['type_namespaces']
                        # if len(namespaces)>0:
                        #     args_seq += '::'.join(namespaces) + '::'


            if include_arg_type and include_arg_name:
                args_seq += ' '

            if include_arg_name:
                if utils.isLoadedClass(arg_dict['type'], byname=True) and wrapper_to_pointer:
                    if arg_dict['type'].count('*') == 0:
                        args_seq += '*' + arg_dict['name'] + '.BEptr'
                    elif arg_dict['type'].count('*') == 1:
                        args_seq += '(*' + arg_dict['name'] + ')' + '.BEptr'
                    else:
                        raise Exception('funcutils.constrArgsBracket cannot presently deal with arguments of type pointer-to-pointer for wrapper classes.')
                else:
                    args_seq += arg_dict['name']

        args_seq += ', '

    args_seq = args_seq.rstrip(', ')
    args_seq = args_seq.strip()
    args_bracket = '(' + args_seq + ')'

    return args_bracket

# ======== END: constrArgsBracket ========



# ======== constrWrapperName ========

def constrWrapperName(func_el):

    # Check if this is an operator function
    is_operator = False
    if func_el.tag == 'OperatorMethod':
        is_operator = True

    func_name = func_el.get('name')

    if is_operator:
        w_func_name = 'operator_' + cfg.operator_names[func_name] + cfg.code_suffix
    else:
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

def constrWrapperArgs(args, add_ref=False):

    #
    # Requires a list of dicts as input, as returned by 'getArgs'.
    #

    import modules.classutils as classutils

    # Copy input list
    w_args = copy.deepcopy(args) 

    # The dict entry 'id' does not make sense for arguments that are translated from 
    # native to abstract type
    for arg_dict in w_args:
        del arg_dict['id']

    for arg_dict in w_args:
        if arg_dict['native'] == True:
            if arg_dict['loaded_class']:

                arg_dict['type'] = classutils.toAbstractType(arg_dict['type'])

                # if len(arg_dict['type_namespaces']) > 0:
                #     # namespaces, type_name = arg_dict['type'].rsplit('::',1)
                #     # arg_dict['type'] = namespaces + '::' + cfg.abstr_class_prefix + type_name
                #     arg_dict['type'] = classutils.getAbstractClassName(arg_dict['type'])
                # else:
                #     arg_dict['type'] = cfg.abstr_class_prefix + arg_dict['type']

                if add_ref:
                    if ('&' not in arg_dict['type']) and ('*' not in arg_dict['type']):
                        arg_dict['type'] = arg_dict['type'] + '&'

            else:
                warnings.warn('The argument "%s" is of a native type "%s" that BOSS is not parsing. The function using this should be ignored.' % (arg_dict['name'], arg_dict['type']))

    return w_args

# ======== END: constrWrapperArgs ========



# ======== constrDeclLine ========

def constrDeclLine(return_type, func_name, args_bracket, keywords=[], is_const=False):

    decl_line = ''
    for keyw in keywords: 
        decl_line += keyw + ' '
    decl_line += return_type + ' ' + func_name + args_bracket

    if is_const:
        decl_line = decl_line + ' const'

    return decl_line

# ======== END: constrDeclLine ========



# # ======== constrWrapperBody ========

# def constrWrapperBody(return_type, func_name, args, return_is_native, keywords=[]):

#     #
#     # Input:
#     # 
#     #   - Return type of original function
#     #   - List of dicts for original arguments
#     #   - Name of original function
#     #   - Boolean stating whether the orignal return type is native

#     w_func_body = ''
#     call_args = []

#     for arg_index in range( len(args) ):
        
#         arg_dict   = args[arg_index]
#         # w_arg_dict = w_args[arg_index]  # Not used?

#         orig_type        = arg_dict['type']
#         n_pointers       = orig_type.count('*')
#         is_ref           = bool(orig_type.count('&'))
#         orig_type_base   = orig_type.replace('*','').replace('&','')
#         name             = arg_dict['name']
#         call_arg         = name

#         # Construct sequence of temp variables needed for conversion of native types
#         if arg_dict['native'] == True:

#             kw_str  = ' '.join(arg_dict['kw'])        
#             kw_str += ' '*bool(len(kw_str))

#             if n_pointers == 0:
#                 temp_var = '_temp_'+name
#                 # w_func_body += ' '*cfg.indent + orig_type_base + ' ' + temp_var + ' = ' + '*(' + name + '.downcast())' + ';\n'
#                 w_func_body += ' '*cfg.indent + kw_str + orig_type_base + '&'*is_ref + ' ' + temp_var + ' = *(reinterpret_cast<' + kw_str + orig_type_base + '*>(&' + name + '));\n'
#                 call_arg = temp_var
#             else:
#                 temp_var = '_temp_'+name
#                 w_func_body += ' '*cfg.indent + kw_str + orig_type_base + '*'*n_pointers + ' ' + temp_var + ' = reinterpret_cast<' + kw_str + orig_type_base + '*'*n_pointers + '>(' + name + ')' + ';\n'
#                 call_arg = temp_var
#                 # for i in range(1, n_pointers+1):
#                 #     temp_var = '_temp_'+name+'_'+str(i)
#                 #     call_arg = temp_var
#                 #     if i == 1:
#                 #         w_func_body += ' '*cfg.indent + orig_type_base + '* ' + temp_var + ' = ' + '(' + n_pointers*'*' + name + ').downcast()' + ';\n'
#                 #     else:
#                 #         prev_temp_var = '_temp_'+name+'_'+str(i-1)
#                 #         w_func_body += ' '*cfg.indent + orig_type_base + i*'*' + ' ' + temp_var + ' = &' + prev_temp_var + ';\n'
#                 # w_func_body += '\n'

#         # Update list of argument names used for original function call
#         call_args.append(call_arg)

#     temp_res_name     = '_temp_result'
#     orig_func_call    = func_name + '(' + ','.join(call_args) + ')'
#     return_type_base  = return_type.replace('*','').replace('&','')
#     n_pointers_return = return_type.count('*')   
#     is_ref            = bool(return_type.count('&'))
#     kw_str            = ' '.join(keywords) + ' '*bool(len(keywords))
    
#     # if return_is_native:
#     #     w_return_type = cfg.abstr_class_prefix + return_type
#     # else:
#     #     w_return_type = return_type

#     # w_return_type_base  = w_return_type.replace('*','').replace('&','')
#     # n_pointers_w_return = w_return_type.count('*')

#     if return_type == 'void':
#         w_func_body += ' '*cfg.indent + orig_func_call + ';\n'
#     else:
#         if n_pointers_return == 0:
#             if is_ref:
#                 w_func_body += ' '*cfg.indent + kw_str + return_type_base + '* _temp_p = new ' + return_type_base + '(' + orig_func_call + ');\n'
#                 w_func_body += ' '*cfg.indent + kw_str + return_type_base + '*& _temp_pr = _temp_p;\n'
#                 w_func_body += ' '*cfg.indent + 'return _temp_pr;\n'

#                 # Pythia8::Particle*  temp_p   = new Pythia8::Particle(front());
#                 # Pythia8::Particle*& temp_pr  = temp_p;
#                 # return temp_pr;

#             else:                
#                 w_func_body += ' '*cfg.indent + 'return new ' + return_type_base + '(' + orig_func_call + ');\n'
#         else:
#             w_func_body += ' '*cfg.indent + 'return ' + orig_func_call + ';\n'

#             # w_func_body += ' '*cfg.indent + 'return reinterpret_cast<' + w_return_type_base + '*'*n_pointers_w_return + '>(' + orig_func_call + ');\n'

#         # w_func_body += ' '*cfg.indent + return_type + ' ' + temp_res_name + ' = ' + func_name + '(' + ','.join(call_args) + ')' ';\n'
#         # if return_is_native:
#         #     w_func_body += ' '*cfg.indent + 'return ' + temp_res_name + '.abstractify()' + ';\n'
#         #     #
#         #     # FIXME: Take 'pointerness' of return type into account
#         #     #
#         # else:
#         #     w_func_body += ' '*cfg.indent + 'return ' + temp_res_name + ';\n'
#     w_func_body = '{\n' + w_func_body + '}'

#     return w_func_body

# # ======== END: constrWrapperBody ========



# ======== constrWrapperBody ========

def constrWrapperBody(return_type, func_name, args, return_is_loaded_class, keywords=[]):

    #
    # Input:
    # 
    #   - Return type of original function
    #   - List of dicts for original arguments
    #   - Name of original function
    #   - Boolean stating whether the orignal return type is native

    # Pointer and reference check
    pointerness, is_ref = utils.pointerAndRefCheck(return_type, byname=True)

    # Generate bracket for calling original function
    args_bracket_notypes = constrArgsBracket(args, include_arg_type=False, cast_to_original=True)

    w_func_body = ''
    w_func_body += '{\n'

    w_func_body += ' '*cfg.indent
    
    if return_type == 'void':
        w_func_body += func_name + args_bracket_notypes + ';\n'

    else:
        w_func_body += 'return '

        use_return_type = return_type
        if is_ref:
            use_return_type = return_type.rstrip('&')

        # The 'new SomeType'-statement should have one less '*' than the return type
        use_return_type.rstrip('*')

        if return_is_loaded_class:
            w_func_body += 'new ' + use_return_type + '(' + func_name + args_bracket_notypes + ');\n'
        else:
            w_func_body += func_name + args_bracket_notypes + ';\n'

    w_func_body += '}'

    return w_func_body

# ======== END: constrWrapperBody ========



# ======== ignoreFunction ========

def ignoreFunction(func_el, limit_pointerness=False):

    return_type_accepted = True
    arg_types_accepted   = True

    # Check function return type
    if 'returns' in func_el.keys():
        if not utils.isAcceptedType(func_el):
            return_type, return_kw, return_id = utils.findType(func_el)
            # warnings.warn('The function "%s" makes use of a non-accepted return type "%s" and will be ignored.' % (func_el.get('name'), return_type))
            print 'INFO: ' + 'The function "%s" makes use of a non-accepted return type "%s" and will be ignored.' % (func_el.get('name'), return_type)
            return_type_accepted = False

    # Check argument types
    args = getArgs(func_el)
    for arg_dict in args:
        arg_type_name = arg_dict['type']
        if not utils.isAcceptedType(arg_type_name, byname=True):
            # warnings.warn('The function "%s" makes use of a non-accepted argument type "%s" and will be ignored.' % (func_el.get('name'), arg_type_name))
            print 'INFO: ' + 'The function "%s" makes use of a non-accepted argument type "%s" and will be ignored.' % (func_el.get('name'), arg_type_name)
            arg_types_accepted = False
            break
        if limit_pointerness == True:
            if utils.isLoadedClass(arg_type_name, byname=True):
                if ('**' in arg_type_name) or ('*&' in arg_type_name):
                    # warnings.warn('The function "%s" makes use of a pointer-to-pointer or reference-to-pointer ("%s") for a loaded class. Such types cannot be handled safely by the BOSS wrapper system and thus this function will be ignored.' % (func_el.get('name'), arg_type_name))
                    print 'INFO: ' + 'The function "%s" makes use of a pointer-to-pointer or reference-to-pointer ("%s") for a loaded class. Such types cannot be handled safely by the BOSS wrapper system and thus this function will be ignored.' % (func_el.get('name'), arg_type_name)
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



# ======== usesLoadedType ========

def usesLoadedType(func_el):

    uses_loaded_type = False

    return_type_name, return_kw, return_id = utils.findType(func_el)
    return_is_loaded = utils.isLoadedClass( cfg.id_dict[return_id] )

    args = getArgs(func_el)
    is_arg_loaded = [arg_dict['loaded_class'] for arg_dict in args]

    if (return_is_loaded) or (True in is_arg_loaded):
        uses_loaded_type = True

    return uses_loaded_type

# ======== END: usesLoadedType ========