###############################
#                             #
#  Function parsing for BOSS  #
#                             #
###############################

import xml.etree.ElementTree as ET
from collections import OrderedDict
import os
import copy

import modules.cfg as cfg
import modules.gb as gb
import modules.utils as utils
import modules.funcutils as funcutils
import modules.classutils as classutils
import modules.infomsg as infomsg


#
# Module-level globals
#



# ====== run ========

# Main function for parsing functions

def run():


    #
    # Loop over all functions 
    #
    
    for func_name_full, func_el in gb.func_dict.items():

        # Clear all info messages
        infomsg.clearInfoMessages()

        # Generate dict with different variations of the function name
        func_name = funcutils.getFunctionNameDict(func_el)


        # Print current function
        print
        print '  Function: ' + func_name['long_templ']
        print '  ----------' + '-'*len(func_name['long_templ'])


        # Check if this function is accepted
        if funcutils.ignoreFunction(func_el):
            continue
       
        # Function namespace
        # namespaces = func_name['long_templ'].split('<',1)[0].split('::')[:-1]
        namespaces = utils.getNamespaces(func_el)
        has_namespace = bool(len(namespaces))


        # Check if this is a template function
        if '<' in func_name['long_templ']:
            is_template = True
        else:
            is_template = False


        # If template function, figure out template variables
        if is_template == True:
            template_bracket, template_types = utils.getTemplateBracket(func_el)
            spec_template_types = utils.getSpecTemplateTypes(func_el)
            print 'TEMPLATE: ', template_bracket, template_types, spec_template_types


        #
        # Generate extra source file with overloaded and wrapper class versions
        #

        # New source file name
        new_source_fname = os.path.join(cfg.extra_output_dir, cfg.function_files_prefix + func_name['short'].lower() + gb.code_suffix + cfg.source_extension)

        # Get include statements
        include_statements = []

        # - Generate include statements based on the types used in the function
        include_statements += utils.getIncludeStatements(func_el, convert_loaded_to='none', input_element='function')
        include_statements += utils.getIncludeStatements(func_el, convert_loaded_to='wrapper_decl', input_element='function', use_full_path=True)
        include_statements += utils.getIncludeStatements(func_el, convert_loaded_to='wrapper_def', input_element='function', use_full_path=True)

        # - Then check if we have a header file for the function in question.
        #   If not, declare the original function as 'extern'
        file_el = gb.id_dict[func_el.get('file')]
        has_function_header = utils.isHeader(file_el)
        if has_function_header:
            header_full_path = file_el.get('name')
            use_path = utils.shortenHeaderPath(header_full_path)
            include_statements.append( '#include "' + use_path + '"')

        # Add include statement for gambit/Backends/function_return_utils.hpp
        include_statements.append( '#include "' + os.path.join(gb.gambit_backend_incl_dir,'function_return_utils.hpp') + '"')

        include_statements = list( OrderedDict.fromkeys(include_statements) )
        include_statements_code = '\n'.join(include_statements) + 2*'\n'


        # If no header file is found for the original function, generate 'extern' declaration
        extern_declaration = ''
        if not has_function_header:
            extern_declaration += funcutils.constrExternFuncDecl(func_el)
            extern_declaration += '\n'


        # If we have access to the function header, we can implement one overloaded versions
        # to account for default value arguments.
        if has_function_header:
            n_overloads = funcutils.numberOfDefaultArgs(func_el)
        else:
            n_overloads = 0


        # Generate code for wrapper class version
        # wrapper_code = classutils.constrWrapperFunction(func_el, indent=cfg.indent, n_indents=0, 
        #                                                         remove_n_args=0, include_full_namespace=False)
        # wrapper_code = generateFunctionWrapperClassVersion(func_el, func_name, namespaces, n_overloads) 
        wrapper_code = generateFunctionWrapperClassVersion(func_el, namespaces, n_overloads) 
        wrapper_code = utils.addIndentation(wrapper_code, len(namespaces)*cfg.indent)
        wrapper_code += '\n'


        # Prepare element in gb.new_code
        if new_source_fname not in gb.new_code.keys():
            gb.new_code[new_source_fname] = {'code_tuples':[], 'add_include_guard':False}

        # Define code string
        n_indents = len(namespaces)
        new_code  = 2*'\n'

        # - Add include statements
        new_code += include_statements_code

        # - Add extern function declaration
        new_code += extern_declaration

        # - Construct the beginning of the namespaces
        new_code += utils.constrNamespace(namespaces, 'open')

        # - Add code for 'wrapper' version
        new_code += wrapper_code

        # - Construct the closing of the namespaces
        new_code += utils.constrNamespace(namespaces, 'close')
        
        new_code += '\n'


        print 
        print '****** CODE ******'
        print 
        print new_code
        print 
        print '****** CODE END ******'
        print


        # Register new code in return_code_dict
        insert_pos = -1   # end of file
        # return_code_dict[new_source_fname]['code_tuples'].append( (insert_pos, new_code) )
        gb.new_code[new_source_fname]['code_tuples'].append( (insert_pos, new_code) )


        #
        # Keep track of functions done
        #
        gb.functions_done.append(func_name)

        # Increase class counter
        gb.n_functions_done += 1

        print

    #
    # End loop over functions
    #

# ====== END: run ========




# ====== generateFunctionWrapperClassVersion ========

# Function for generating a source file containing wrapper
# functions that make use of the wrapper classes.

def generateFunctionWrapperClassVersion(func_el, namespaces, n_overloads):

    new_code = ''

    # Check if this function makes use of any loaded types
    uses_loaded_type = funcutils.usesLoadedType(func_el)

    # Function name
    func_name = func_el.get('name')

    # Name of wrapper function
    wr_func_name = func_name + gb.code_suffix


    # Determine return type
    return_type_dict = utils.findType(func_el)
    return_el = return_type_dict['el']
    pointerness    = return_type_dict['pointerness']
    is_ref         = return_type_dict['is_reference']
    return_type_kw = return_type_dict['cv_qualifiers']
    
    return_kw_str  = ' '.join(return_type_kw) + ' '*bool(len(return_type_kw))
    
    return_is_loaded    = utils.isLoadedClass(return_el)

    return_type   = return_type_dict['name'] + '*'*pointerness + '&'*is_ref


    # If return-by-value, then a const qualifier on the return value is meaningless
    # (will result in a compiler warning)
    if (pointerness == 0) and (is_ref == False):
        if 'const' in return_type_kw:
            return_type_kw.remove('const')


    # Arguments
    args = funcutils.getArgs(func_el)


    # One function for each set of default arguments
    n_overloads = funcutils.numberOfDefaultArgs(func_el)
    for remove_n_args in range(n_overloads+1):

        # Check that the function is acceptable
        if funcutils.ignoreFunction(func_el, limit_pointerness=True, remove_n_args=remove_n_args):
            continue

        if remove_n_args == 0:
            use_args = args
        else:
            use_args = args[:-remove_n_args]

        # Argument bracket
        args_bracket = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=True, include_namespace=True, use_wrapper_class=True, use_wrapper_base_class=False)                

        # Name of original function to call
        call_func_name = func_name

        # Convert return type if loaded class
        if utils.isLoadedClass(return_el):
            wrapper_return_type = classutils.toWrapperType(return_type, remove_reference=True)
        else:
            wrapper_return_type = return_type

        # Write declaration line
        new_code += return_kw_str + wrapper_return_type + ' ' + wr_func_name + args_bracket + '\n'

        # Write function body
        indent = ' '*cfg.indent
        new_code += '{\n'

        if return_type == 'void':
            new_code += indent
        else:
            new_code += indent + 'return '

        # args_bracket_notypes = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=False, wrapper_to_pointer=True)
        args_bracket_notypes = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=False, cast_to_original=True, wrapper_to_pointer=True)

        if return_is_loaded: 

            abs_return_type_simple = classutils.toAbstractType(return_type, include_namespace=True, remove_reference=True, remove_pointers=True)
            wrapper_return_type_simple = wrapper_return_type.replace('*','').replace('&','')

            if is_ref:  # Return-by-reference
                new_code += 'reference_returner< ' + wrapper_return_type_simple + ', ' + abs_return_type_simple +  ' >( ' + call_func_name + args_bracket_notypes + ' );\n'

            elif (not is_ref) and (pointerness > 0):  # Return-by-pointer
                new_code += 'pointer_returner< ' + wrapper_return_type_simple + ', ' + abs_return_type_simple +  ' >( ' + call_func_name + args_bracket_notypes + ' );\n'
            
            else:  # Return-by-value
                new_code += wrapper_return_type + '( ' + call_func_name + args_bracket_notypes + ' );\n'
        
        else:                
            new_code += call_func_name + args_bracket_notypes + ';\n'

        new_code += '}\n'
        new_code += '\n'


    return new_code

# ====== END: generateFunctionWrapperClassVersion ========



# # ====== generateFunctionWrapperClassVersion ========

# # Function for generating a source file containing wrapper
# # functions that make use of the wrapper classes.

# def generateFunctionWrapperClassVersion(func_el, func_name, namespaces, n_overloads):

#     new_code = ''


#     #
#     # Get info on function
#     #


#     # Identify arguments, translate argument type of loaded classes
#     # and construct the argument bracket
#     args = funcutils.getArgs(func_el)
#     w_args = funcutils.constrWrapperArgs(args, add_ref=True)

#     # Identify return type 
#     return_type_dict = utils.findType( gb.id_dict[func_el.get('returns')] )
#     return_el     = return_type_dict['el']
#     pointerness   = return_type_dict['pointerness']
#     is_ref        = return_type_dict['is_reference']
#     return_kw     = return_type_dict['cv_qualifiers']
    
#     return_kw_str = ' '.join(return_kw) + ' '*bool(len(return_kw))

#     return_type   = return_type_dict['name'] + '*'*pointerness + '&'*is_ref


#     #
#     # Wrapper function
#     #

#     wrapper_code = '// Wrapper function(s)\n'

#     # Wrapper function name
#     wr_func_name = func_name['short'] + gb.code_suffix

#     # Check constness
#     if ('const' in func_el.keys()) and (func_el.get('const')=='1'):
#         is_const = True
#     else:
#         is_const = False

#     # One function for each set of default arguments
#     for remove_n_args in range(n_overloads+1):

#         if remove_n_args == 0:
#             use_args   = args
#             use_w_args = w_args
#         else:
#             use_args   = args[:-remove_n_args]
#             use_w_args = w_args[:-remove_n_args]

#         args_bracket_wrapper = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=True, include_namespace=True, use_wrapper_class=True)
#         args_bracket_wrapper_notypes = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=False, cast_to_original=True, wrapper_to_pointer=True)

#         # Name of function to call
#         call_func_name = func_name['short']

#         # Convert return type if loaded class
#         if utils.isLoadedClass(return_el):
#             wrapper_return_type = classutils.toWrapperType(return_type, remove_reference=True)
#         else:
#             wrapper_return_type = return_type

#         # Write declaration line
#         wrapper_code += return_kw_str + wrapper_return_type + ' ' + wr_func_name + args_bracket_wrapper + is_const*' const' + '\n'

#         # Write function body
#         indent = ' '*cfg.indent
#         wrapper_code += '{\n'

#         if return_type == 'void':
#             wrapper_code += indent
#         else:
#             wrapper_code += indent + 'return '

#         wrapper_code += call_func_name + args_bracket_wrapper_notypes + ';\n'

#         wrapper_code += '}\n'
#         wrapper_code += '\n'

#     new_code += wrapper_code

#     return new_code

# # ====== END: generateFunctionWrapperClassVersion ========

