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



#
# Main function for parsing functions
#

def run():

    # Prepare returned dict
    return_code_dict = OrderedDict()

    #
    # Loop over all functions 
    #
    
    for func_name_full, func_el in gb.func_dict.items():

        # Print current function
        print
        print '~~~~ Current function: ' + func_name_full + ' ~~~~'
        print


        # Check if this function is accepted
        if funcutils.ignoreFunction(func_el):
            continue
       
        # Function name and namespace
        func_name = func_el.get('name')
        namespaces = func_name_full.split('<',1)[0].split('::')[:-1]
        has_namespace = bool(len(namespaces))


        # Check if this is a template function
        if '<' in func_name_full:
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
        new_source_fname = os.path.join(cfg.extra_output_dir, func_name.lower() + gb.code_suffix + cfg.source_extension)

        # Get include statements
        include_statements = []

        # - Generate include statements based on the types used in the function
        include_statements += utils.getIncludeStatements(func_el, convert_loaded_to='none', input_element='function', add_extra_include_path=True)
        include_statements += utils.getIncludeStatements(func_el, convert_loaded_to='wrapper_decl', input_element='function', add_extra_include_path=False, use_full_path=True)
        include_statements += utils.getIncludeStatements(func_el, convert_loaded_to='wrapper_def', input_element='function', add_extra_include_path=False, use_full_path=True)

        # - Then check if we have a header file for the function in question.
        #   If not, declare the original function as 'extern'
        file_el = gb.id_dict[func_el.get('file')]
        has_function_header = utils.isHeader(file_el)
        if has_function_header:
            header_full_path = file_el.get('name')
            header_base_name = os.path.basename(header_full_path)
            include_statements.append('#include "' + os.path.join(cfg.add_path_to_includes, header_base_name) + '"')

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
        wrapper_code = generateFunctionWrapperClassVersion(func_el, namespaces, n_overloads) 
        wrapper_code = utils.addIndentation(wrapper_code, len(namespaces)*cfg.indent)
        wrapper_code += '\n'

        # Prepare element in return_code_dict
        if new_source_fname not in return_code_dict.keys():
            return_code_dict[new_source_fname] = {'code_tuples':[], 'add_include_guard':False}

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


        # Register new code in return_code_dict
        insert_pos = -1   # end of file
        return_code_dict[new_source_fname]['code_tuples'].append( (insert_pos, new_code) )

        print
        print '~~~~ Function ' + func_name_full + ' done ~~~~'
        print
        print


    #
    # Return result
    #

    return return_code_dict




#
# Function for generating a source file containing wrapper functions (that make use of the wrapper classes)
#
def generateFunctionWrapperClassVersion(func_el, namespaces, n_overloads):

    new_code = ''


    #
    # Get info on function
    #

    # Identify arguments, translate argument type of loaded classes
    # and construct the argument bracket
    args = funcutils.getArgs(func_el)
    w_args = funcutils.constrWrapperArgs(args, add_ref=True)

    # Identify return type 
    return_type_dict = utils.findType( gb.id_dict[func_el.get('returns')] )
    return_el     = return_type_dict['el']
    pointerness   = return_type_dict['pointerness']
    is_ref        = return_type_dict['is_reference']
    return_kw     = return_type_dict['cv_qualifiers']
    
    return_kw_str = ' '.join(return_kw) + ' '*bool(len(return_kw))

    return_type   = return_type_dict['name'] + '*'*pointerness + '&'*is_ref


    #
    # Wrapper function
    #

    wrapper_code = '// Wrapper function(s)\n'

    # Function name
    func_name = func_el.get('name') + gb.code_suffix

    # Check constness
    if ('const' in func_el.keys()) and (func_el.get('const')=='1'):
        is_const = True
    else:
        is_const = False

    # One function for each set of default arguments
    for remove_n_args in range(n_overloads+1):

        if remove_n_args == 0:
            use_args   = args
            use_w_args = w_args
        else:
            use_args   = args[:-remove_n_args]
            use_w_args = w_args[:-remove_n_args]

        args_bracket_wrapper = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=True, include_namespace=True, use_wrapper_class=True)
        args_bracket_wrapper_notypes = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=False, cast_to_original=True, wrapper_to_pointer=True)

        # Name of function to call
        call_func_name = func_el.get('name')

        # Convert return type if loaded class
        if return_is_loaded_class:
            wrapper_return_type = classutils.toWrapperType(return_type, remove_reference=True)
        else:
            wrapper_return_type = return_type

        # Write declaration line
        wrapper_code += return_kw_str + wrapper_return_type + ' ' + func_name + args_bracket_wrapper + is_const*' const' + '\n'

        # Write function body
        indent = ' '*cfg.indent
        wrapper_code += '{\n'

        if return_type == 'void':
            wrapper_code += indent
        else:
            wrapper_code += indent + 'return '

        wrapper_code += call_func_name + args_bracket_wrapper_notypes + ';\n'

        wrapper_code += '}\n'
        wrapper_code += '\n'

    new_code += wrapper_code

    return new_code
