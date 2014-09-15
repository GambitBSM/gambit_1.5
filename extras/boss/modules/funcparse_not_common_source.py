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
import modules.utils as utils
import modules.funcutils as funcutils
import modules.classutils as classutils

#
# Module-level globals
#

wrapper_source_include_statements = []

#
# Main function for parsing functions
#

def run():

    # Prepare returned dict
    return_code_dict = OrderedDict()

    #
    # Loop over all functions 
    #
    
    for func_name_full, func_el in cfg.func_dict.items():

        # Check if this function is accepted
        if funcutils.ignoreFunction(func_el):
            continue

        # Print current function
        print
        print '~~~~ Current function: ' + func_name_full + ' ~~~~'
        print
        

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
        # Generate extra source file with overloaded and 'abstract' versions
        #

        # New source file name
        new_source_fname = os.path.join(cfg.extra_output_dir, func_name.lower() + cfg.code_suffix + cfg.source_extension)

        # Get include statements
        include_statements = []

        # - Generate include statements based on the types used in the function
        include_statements += utils.getIncludeStatements(func_el, convert_loaded_to='none', input_element='function', add_extra_include_path=True)

        # - Then check if we have a header file for the function in question.
        #   If not, declare the original function as 'extern'
        file_el = cfg.id_dict[func_el.get('file')]
        has_function_header = utils.isHeader(file_el)
        if has_function_header:
            header_full_path = file_el.get('name')
            header_base_name = os.path.basename(header_full_path)
            include_statements.append('#include "' + os.path.join(cfg.add_path_to_includes, header_base_name) + '"')

        include_statements = list( set(include_statements) )
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


        abstr_code = ''
        for remove_n_args in range(n_overloads+1):

            # Generate code for abstract version
            abstr_code += classutils.constrWrapperFunction(func_el, indent=cfg.indent, n_indents=len(namespaces), remove_n_args=remove_n_args)
            abstr_code += '\n'

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

        # - Add code for 'abstract' version
        new_code += abstr_code

        # - Construct the closing of the namespaces
        new_code += utils.constrNamespace(namespaces, 'close')
        
        new_code += '\n'


        # Register new code in return_code_dict
        insert_pos = -1   # end of file
        return_code_dict[new_source_fname]['code_tuples'].append( (insert_pos, new_code) )


        # Generate a version of the function that uses the wrapper classes. Add this to a common source file
        short_wrapper_function_fname = cfg.wrapper_header_prefix + 'functions.cpp'
        wrapper_function_fname = os.path.join(cfg.extra_output_dir, short_wrapper_function_fname)
        wrapper_function_file_content = generateWrapperFunctionSource(func_el, namespaces, n_overloads)

        if wrapper_function_fname not in return_code_dict.keys():
            return_code_dict[wrapper_function_fname] = {'code_tuples':[], 'add_include_guard':False}
        return_code_dict[wrapper_function_fname]['code_tuples'].append( (-1, wrapper_function_file_content) )


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
def generateWrapperFunctionSource(func_el, namespaces, n_overloads):

    new_code = ''

    #
    # Generate include statements
    #
    
    include_statements = []
    include_statements.append('#include <cstddef>')
    include_statements.append('#include "nullptr_check' + cfg.header_extension + '"')
    
    include_statements += utils.getIncludeStatements(func_el, convert_loaded_to='wrapper', input_element='function')
    
    include_statements_code = ''

    for inc_statement in include_statements:
        if inc_statement not in wrapper_source_include_statements:
            include_statements_code += inc_statement + '\n'
            wrapper_source_include_statements.append(inc_statement)

    new_code += include_statements_code


    #
    # Get info on function
    #

    # Identify arguments, translate argument type of loaded classes
    # and construct the argument bracket
    args = funcutils.getArgs(func_el)
    w_args = funcutils.constrWrapperArgs(args, add_ref=True)

    # Identify return type 
    return_type, return_kw, return_id = utils.findType( cfg.id_dict[func_el.get('returns')] )
    return_el = cfg.id_dict[return_id]

    return_is_loaded_class = utils.isLoadedClass(return_el)
    pointerness, is_ref = utils.pointerAndRefCheck(return_type, byname=True)

    # Return keywords
    return_kw_str = ' '.join(return_kw)        
    return_kw_str += ' '*bool(len(return_kw))


    #
    # Create function pointers to be filled by dynamic loading
    #

    # Choose function pointer return type
    if return_is_loaded_class:
        if pointerness == 0:
            func_ptr_return_type = classutils.toAbstractType(return_type, include_namespace=True, add_pointer=True, remove_reference=True)
        else:
            func_ptr_return_type = classutils.toAbstractType(return_type, include_namespace=True)
    else:
        func_ptr_return_type = return_type


    # One function pointer for each set of default arguments
    func_ptr_code = ''
    for remove_n_args in range(n_overloads+1):

        if remove_n_args == 0:
            use_w_args = w_args
        else:
            use_w_args = w_args[:-remove_n_args]

        args_bracket = funcutils.constrArgsBracket(use_w_args, include_arg_name=False, include_arg_type=True, include_namespace=True)

        # Function pointer name
        func_ptr_name = func_el.get('name') + '_ptr'
        if remove_n_args > 0:
            func_ptr_name += '_overload_' + str(remove_n_args)

        # Construct function pointer code
        func_ptr_code += func_ptr_return_type + ' ' + return_kw_str + '(*' + func_ptr_name + ')' + args_bracket + ' = NULL;\n'

    if func_ptr_code != '':
        new_code += '\n'
        new_code += "// Function pointer(s) to be filled by dynamic loading\n"
        new_code += func_ptr_code + '\n'


    #
    # Wrapper function
    #

    wrapper_code = '// Wrapper function(s)\n'

    # Function name
    func_name = func_el.get('name')

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
        args_bracket_wrapper_notypes  = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=False, wrapper_to_pointer=True)
        args_bracket_func_ptr_nonames = funcutils.constrArgsBracket(use_w_args, include_arg_name=False, include_arg_type=True, include_namespace=True)

        # Name of function pointer to call
        call_func_ptr_name = func_el.get('name') + '_ptr'
        if remove_n_args > 0:
            call_func_ptr_name += '_overload_' + str(remove_n_args)

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


        # wrapper_code += 'nullptr_check(' + call_func_ptr_name + ')' + args_bracket_notypes + ';\n'

        wrapper_code += 'nullptr_check< '+ func_ptr_return_type + ' ' + return_kw_str + '(*)' + args_bracket_func_ptr_nonames + ' >(' + call_func_ptr_name + ')' + args_bracket_wrapper_notypes + ';\n'

        wrapper_code += '}\n'
        wrapper_code += '\n'

    new_code += wrapper_code


    new_code += 2*'\n'
    return new_code



        # # Function return type
        # return_type, return_kw, return_id = utils.findType( cfg.id_dict[func_el.get('returns')] )
        # return_el = cfg.id_dict[return_id]
        # return_is_native = utils.isNative(return_el)


        # # Function arguments (get list of dicts with argument info)
        # args = funcutils.getArgs(func_el)


        # # Construct wrapper function name
        # w_func_name = funcutils.constrWrapperName(func_el)


        # # Choose wrapper return type
        # return_type_base = return_type.replace('*','').replace('&','')
        # if (return_type == 'void') or (return_type.count('*') > 0):
        #     w_return_type = return_type
        # else:
        #     w_return_type = return_type.replace(return_type_base, return_type_base+'*')
        # # if return_is_native:
        # #     w_return_type = cfg.abstr_class_prefix + return_type
        # # else:
        # #     w_return_type = return_type


        # # Construct list of arguments for wrapper function
        # w_args = funcutils.constrWrapperArgs(args)

        # # Construct bracket with input arguments for wrapper function
        # w_args_bracket = funcutils.constrArgsBracket(w_args)

        # # Construct declaration line for wrapper function
        # w_func_line = funcutils.constrDeclLine(w_return_type, w_func_name, w_args_bracket, keywords=return_kw)


        # # Construct function body for wrapper function
        # w_func_body = funcutils.constrWrapperBody(return_type, func_name, args, return_is_native)


        # Combine new code

        # Get info
        # source_file_el   = cfg.id_dict[func_el.get('file')]
        # source_file_name = source_file_el.get('name')
