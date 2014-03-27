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

#
# Module-level globals
#


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
        print 'Current function: ' + func_name_full
        
        # NOT NEEDED WHEN PARSING ONLY LISTED FUNCTIONS:
        #
        # # Check if this function is native to the source code
        # if utils.isNative(func_el):
        #     if 'extern' in func_el.keys() and func_el.get('extern') == "1":
        #         print ' '*15 + '--> Function %s RECOGNIZED but IGNORED!' % func_name_full                
        #         continue
        #     else:
        #         print ' '*15 + '--> Function %s ACCEPTED!' % func_name_full
        # else:
        #     # print ' '*15 + '--> Function IGNORED!'
        #     continue


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


        # Function return type
        return_type, return_kw, return_id = utils.findType( cfg.id_dict[func_el.get('returns')] )
        return_el = cfg.id_dict[return_id]
        return_is_native = utils.isNative(return_el)


        # Function arguments (get list of dicts with argument info)
        args = funcutils.getArgs(func_el)


        # Construct wrapper function name
        w_func_name = funcutils.constrWrapperName(func_el)


        # Choose wrapper return type
        return_type_base = return_type.replace('*','').replace('&','')
        if (return_type == 'void') or (return_type.count('*') > 0):
            w_return_type = return_type
        else:
            w_return_type = return_type.replace(return_type_base, return_type_base+'*')
        # if return_is_native:
        #     w_return_type = cfg.abstr_class_prefix + return_type
        # else:
        #     w_return_type = return_type


        # Construct list of arguments for wrapper function
        w_args = funcutils.constrWrapperArgs(args)


        # Construct bracket with input arguments for wrapper function
        w_args_bracket = funcutils.constrArgsBracket(w_args)


        # Construct declaration line for wrapper function
        w_func_line = funcutils.constrDeclLine(w_return_type, w_func_name, w_args_bracket, keywords=return_kw)


        # Construct function body for wrapper function
        w_func_body = funcutils.constrWrapperBody(return_type, func_name, args, return_is_native)


        # Combine new code

        # Get info
        source_file_el   = cfg.id_dict[func_el.get('file')]
        source_file_name = source_file_el.get('name')

        # Prepare element in return_code_dict
        if source_file_name not in return_code_dict.keys():
            return_code_dict[source_file_name] = []

        # Define code string
        n_indents = 0
        add_code  = 2*'\n'

        # Construct the beginning of the namespaces
        for ns in namespaces:
            add_code += ' '*n_indents*cfg.indent + 'namespace ' + ns + '\n'
            add_code += ' '*n_indents*cfg.indent + '{' + '\n'
            n_indents += 1

        add_code += utils.addIndentation(w_func_line, n_indents*cfg.indent) + '\n'
        add_code += utils.addIndentation(w_func_body, n_indents*cfg.indent) + '\n'

        # - Construct the closing of the namespaces
        for ns in namespaces:
            n_indents -= 1
            add_code += ' '*n_indents*cfg.indent + '}' + '\n'
        
        add_code += '\n'


        # Register new code in return_code_dict
        insert_pos = -1   # end of file
        return_code_dict[source_file_name].append( (insert_pos, add_code) )


    #
    # Return result
    #

    return return_code_dict
