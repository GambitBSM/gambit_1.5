############################
#                          #
#  Class parsing for BOSS  #
#                          #
############################

import xml.etree.ElementTree as ET
from collections import OrderedDict
import os
import warnings

import modules.cfg as cfg
import modules.gb as gb
import modules.utils as utils
import modules.classutils as classutils
import modules.funcutils as funcutils

#
# Module-level globals
#

classes_done          = []
template_done         = []
templ_spec_done       = []
added_parent          = []

#
# Main function for parsing classes
#

def run():

    # Prepare returned dict
    return_code_dict = OrderedDict()

    # Prepare dict to keep track of include statements
    includes = OrderedDict()

    #
    # Loop over all classes 
    #
    
    for class_name_long, class_el in gb.class_dict.items():

        # Print current class
        print
        print '~~~~ Current class: ' + class_name_long + ' ~~~~'
        print

        # Check if we've done this class already
        if class_name_long in classes_done:
            print ' '*15 + '--> Class already done'
            continue

        # Generate dicts with different variations of the class name
        class_name       = classutils.getClassNameDict(class_el)
        abstr_class_name = classutils.getClassNameDict(class_el, abstract=True)

        # Check if this is a template class
        if '<' in class_name['long_templ']:
            is_template = True
        else:
            is_template = False

        # Check if there are any pure virtual methods in this class.
        # If yes, print a warning and skip this class. (We cannot load such classes.)
        check_member_elements = utils.getMemberElements(class_el)        
        pure_virtual_members = []
        for mem_el in check_member_elements:
            if mem_el.tag in ['Constructor', 'Destructor', 'Method', 'OperatorMethod']:
                if ('pure_virtual' in mem_el.keys()) and (mem_el.get('pure_virtual')=='1'):
                    pure_virtual_members.append(mem_el.get('name'))
        if len(pure_virtual_members) > 0:
            print 'INFO: ' + 'The following member functions in class "%s" are pure virtual: %s' % (class_name['long_templ'], ', '.join(pure_virtual_members))
            print 'INFO: ' + 'Class "%s" cannot be loaded and will be ignored.' % (class_name['long_templ'])
            gb.class_dict.pop(class_name_long)
            continue

        # Make list of all types in use in this class
        all_types_in_class = utils.getAllTypesInClass(class_el, include_parents=True)

        # Set a bunch of generally useful variables 
        original_class_file_el        = gb.id_dict[class_el.get('file')]
        original_class_file_name      = original_class_file_el.get('name')
        original_class_file_name_base = os.path.basename(original_class_file_name)
        original_class_file_dir       = os.path.split(original_class_file_name)[0]
        extras_src_file_name          = os.path.join(cfg.extra_output_dir, class_name['short'] + '_extras' + gb.code_suffix + cfg.source_extension)

        short_abstr_class_fname = gb.new_header_files[class_name['long']]['abstract']
        abstr_class_fname       = os.path.join(cfg.extra_output_dir, short_abstr_class_fname)

        namespaces    = class_name['long'].split('::')[:-1]
        has_namespace = bool(len(namespaces))

        has_copy_constructor, copy_constructor_id         = classutils.checkCopyConstructor(class_el, return_id=True)
        has_assignment_operator, assignment_is_artificial = classutils.checkAssignmentOperator(class_el)
        
        if has_assignment_operator and assignment_is_artificial:
            construct_assignment_operator = True
        else:
            construct_assignment_operator = False


        # Prepare entries in return_code_dict and includes
        if abstr_class_fname not in return_code_dict.keys():
            return_code_dict[abstr_class_fname] = {'code_tuples':[], 'add_include_guard':True}
        if original_class_file_name not in return_code_dict.keys():
            return_code_dict[original_class_file_name] = {'code_tuples':[], 'add_include_guard':False}
        if original_class_file_name not in includes.keys():
            includes[original_class_file_name] = []
        if extras_src_file_name not in return_code_dict.keys():
            return_code_dict[extras_src_file_name] = {'code_tuples':[], 'add_include_guard':False}


        # Treat the first specialization of a template class differently
        if is_template and class_name['long'] not in template_done:
            template_bracket, template_types = utils.getTemplateBracket(class_el)
            
            empty_templ_class_decl = ''
            empty_templ_class_decl += classutils.constrEmptyTemplClassDecl(abstr_class_name['short'], namespaces, template_bracket, indent=cfg.indent)
            empty_templ_class_decl += classutils.constrTemplForwDecl(class_name['short'], namespaces, template_bracket, indent=cfg.indent)

            return_code_dict[abstr_class_fname]['code_tuples'].append( (0, empty_templ_class_decl) )


        # Get template arguments for specialization, 
        # and check that they are acceptable
        if is_template and class_name['long'] not in templ_spec_done:
            spec_template_types = utils.getSpecTemplateTypes(class_el)
            for template_type in spec_template_types:
                if (template_type not in gb.accepted_types):
                    raise Exception("The template specialization type '" + template_type + "' for class " + class_name['long'] + " is not among accepted types.")


        #
        # Construct code for the abstract class header file and register it
        #
        
        class_decl = ''

        # - Add include statements
        include_statements  = []
        include_statements  = ['#include "abstractbase.hpp"']
        include_statements += ['#include "' + gb.frwd_decls_abs_fname + cfg.header_extension + '"']
        include_statements += ['#include "' + gb.frwd_decls_wrp_fname + cfg.header_extension + '"']
        # include_statements += classutils.getIncludeStatements(all_types_in_class, 'abstract', exclude_types=[class_name])
        include_statements += utils.getIncludeStatements(class_el, convert_loaded_to='abstract', exclude_types=[class_name])
        include_statements_code = '\n'.join(include_statements) + 2*'\n'
        class_decl += include_statements_code


        # - Add forward declarations at the top of the header file  
        #
        #   FIXME: For now, this includes forward declarations of *all* loaded native classes.
        #          Should limit this to only the forward declarations needed for the given abstract class.
        #   FIXME: Anything depending on the forward declaration of the original class should go into a 
        #          separate header (abstract_SomeClass_extra.hpp).
        # 
        # class_decl += '// Forward declarations:\n'
        # class_decl += utils.constrForwardDecls()
        # class_decl += '\n'


        # - Add the the code for the abstract class
        if (is_template == True) and (class_name['long'] in templ_spec_done):
            pass
        elif (is_template == True) and (class_name['long'] not in templ_spec_done):
            class_decl += classutils.constrAbstractClassDecl(class_el, class_name['short'], abstr_class_name['short'], namespaces, 
                                                             indent=cfg.indent, template_types=spec_template_types, 
                                                             has_copy_constructor=has_copy_constructor, construct_assignment_operator=construct_assignment_operator)
            class_decl += '\n'
        else:
            class_decl += classutils.constrAbstractClassDecl(class_el, class_name['short'], abstr_class_name['short'], namespaces, indent=cfg.indent, 
                                                             has_copy_constructor=has_copy_constructor, construct_assignment_operator=construct_assignment_operator)
            class_decl += '\n'

        # - Register code
        return_code_dict[abstr_class_fname]['code_tuples'].append( (-1, class_decl) )


        #
        # Add abstract class to inheritance list of original class
        #

        line_number = int(class_el.get('line'))

        # Get file content, but replace all comments with whitespace
        f = open(original_class_file_name, 'r')
        file_content = f.read()
        f.close()
        file_content_nocomments = utils.removeComments(file_content, insert_blanks=True)

        # Find index of the \n in line number line_number
        newline_pos = utils.findNewLinePos(file_content_nocomments, line_number)

        # First, find position of class name
        search_limit = newline_pos
        while search_limit > -1:
            pos = file_content_nocomments[:search_limit].rfind(class_name['short'])
            pre_char  = file_content_nocomments[pos-1]
            post_char = file_content_nocomments[pos+len(class_name['short'])]
            if (pre_char in [' ','\n','\t']) and (post_char in [' ', ':', '\n', '<', '{']):
                break
            else:
                search_limit = pos

        class_name_pos = pos


        # Special preparations for template classes:
        if is_template:
        
            # - Determine whether this is the source for the general template 
            #   or for a specialization (look for '<' after class name)
            temp_pos = class_name_pos + len(class_name['short'])
            while True:
                next_char = file_content_nocomments[temp_pos]
                if next_char not in [' ', '\t', '\n']:
                    break
                else:
                    temp_pos += 1
            if next_char == '<':
                src_is_specialization = True
            else:
                src_is_specialization = False

            # - Prepare the template bracket string
            if src_is_specialization:
                add_template_bracket = '<' + ','.join(spec_template_types) + '>'
            else:
                add_template_bracket = '<' + ','.join(template_types) + '>'


        # If no previous parent classes:
        if (class_el.get('bases') == "") and (class_name['long'] not in added_parent):

            # - Calculate insert position
            insert_pos = class_name_pos + len(class_name['short'])
            if is_template and src_is_specialization:
                insert_pos += len(add_template_bracket)

            # - Generate code
            add_code = ' : public virtual ' + abstr_class_name['short']
            if is_template == True:
                add_code += add_template_bracket

        # If there are previous parent classes
        else:

            # - Get colon position
            if is_template and src_is_specialization:
                temp_pos = class_name_pos + len(class_name['short']) + len(add_template_bracket)
            else:
                temp_pos = class_name_pos + len(class_name['short'])
            colon_pos = temp_pos + file_content_nocomments[temp_pos:newline_pos].find(':')

            # - Calculate insert position
            insert_pos = colon_pos + 1

            # - Generate code
            add_code = ' public virtual ' + abstr_class_name['short']
            if is_template == True:
                add_code += add_template_bracket
            add_code += ','

        # - Register new code
        return_code_dict[original_class_file_name]['code_tuples'].append( (insert_pos, add_code) )

        # - Update added_parent dict
        added_parent.append(class_name['long'])


        # Generate code for #include statement in orginal header/source file 
        # include_line = '#include "' + os.path.join(cfg.add_path_to_includes, short_abstr_class_fname) + '"'
        # include_line = '#include "' + short_abstr_class_fname + '"'
        include_line = '#include "' + os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, short_abstr_class_fname ) + '"'
        if include_line in includes[original_class_file_name]:
            pass
        else:
            # - Find position
            if is_template == True:
                insert_pos = file_content_nocomments[:class_name_pos].rfind('template')
            else:
                insert_pos = max(file_content_nocomments[:class_name_pos].rfind('class'), file_content_nocomments[:class_name_pos].rfind('struct'))
            # - Adjust for the indentation
            use_indent = ''
            while insert_pos > 0:
                char = file_content[insert_pos-1]
                if char in [' ','\t']:
                    use_indent += char
                    insert_pos -= 1
                else:
                    break

            # - Construct code
            include_code = ''

            include_code += use_indent
            for ns in namespaces:
                include_code += '} '
            include_code += '\n'*has_namespace
            include_code += use_indent + include_line + '\n'
            include_code += use_indent + '#include "abstracts_typedefs.hpp"\n'
            include_code += use_indent
            for ns in namespaces:
                include_code += 'namespace ' + ns + ' { '
            include_code += '\n'*has_namespace

            # - Register code
            return_code_dict[original_class_file_name]['code_tuples'].append( (insert_pos, include_code) )

            # - Register include line
            includes[original_class_file_name].append(include_line)


        # Generate wrappers for all member functions that make use of native types and/or default arguments,
        # and reference returning functions for all public member variables that have an accepted type.
        # Put declarations in the original class header and implementations in a separate source file.
        
        # - Create lists of all 'non-artificial' members of the class
        member_methods   = []
        member_variables = []
        member_operators = []
        if 'members' in class_el.keys():
            for mem_id in class_el.get('members').split():
                el = gb.id_dict[mem_id]
                if not 'artificial' in el.keys():
                    if (el.tag == 'Method') and (not funcutils.ignoreFunction(el)):
                        member_methods.append(el)
                    elif (el.tag == 'OperatorMethod') and (not funcutils.ignoreFunction(el)):
                        if funcutils.usesNativeType(el):
                            member_operators.append(el)
                    elif (el.tag in ('Field', 'Variable')) and (el.get('access') == 'public'):
                        if utils.isAcceptedType(el):
                            member_variables.append(el)

                    # append_method = False
                    # return_type_name, return_kw, return_id = utils.findType(el)
                    # return_is_native = utils.isNative( gb.id_dict[return_id] )
                    # args = funcutils.getArgs(el)
                    # is_arg_native = [arg_dict['native'] for arg_dict in args]
                    # if return_is_native or True in is_arg_native:
                    #     member_methods.append(el)

        # - Determine insert position
        rel_pos_start, rel_pos_end = utils.getBracketPositions(file_content_nocomments[class_name_pos:], delims=['{','}'])
        class_body_start = class_name_pos + rel_pos_start
        class_body_end   = class_name_pos + rel_pos_end
        insert_pos = class_body_end

        # - Generate code for wrapper functions for each each member function
        #   A declaration goes into the original class header, 
        #   while implementations are put in a new source file.

        declaration_code     = '\n'
        implementation_code  = '\n'
        current_access = None
        for method_el in member_methods:

            # We need to generate as many overloaded versions as there are arguments with default values
            n_overloads = funcutils.numberOfDefaultArgs(method_el)
            
            # Check for native types
            uses_native_type = funcutils.usesNativeType(method_el)

            # If no native types are used and no arguments have default values, 
            # we don't need a wrapper
            if (not uses_native_type) and (n_overloads == 0):
                continue

            # Generate wrapper code
            for remove_n_args in range(n_overloads+1):
                
                if (remove_n_args==0) and (not uses_native_type):
                    continue

                # The declaration is put inside the original class
                method_access = method_el.get('access')
                if method_access != current_access:
                    declaration_code += ' '*(len(namespaces)+1)*cfg.indent + method_access +':\n'
                    current_access = method_access
                declaration_code += classutils.constrWrapperFunction(method_el, indent=cfg.indent, n_indents=len(namespaces)+2, 
                                                                 remove_n_args=remove_n_args, only_declaration=True)
                declaration_code += '\n'

                
                # The implementation goes into a new source file
                implementation_code += classutils.constrWrapperFunction(method_el, indent=cfg.indent, n_indents=0, 
                                                                        remove_n_args=remove_n_args, include_full_namespace=True)
                implementation_code += 2*'\n'

        # - Register code
        return_code_dict[original_class_file_name]['code_tuples'].append( (insert_pos, declaration_code) )            
        return_code_dict[extras_src_file_name]['code_tuples'].append( (-1, implementation_code) )            


        # - Generate code for each member operator
        operator_declaration_code    = '\n'
        operator_implementation_code = '\n'
        for operator_el in member_operators:
            operator_access = operator_el.get('access')
            if operator_access != current_access:
                operator_declaration_code += ' '*(len(namespaces)+1)*cfg.indent + operator_access +':\n'
                current_access = operator_access

            # If default arguments are used, we need several overloads
            n_overloads = funcutils.numberOfDefaultArgs(operator_el)
            for remove_n_args in range(n_overloads+1):

                # Put declaration in original class
                operator_declaration_code += classutils.constrWrapperFunction(operator_el, indent=cfg.indent, n_indents=len(namespaces)+2, 
                                                                              remove_n_args=remove_n_args, only_declaration=True)
                operator_declaration_code += '\n'


                # Put implementation in a new source file
                operator_implementation_code += classutils.constrWrapperFunction(operator_el, indent=cfg.indent, n_indents=0, 
                                                                                 remove_n_args=remove_n_args, include_full_namespace=True)
                operator_implementation_code += 2*'\n'


        # - Register code
        return_code_dict[original_class_file_name]['code_tuples'].append( (insert_pos, operator_declaration_code) )            
        return_code_dict[extras_src_file_name]['code_tuples'].append( (-1, operator_implementation_code) )            


        # - Generate a reference-returning method for each (public) member variable:
        ref_func_declaration_code    = ''
        ref_func_implementation_code = ''
        if len(member_variables) > 0:
            n_indents = len(namespaces)
            ref_func_declaration_code += '\n'
            ref_func_declaration_code += ' '*cfg.indent*(n_indents+1) + 'public:\n'
            for var_el in member_variables:

                # Put declaration in original code
                ref_func_declaration_code += classutils.constrVariableRefFunction(var_el, virtual=False, indent=cfg.indent, n_indents=n_indents+2, 
                                                                                  only_declaration=True)
                ref_func_declaration_code += '\n'

                # Put implementation in a new source file
                ref_func_implementation_code += classutils.constrVariableRefFunction(var_el, virtual=False, indent=cfg.indent, n_indents=0,
                                                                                     include_full_namespace=True) 
                ref_func_implementation_code += '\n'


        # - Register code
        if ref_func_declaration_code != '':
            return_code_dict[original_class_file_name]['code_tuples'].append( (insert_pos, ref_func_declaration_code) )            
            return_code_dict[extras_src_file_name]['code_tuples'].append( (-1, ref_func_implementation_code) )            


        # - Generate pointer-based copy and assignment functions
        n_indents = len(namespaces)
        ptr_declaration_code = '\n'
        ptr_implementation_code = '\n'

        ptr_declaration_code += ' '*cfg.indent*(n_indents+1) + 'public:\n'
        ptr_declaration_code += classutils.constrPtrCopyFunc(class_el, abstr_class_name['short'], class_name['short'], virtual=False, indent=cfg.indent, n_indents=n_indents+2, only_declaration=True)
        ptr_declaration_code += classutils.constrPtrAssignFunc(class_el, abstr_class_name['short'], class_name['short'], virtual=False, indent=cfg.indent, n_indents=n_indents+2, only_declaration=True)
        
        ptr_implementation_code += classutils.constrPtrCopyFunc(class_el, abstr_class_name['short'], class_name['short'], virtual=False, indent=cfg.indent, n_indents=0, include_full_namespace=True)
        ptr_implementation_code += classutils.constrPtrAssignFunc(class_el, abstr_class_name['short'], class_name['short'], virtual=False, indent=cfg.indent, n_indents=0, include_full_namespace=True)

        # - Register code
        return_code_dict[original_class_file_name]['code_tuples'].append( (insert_pos, ptr_declaration_code) )
        return_code_dict[extras_src_file_name]['code_tuples'].append( (-1, ptr_implementation_code) )


        # - Generate include statements for the new source file
        include_statements = utils.getIncludeStatements(class_el, convert_loaded_to='abstract', add_extra_include_path=False, input_element='class', use_full_path=True)
        include_statements.append ('#include "abstracts_typedefs.hpp"')
        include_statements.append ('#include "wrappers_typedefs.hpp"')
        if utils.isHeader(original_class_file_el):
            include_statements.append( '#include "' + os.path.join(cfg.add_path_to_includes, original_class_file_name_base) + '"')
        include_statements = list( OrderedDict.fromkeys(include_statements) )
        include_statements_code = '\n'.join(include_statements) + '\n'

        # - Register the code
        return_code_dict[extras_src_file_name]['code_tuples'].append( (0, include_statements_code) )            


        #
        # Generate factory function(s)
        #

        factory_file_content  = ''
        # if is_template and class_name['long'] in template_done:
        #     pass
        # else:
        #     original_class_file_name_base = os.path.basename(original_class_file_name)
        #     factory_file_content += '#include "' + os.path.join(cfg.add_path_to_includes, original_class_file_name_base) + '"\n'
        #     factory_file_content += '\n'
        if is_template:
            factory_file_content += classutils.constrFactoryFunction(class_el, class_name, indent=cfg.indent, template_types=spec_template_types, skip_copy_constructors=True, use_wrapper_return=False, use_wrapper_args=True)
        else:
            factory_file_content += classutils.constrFactoryFunction(class_el, class_name, indent=cfg.indent, skip_copy_constructors=True, use_wrapper_return=False, use_wrapper_args=True)
        factory_file_content += '\n'

        # - Generate factory file name
        dir_name = cfg.extra_output_dir
        factory_file_name = os.path.join(dir_name, cfg.factory_file_prefix + class_name['short'] + cfg.source_extension)

        # - Register code
        if factory_file_name not in return_code_dict.keys():
            return_code_dict[factory_file_name] = {'code_tuples':[], 'add_include_guard':False}
        return_code_dict[factory_file_name]['code_tuples'].append( (-1, factory_file_content) )


        #
        # Generate a header containing the GAMBIT wrapper class
        #

        # short_wrapper_class_fname = cfg.wrapper_header_prefix + class_name['short'] + cfg.header_extension
        # short_wrapper_class_fname = gb.new_header_files[class_name['long']]['wrapper']
        
        wrapper_decl_header_fname = gb.new_header_files[class_name['long']]['wrapper_decl']
        wrapper_def_header_fname  = gb.new_header_files[class_name['long']]['wrapper_def']
        wrapper_header_fname      = gb.new_header_files[class_name['long']]['wrapper']

        wrapper_decl_header_path = os.path.join(cfg.extra_output_dir, wrapper_decl_header_fname)
        wrapper_def_header_path  = os.path.join(cfg.extra_output_dir, wrapper_def_header_fname)
        wrapper_header_path      = os.path.join(cfg.extra_output_dir, wrapper_header_fname)
        
        # - Get code for the declaration and implementation headers
        wrapper_decl_header_content, wrapper_def_header_content = generateWrapperHeader(class_el, class_name, abstr_class_name, namespaces, 
                                                                                        short_abstr_class_fname, 
                                                                                        construct_assignment_operator, has_copy_constructor,
                                                                                        copy_constructor_id=copy_constructor_id)

        # - Code for the overall header file
        wrapper_header_content  = '\n'
        wrapper_header_content += '#include "' + wrapper_decl_header_fname + '"\n'
        wrapper_header_content += '#include "' + wrapper_def_header_fname + '"\n'
        wrapper_header_content += '\n'


        # - Register code
        if wrapper_decl_header_path not in return_code_dict.keys():
            return_code_dict[wrapper_decl_header_path] = {'code_tuples':[], 'add_include_guard':True}
        return_code_dict[wrapper_decl_header_path]['code_tuples'].append( (0, wrapper_decl_header_content) )

        if wrapper_def_header_path not in return_code_dict.keys():
            return_code_dict[wrapper_def_header_path] = {'code_tuples':[], 'add_include_guard':True}
        return_code_dict[wrapper_def_header_path]['code_tuples'].append( (0, wrapper_def_header_content) )

        if wrapper_header_path not in return_code_dict.keys():
            return_code_dict[wrapper_header_path] = {'code_tuples':[], 'add_include_guard':True}
        return_code_dict[wrapper_header_path]['code_tuples'].append( (0, wrapper_header_content) )



        #
        # Construct a function for deleting a pointer-to-wrapper and place it in the designated source & header file
        #

        wrapper_class_name = classutils.toWrapperType(class_name['long'], include_namespace=True)

        # - Include statement for the header file
        # wrapper_decl_header_path = os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, gb.new_header_files[class_name['long']]['wrapper_decl'] )
        # wrapper_include_statement_decl = '#include "' + wrapper_decl_header_path + '"\n'
        wrapper_include_statement_decl = '#include "' + gb.new_header_files[class_name['long']]['wrapper_decl_fullpath'] + '"\n'

        # - Function declaration
        w_deleter_decl  = '\n'
        w_deleter_decl += 'void wrapper_deleter(' + wrapper_class_name + '*);\n'

        # - Function implementation
        w_deleter_impl  = '\n'
        w_deleter_impl += 'void wrapper_deleter(' + wrapper_class_name + '* wptr)\n'
        w_deleter_impl += '{\n'
        w_deleter_impl += ' '*cfg.indent + 'delete wptr;\n'
        w_deleter_impl += '}\n'

        # - Register code
        w_deleter_header_path = os.path.join(cfg.extra_output_dir, gb.wrapper_deleter_fname + cfg.header_extension)
        w_deleter_source_path = os.path.join(cfg.extra_output_dir, gb.wrapper_deleter_fname + cfg.source_extension)

        if w_deleter_header_path not in return_code_dict.keys():
            return_code_dict[w_deleter_header_path] = {'code_tuples':[], 'add_include_guard':True}

            return_code_dict[w_deleter_header_path]['code_tuples'].append( (0, '#include "wrappers_typedefs.hpp"\n') )

        return_code_dict[w_deleter_header_path]['code_tuples'].append( (0, wrapper_include_statement_decl) )        
        return_code_dict[w_deleter_header_path]['code_tuples'].append( (-1, w_deleter_decl) )        

        if w_deleter_source_path not in return_code_dict.keys():
            return_code_dict[w_deleter_source_path] = {'code_tuples':[], 'add_include_guard':False}
        return_code_dict[w_deleter_source_path]['code_tuples'].append( (-1, w_deleter_impl) )        




        # #
        # # Add include statement to 'abstracts_includes.hpp'
        # #

        # # abstr_include_statement = '#include "' + os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, short_abstr_class_fname ) + '"\n'
        # abstr_include_statement = '#include "' + gb.new_header_files[class_name['long']]['abstract_fullpath'] + '"\n'

        # abstracts_includes_header_path = os.path.join(cfg.extra_output_dir, 'abstracts_includes.hpp')
        # if abstracts_includes_header_path not in return_code_dict.keys():
        #     return_code_dict[abstracts_includes_header_path] = {'code_tuples':[], 'add_include_guard':False}
        # return_code_dict[abstracts_includes_header_path]['code_tuples'].append( (0, abstr_include_statement) )


        # #
        # # Add include statements to 'wrappers_includes.hpp'
        # #

        # # wrapper_decl_header_path = os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, gb.new_header_files[class_name['long']]['wrapper_decl'] )
        # # wrapper_include_statement_decl = '#include "' + wrapper_decl_header_path + '"\n'
        # wrapper_include_statement_decl = '#include "' + gb.new_header_files[class_name['long']]['wrapper_decl_fullpath'] + '"\n'

        # wrappers_decl_includes_header_path = os.path.join(cfg.extra_output_dir, 'wrappers_decl_includes.hpp')
        # if wrappers_decl_includes_header_path not in return_code_dict.keys():
        #     return_code_dict[wrappers_decl_includes_header_path] = {'code_tuples':[], 'add_include_guard':False}
        # return_code_dict[wrappers_decl_includes_header_path]['code_tuples'].append( (0, wrapper_include_statement_decl) )        


        # # wrapper_def_header_path = os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, gb.new_header_files[class_name['long']]['wrapper_def'] )
        # # wrapper_include_statement_def = '#include "' + wrapper_def_header_path + '"\n'
        # wrapper_include_statement_def = '#include "' + gb.new_header_files[class_name['long']]['wrapper_def_fullpath'] + '"\n'

        # wrappers_def_includes_header_path = os.path.join(cfg.extra_output_dir, 'wrappers_def_includes.hpp')
        # if wrappers_def_includes_header_path not in return_code_dict.keys():
        #     return_code_dict[wrappers_def_includes_header_path] = {'code_tuples':[], 'add_include_guard':False}
        # return_code_dict[wrappers_def_includes_header_path]['code_tuples'].append( (0, wrapper_include_statement_def) )        



        #
        # Add typedef to 'abstracts_typedefs.hpp'
        #

        indent = ' '*cfg.indent*len(namespaces)
        abstr_typedef_code  = ''
        abstr_typedef_code += utils.constrNamespace(namespaces, 'open', indent=cfg.indent)

        temp_namespace_list = [gb.gambit_backend_namespace] + namespaces
        abstr_typedef_code += indent + 'typedef ' + '::'.join(temp_namespace_list) + '::' + abstr_class_name['short'] + ' ' + abstr_class_name['short'] + ';\n'

        abstr_typedef_code += utils.constrNamespace(namespaces, 'close', indent=cfg.indent)
        abstr_typedef_code += '\n'

        frw_decl_include_statement       = '#include "' + os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, gb.frwd_decls_abs_fname + cfg.header_extension) + '"\n'
        identification_include_statement = '#include "' + os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, 'identification.hpp') + '"\n\n'
        undef_include_statement          = '#include "backend_undefs.hpp"\n'

        abstracts_typedefs_header_path = os.path.join(cfg.extra_output_dir, 'abstracts_typedefs.hpp')
        if abstracts_typedefs_header_path not in return_code_dict.keys():
            return_code_dict[abstracts_typedefs_header_path] = {'code_tuples':[], 'add_include_guard':False}

            return_code_dict[abstracts_typedefs_header_path]['code_tuples'].append( (0,  frw_decl_include_statement) ) 
            return_code_dict[abstracts_typedefs_header_path]['code_tuples'].append( (len(frw_decl_include_statement), identification_include_statement) ) 
            return_code_dict[abstracts_typedefs_header_path]['code_tuples'].append( (-1, undef_include_statement) ) 

        return_code_dict[abstracts_typedefs_header_path]['code_tuples'].append( (-len(undef_include_statement), abstr_typedef_code) )


        #
        # Add typedef to 'wrappers_typdefs.hpp'
        #

        short_wrapper_class_name = classutils.toWrapperType(class_name['short'])

        wrapper_typedef_code  = ''
        wrapper_typedef_code += utils.constrNamespace(namespaces,'open')

        temp_namespace_list = [gb.gambit_backend_namespace] + namespaces
        wrapper_typedef_code += indent + 'typedef ' + '::'.join(temp_namespace_list) + '::' + class_name['short'] + ' ' + short_wrapper_class_name + ';\n'

        wrapper_typedef_code += utils.constrNamespace(namespaces,'close')
        wrapper_typedef_code += '\n'

        frw_decl_include_statement       = '#include "' + os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, gb.frwd_decls_wrp_fname + cfg.header_extension) + '"\n'
        identification_include_statement = '#include "' + os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, 'identification.hpp') + '"\n\n'
        undef_include_statement          = '#include "backend_undefs.hpp"\n'

        wrapper_typedefs_path = os.path.join(cfg.extra_output_dir, 'wrappers_typedefs.hpp')

        if wrapper_typedefs_path not in return_code_dict.keys():
            return_code_dict[wrapper_typedefs_path] = {'code_tuples':[], 'add_include_guard':False}

            return_code_dict[wrapper_typedefs_path]['code_tuples'].append( (0,  frw_decl_include_statement) ) 
            return_code_dict[wrapper_typedefs_path]['code_tuples'].append( (len(frw_decl_include_statement), identification_include_statement) ) 
            return_code_dict[wrapper_typedefs_path]['code_tuples'].append( (-1, undef_include_statement) ) 

        return_code_dict[wrapper_typedefs_path]['code_tuples'].append( (-len(undef_include_statement), wrapper_typedef_code) )



        #
        # Keep track of classes done
        #

        classes_done.append(class_name['long'])
        if is_template: 
            if class_name['long'] not in template_done:
                template_done.append(class_name['long'])
            if class_name['long'] not in templ_spec_done:
                templ_spec_done.append(class_name['long'])
        
        print
        print '~~~~ Class ' + class_name_long + ' done ~~~~'
        print
        print

    #
    # END: Loop over all classes in gb.class_dict
    #    

    # print
    # print 'CLASSES DONE                  : ', classes_done
    # print 'TEMPLATES DONE                : ', template_done
    # print 'TEMPLATE SPECIALIZATIONS DONE : ', templ_spec_done
    # print


    # Increase class counter
    gb.n_classes_done += 1

    #
    # Return result 
    #

    return return_code_dict

#
# END: run()
#

#
# -------------------------------------------------------------------
#

#
# Function for generating a header file with a GAMBIT wrapper class
#
def generateWrapperHeader(class_el, class_name, abstr_class_name, namespaces, 
                          short_abstr_class_fname,
                          construct_assignment_operator, has_copy_constructor,
                          copy_constructor_id=''):

    # Useful variables
    indent = ' '*cfg.indent
    # wrapper_base_class_name = 'WrapperBase<' + abstr_class_name['long'] + '>'

    # Useful lists
    class_variables    = []
    class_functions    = utils.getMemberFunctions(class_el, include_artificial=False, include_inherited=cfg.wrap_inherited_members, 
                                                            only_accepted=True, limit_pointerness=True, include_operators=True)
    class_constructors = []
    class_members      = utils.getMemberElements(class_el, include_artificial=False)
    class_members_full = utils.getMemberElements(class_el, include_artificial=True)
    for mem_el in class_members:
        # if (mem_el.tag == 'Method') and (mem_el.get('access') == 'public'):
        #     if funcutils.ignoreFunction(mem_el, limit_pointerness=True):
        #         warnings.warn('The member "%s" in class "%s" makes use of a non-accepted type and will be ignored.' % (mem_el.get('name'), class_name['short']))
        #     else:
        #         class_functions.append(mem_el)

        if (mem_el.tag in ('Field', 'Variable')) and (mem_el.get('access') == 'public'):
            if utils.isAcceptedType(mem_el):
                class_variables.append(mem_el)
            else:
                # warnings.warn('The member "%s" in class "%s" makes use of a non-accepted type and will be ignored.' % (mem_el.get('name'), class_name['short']))
                print 'INFO: ' + 'The member "%s" in class "%s" makes use of a non-accepted type and will be ignored.' % (mem_el.get('name'), class_name['short'])
        else:
            pass
    for mem_el in class_members_full:

        # Skip the copy constructor
        if has_copy_constructor and (mem_el.get('id') == copy_constructor_id):
            continue

        # Store constructor if acceptable
        if (mem_el.tag == 'Constructor') and (mem_el.get('access') == 'public'):
            if funcutils.ignoreFunction(mem_el, limit_pointerness=True):
                # warnings.warn('The constructor "%s" in class "%s" makes use of a non-accepted type and will be ignored.' % (mem_el.get('name'), class_name['short']))
                print 'INFO: ' + 'The constructor "%s" in class "%s" makes use of a non-accepted type and will be ignored.' % (mem_el.get('name'), class_name['short'])
            else:
                class_constructors.append(mem_el)


    # Create a list of dicts with info on the (loaded) parent classes
    loaded_parent_classes = utils.getParentClasses(class_el, only_loaded_classes=True)


    #
    # Start code generation
    #

    decl_code = classutils.constrWrapperDecl(class_name, abstr_class_name, loaded_parent_classes, class_variables, class_functions, class_constructors, has_copy_constructor, indent=indent)

    def_code  = classutils.constrWrapperDef(class_name, abstr_class_name, loaded_parent_classes, class_variables, class_functions, class_constructors, has_copy_constructor, indent=indent, do_inline=True)


    # Insert tags for the GAMBIT namespace
    decl_code = '\n__START_GAMBIT_NAMESPACE__\n' + decl_code + '\n__END_GAMBIT_NAMESPACE__\n'
    def_code  = '\n__START_GAMBIT_NAMESPACE__\n' + def_code  + '\n__END_GAMBIT_NAMESPACE__\n'

    # Insert include statements needed by GAMBIT 
    decl_code = '#include "identification.hpp"\n' + decl_code + '\n#include "backend_undefs.hpp"\n'
    def_code  = '#include "identification.hpp"\n' + def_code + '\n#include "backend_undefs.hpp"\n'


    
    #
    # Add #include statements for the declaration code
    #

    decl_code_include_statements = []

    # - Header where NULL is defined
    decl_code_include_statements.append( '#include <cstddef>' )

    # - Header where the function nullptr_check is defined
    # decl_code_include_statements.append( '#include "nullptr_check' + cfg.header_extension + '"')

    # - Header with forward declarations to all wrapper classes
    decl_code_include_statements.append( '#include "' + gb.frwd_decls_wrp_fname + cfg.header_extension + '"')

    # - Base class for all wrapper classes
    decl_code_include_statements.append( '#include "' + 'wrapperbase.hpp"')
    # decl_code_include_statements.append( '#include "' + gb.new_header_files['WrapperBase']['wrapper'] + '"')

    # - Abstract class for the original class
    decl_code_include_statements.append( '#include "' + gb.new_header_files[class_name['long']]['abstract'] + '"' )
    # decl_code_include_statements.append( '#include "' + os.path.join( cfg.add_path_to_includes, gb.new_header_files[class_name['long']]['abstract']) + '"' )

    # - Wrapper parent classes
    for parent_dict in loaded_parent_classes:
        decl_code_include_statements.append('#include "' + gb.new_header_files[ parent_dict['class_name']['long'] ]['wrapper_decl'] + '"')

    # - Any other types (excluding the current wrapper class)
    # decl_code_include_statements += classutils.getIncludeStatements(all_types_in_class, 'wrapper', add_extra_include_path=False)
    # decl_code_include_statements += utils.getIncludeStatements(class_el, convert_loaded_to='wrapper', add_extra_include_path=False, exclude_types=[class_name])
    decl_code_include_statements += utils.getIncludeStatements(class_el, convert_loaded_to='wrapper_decl', add_extra_include_path=False, exclude_types=[class_name], skip_forward_declared=True)

    # Remove duplicates and construct code
    decl_code_include_statements = list( OrderedDict.fromkeys(decl_code_include_statements) )
    decl_include_statements_code = '\n'.join(decl_code_include_statements) + 2*'\n'
    decl_code = decl_include_statements_code + decl_code



    #
    # Add #include statements for the declaration code
    #

    def_code_include_statements = []

    # - Any other types (excluding the current wrapper class)
    # def_code_include_statements += classutils.getIncludeStatements(all_types_in_class, 'wrapper', add_extra_include_path=False)
    # def_code_include_statements += utils.getIncludeStatements(class_el, convert_loaded_to='wrapper', add_extra_include_path=False, exclude_types=[class_name])
    def_code_include_statements += utils.getIncludeStatements(class_el, convert_loaded_to='wrapper_decl', add_extra_include_path=False, exclude_types=[class_name], skip_forward_declared=False)

    # Remove duplicates and construct code
    def_code_include_statements = list( OrderedDict.fromkeys(def_code_include_statements) )
    def_include_statements_code = '\n'.join(def_code_include_statements) + 2*'\n'
    def_code = def_include_statements_code + def_code


    # Return code
    return decl_code, def_code

#
# END: generateWrapperHeader
#   
