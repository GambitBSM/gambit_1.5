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
    
    for class_name_long, class_el in cfg.class_dict.items():

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
            cfg.class_dict.pop(class_name_long)
            continue

        # Make list of all types in use in this class
        all_types_in_class = utils.getAllTypesInClass(class_el, include_parents=True)

        # Set a bunch of generally useful variables 
        original_class_file_el        = cfg.id_dict[class_el.get('file')]
        original_class_file_name      = original_class_file_el.get('name')
        original_class_file_name_base = os.path.basename(original_class_file_name)
        original_class_file_dir       = os.path.split(original_class_file_name)[0]
        extras_src_file_name          = os.path.join(cfg.extra_output_dir, class_name['short'] + '_extras' + cfg.code_suffix + cfg.source_extension)

        short_abstr_class_fname = cfg.new_header_files[class_name['long']]['abstract']
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
                if (template_type not in cfg.accepted_types):
                    raise Exception("The template specialization type '" + template_type + "' for class " + class_name['long'] + " is not among accepted types.")


        #
        # Construct code for the abstract class header file and register it
        #
        
        class_decl = ''

        # - Add include statements
        include_statements  = []
        include_statements  = ['#include "AbstractBase' + cfg.header_extension + '"']
        include_statements += ['#include "' + cfg.frwd_decls_abs_fname + cfg.header_extension + '"']
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
        include_line = '#include "' + short_abstr_class_fname + '"'
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
                el = cfg.id_dict[mem_id]
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
                    # return_is_native = utils.isNative( cfg.id_dict[return_id] )
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
        include_statements = utils.getIncludeStatements(class_el, convert_loaded_to='abstract', add_extra_include_path=True, input_element='class')
        include_statements = utils.getIncludeStatements(class_el, convert_loaded_to='abstract', add_extra_include_path=True, input_element='class')
        if utils.isHeader(original_class_file_el):
            include_statements.append( '#include "' + os.path.join(cfg.add_path_to_includes, original_class_file_name_base) + '"')
        include_statements = list(set(include_statements))
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
            factory_file_content += classutils.constrFactoryFunction(class_el, class_name, indent=cfg.indent, template_types=spec_template_types, skip_copy_constructors=True, use_wrapper_class=True)
        else:
            factory_file_content += classutils.constrFactoryFunction(class_el, class_name, indent=cfg.indent, skip_copy_constructors=True, use_wrapper_class=True)
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
        short_wrapper_class_fname = cfg.new_header_files[class_name['long']]['wrapper']
        wrapper_class_fname = os.path.join(cfg.extra_output_dir, short_wrapper_class_fname)
        
        wrapper_class_file_content = generateWrapperHeader(class_el, class_name, abstr_class_name, namespaces, 
                                                           short_abstr_class_fname, short_wrapper_class_fname,
                                                           construct_assignment_operator, has_copy_constructor,
                                                           copy_constructor_id=copy_constructor_id)

        # - Update dict of wrapper class header files
        #wrapper_class_headers[class_name['short']] = wrapper_class_fname

        # - Register code
        if wrapper_class_fname not in return_code_dict.keys():
            return_code_dict[wrapper_class_fname] = {'code_tuples':[], 'add_include_guard':True}
        return_code_dict[wrapper_class_fname]['code_tuples'].append( (0, wrapper_class_file_content) )


        #
        # Add wrapper header to 'master header' (header that includes all wrapper classes)
        #

        master_header_include = '#include "' + short_wrapper_class_fname + '"\n'

        all_wrapper_header_path = os.path.join(cfg.extra_output_dir, cfg.all_wrapper_fname + cfg.header_extension)          
        if all_wrapper_header_path not in return_code_dict.keys():
            return_code_dict[all_wrapper_header_path] = {'code_tuples':[], 'add_include_guard':True}
        return_code_dict[all_wrapper_header_path]['code_tuples'].append( (0, master_header_include) )        


        #
        # Construct a function for deleting a pointer-to-wrapper and place it in the designated source & header file
        #

        short_wrapper_class_name = classutils.toWrapperType(class_name['short'])

        # - Function declaration
        w_deleter_decl  = '\n'
        w_deleter_decl += 'void wrapper_deleter(' + short_wrapper_class_name + '*);\n'

        # - Function implementation
        w_deleter_impl  = '\n'
        w_deleter_impl += 'void wrapper_deleter(' + short_wrapper_class_name + '* wptr)\n'
        w_deleter_impl += '{\n'
        w_deleter_impl += ' '*cfg.indent + 'delete wptr;\n'
        w_deleter_impl += '}\n'

        # - Register code
        w_deleter_header_path = os.path.join(cfg.extra_output_dir, cfg.wrapper_deleter_fname + cfg.header_extension)
        w_deleter_source_path = os.path.join(cfg.extra_output_dir, cfg.wrapper_deleter_fname + cfg.source_extension)

        if w_deleter_header_path not in return_code_dict.keys():
            return_code_dict[w_deleter_header_path] = {'code_tuples':[], 'add_include_guard':True}
        return_code_dict[w_deleter_header_path]['code_tuples'].append( (-1, w_deleter_decl) )        

        if w_deleter_source_path not in return_code_dict.keys():
            return_code_dict[w_deleter_source_path] = {'code_tuples':[], 'add_include_guard':False}
        return_code_dict[w_deleter_source_path]['code_tuples'].append( (-1, w_deleter_impl) )        



        #
        # Add typedef to the wrapper_typedefs file (example: typedef SomeClass_gambit SomeClass; )
        #

        wrapper_class_name = classutils.toWrapperType(class_name['short'])
        typedef_code = 'typedef ' + wrapper_class_name + ' ' + class_name['short'] + ';\n'

        wrapper_typedefs_path = cfg.wrapper_typedefs_fname + cfg.header_extension

        if wrapper_typedefs_path not in return_code_dict.keys():
            return_code_dict[wrapper_typedefs_path] = {'code_tuples':[], 'add_include_guard':True}
        return_code_dict[wrapper_typedefs_path]['code_tuples'].append( (-1, typedef_code) )


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
    # END: Loop over all classes in cfg.class_dict
    #    

    # print
    # print 'CLASSES DONE                  : ', classes_done
    # print 'TEMPLATES DONE                : ', template_done
    # print 'TEMPLATE SPECIALIZATIONS DONE : ', templ_spec_done
    # print


    #
    # Return result
    #

    return return_code_dict

#
# END: run()
#



#
# Function for generating a header file with a GAMBIT wrapper class
#
def generateWrapperHeader(class_el, class_name, abstr_class_name, namespaces, 
                          short_abstr_class_fname, short_wrapper_class_fname,
                          construct_assignment_operator, has_copy_constructor,
                          copy_constructor_id=''):

    # Useful variables
    indent = ' '*cfg.indent
    wrapper_base_class_name = 'WrapperBase<' + abstr_class_name['long'] + '>'

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

    code = ''


    #
    # --- UPDATE: The code that declares factory function pointers is currently commented out
    #
    # # Create one factory function pointer for each constructor.
    # # Initialize pointers to NULL - they will later be filled by dynamic loading
    # # We append a number (e.g. '_1') to the function pointer name to make sure each overload of a constructor
    # # gets a unique.
    # temp_code = ''
    # for i, constr_el in enumerate(class_constructors):

    #     # We need pointers for all the overloaded factory functions (generated due to default value arguments)
    #     n_overloads = funcutils.numberOfDefaultArgs(constr_el)

    #     # Identify arguments, translate argument type of loaded classes
    #     # and construct the argument bracket
    #     args = funcutils.getArgs(constr_el)
    #     w_args = funcutils.constrWrapperArgs(args, add_ref=True)

    #     # One factory function pointer for each set of default arguments
    #     for remove_n_args in range(n_overloads+1):

    #         if remove_n_args == 0:
    #             use_w_args = w_args
    #         else:
    #             use_w_args = w_args[:-remove_n_args]

    #         args_bracket = funcutils.constrArgsBracket(use_w_args, include_arg_name=False, include_arg_type=True, include_namespace=True)

    #         # Factory pointer name
    #         factory_ptr_name = 'Factory_' + class_name['short'] + '_' + str(i)
    #         if remove_n_args > 0:
    #             factory_ptr_name += '_overload_' + str(remove_n_args)

    #         # Construct factory pointer code
    #         temp_code += abstr_class_name['long'] + '* (*' + factory_ptr_name + ')' + args_bracket + ' = NULL;\n'

    # if temp_code != '':
    #     code += '\n'
    #     code += "// Factory function pointers to be filled by dynamic loading\n"
    #     code += temp_code


    #
    # Generate code for wrapper class
    #

    short_wrapper_class_name = class_name['short'] + cfg.code_suffix

    # Construct line declaring inheritance from the wrapper base class (WrapperBase)
    # and from other wrapper classes
    inheritance_line = 'public ' + wrapper_base_class_name
    for parent_dict in loaded_parent_classes:

        inheritance_line += 'virtual '*parent_dict['virtual'] + parent_dict['access'] + ' ' + parent_dict['wrapper_name'] + ', '

    if inheritance_line != '':
        inheritance_line = ' : ' + inheritance_line.rstrip(', ')

    # Class declaration line
    code += '\n'
    code += 'class ' + short_wrapper_class_name + inheritance_line + '\n'

    # Class body
    code += '{\n'

    # # Start private block
    # code += indent + 'private:\n'
    
    # # Private variable: bool member_element
    # code += 2*indent + 'bool member_variable;\n'

    # Start public block
    code += indent + 'public:\n'

    # Variables:
    code += 2*indent + '// Member variables: \n'

    # # Public variable: pointer to abstract class: Abstract_ClassX* BEptr
    # code += 2*indent + abstr_class_name['long'] + '* BEptr;\n'

    # Add references to all public variables
    for var_el in class_variables:

        # Variable name
        var_name = var_el.get('name')

        # Determine variable type
        var_type, var_type_kw, var_type_id = utils.findType(var_el)
        var_is_loaded_class = utils.isLoadedClass(var_type, byname=True)

        if var_is_loaded_class:
            use_var_type = classutils.toWrapperType(var_type)
        else:
            use_var_type = var_type

        # Remove '&' from use_var_type if it exists
        use_var_type = use_var_type.replace('&','')

        # Write line
        # FIXME: Should we include keywords here (const, static, ...)?
        if var_is_loaded_class:
            code += 2*indent + use_var_type + ' ' + var_name + ';\n'
        else:
            code += 2*indent + use_var_type + '& ' + var_name + ';\n'


    # Functions:
    code += '\n'
    code += 2*indent + '// Member functions: \n'

    # Add wrappers for all member functions, including operator functions
    # and overloaded versions of functions with default value arguments

    # class_functions = funcutils.multiplyDefaultArgFunctions(class_functions)
    # done_members = []

    for func_el in class_functions:

        # Check if this is an operator function
        is_operator = False
        if func_el.tag == 'OperatorMethod':
            is_operator = True

        # Check if this function makes use of any loaded types
        uses_loaded_type = funcutils.usesLoadedType(func_el)

        # # Check if we're generating an overload of a previous member function
        # overload_n = done_members.count(func_el)
        # is_overload = bool(overload_n)

        # Function name
        if is_operator:
            func_name = 'operator' + func_el.get('name')
        else:
            func_name = func_el.get('name')

        # Check constness
        if ('const' in func_el.keys()) and (func_el.get('const')=='1'):
            is_const = True
        else:
            is_const = False

        # Determine return type
        return_type, return_type_kw, return_type_id = utils.findType(func_el)
        return_type_el = cfg.id_dict[return_type_id]

        pointerness, is_ref = utils.pointerAndRefCheck(return_type, byname=True)
        return_is_loaded    = utils.isLoadedClass(return_type, byname=True)

        # Convert return type if loaded class
        if return_is_loaded:
            # use_return_type = classutils.toWrapperType(return_type, remove_reference=True)
            use_return_type = classutils.toWrapperType(return_type)
        else:
            use_return_type = return_type

        # If return-by-value, then a const qualifier on the return value is meaningless
        # (will result in a compiler warning)
        if utils.pointerAndRefCheck(use_return_type, byname=True) == (0, False):
            if 'const' in return_type_kw:
                return_type_kw.remove('const')

        # Return keywords
        return_kw_str = ' '.join(return_type_kw)        
        return_kw_str += ' '*bool(len(return_type_kw))

        # Arguments
        args = funcutils.getArgs(func_el)


        # One function for each set of default arguments
        n_overloads = funcutils.numberOfDefaultArgs(func_el)
        for remove_n_args in range(n_overloads+1):

            if remove_n_args == 0:
                use_args = args
            else:
                use_args = args[:-remove_n_args]

            # For operators, take wrapper class types as input. For functions, take WrapperBase as input.
            if is_operator:
                args_bracket = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=True, include_namespace=True, use_wrapper_class=True)
            else:
                args_bracket = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=True, include_namespace=True, use_wrapper_class=True, use_wrapper_base_class=True)                

            # Name of function to call (in abstract class)
            # if (remove_n_args > 0) and (not is_operator):
            #     call_func_name = func_name + cfg.code_suffix
            # else:
            #     call_func_name = func_name
            if is_operator:
                call_func_name = 'operator_' + cfg.operator_names[func_el.get('name')] + cfg.code_suffix
            else:
                # call_func_name = func_name + cfg.code_suffix
                if uses_loaded_type or (remove_n_args>0):
                    call_func_name = func_name + cfg.code_suffix 
                else:
                    call_func_name = func_name


            # Write declaration line
            code += 2*indent + return_kw_str + use_return_type + ' ' + func_name + args_bracket + is_const*' const' + '\n'

            # Write function body
            code += 2*indent + '{\n'

            if return_type == 'void':
                code += 3*indent
            else:
                code += 3*indent + 'return '

            args_bracket_notypes = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=False, wrapper_to_pointer=True)

            if return_is_loaded: 
                if is_ref:
                    # return reference_returner<X_GAMBIT, Abstract_X>( BEptr->return_ref_this_GAMBIT() );
                    # code += use_return_type + '( BEptr->' + call_func_name + args_bracket_notypes + ' , true);\n'
                    abs_return_type_for_templ = classutils.toAbstractType(return_type, include_namespace=True, remove_reference=True, remove_pointers=True)
                    w_return_type_for_templ = classutils.toWrapperType(return_type, remove_reference=True, remove_pointers=True)
                    code += 'reference_returner< ' + w_return_type_for_templ + ', ' + abs_return_type_for_templ +  ' >( BEptr->' + call_func_name + args_bracket_notypes + ' );\n'
                elif (not is_ref) and (pointerness > 0):
                    abs_return_type_for_templ = classutils.toAbstractType(return_type, include_namespace=True, remove_reference=True, remove_pointers=True)
                    w_return_type_for_templ = classutils.toWrapperType(return_type, remove_reference=True, remove_pointers=True)
                    code += 'pointer_returner< ' + w_return_type_for_templ + ', ' + abs_return_type_for_templ +  ' >( BEptr->' + call_func_name + args_bracket_notypes + ' );\n'
                else:
                    code += use_return_type + '( BEptr->' + call_func_name + args_bracket_notypes + ' );\n'
            else:                
                code += 'BEptr->' + call_func_name + args_bracket_notypes + ';\n'

            code += 2*indent + '}\n'
            code += '\n'

        # # Keep track of functions done
        # done_members.append(func_el)

    # # Add special member function: _set_member(bool) - set the variable 'member_variable' to true/false
    # code += '\n'
    # code += 2*indent + '// Special member function to set member_variable: \n'
    # code += 2*indent + 'void _set_member_variable(bool in) { member_variable = in; }\n'    


    #
    # Add all constructors here...
    #

    # First generate some code common to all constructors
    common_init_list_code = ''
    for var_el in class_variables:
        var_name = var_el.get('name')
        var_type, var_type_kw, var_type_id = utils.findType(var_el)
        if utils.isLoadedClass(var_type, byname=True):
            common_init_list_code += 3*indent + var_name + '(&(BEptr->' + var_name + '_ref' + cfg.code_suffix + '())),\n'
        else:
            common_init_list_code += 3*indent + var_name + '(BEptr->' + var_name + '_ref' + cfg.code_suffix + '()),\n'
    if common_init_list_code != '':
        common_init_list_code = common_init_list_code.rstrip(',\n') + '\n'

    common_constructor_body = ''
    common_constructor_body += 2*indent + '{\n'
    common_constructor_body += 3*indent + 'BEptr->wrapper' + cfg.code_suffix + '(this);\n'

    for var_el in class_variables:
        var_name = var_el.get('name')
        var_type, var_type_kw, var_type_id = utils.findType(var_el)
        if utils.isLoadedClass(var_type, byname=True):
            common_constructor_body += 3*indent + var_name + '._setMemberVariable(true);\n'

    common_constructor_body += 2*indent + '}\n'



    #
    # --- UPDATE: The code for generating wrapper class constructors is currently commented out
    #

    # # Add wrappers for all original constructors except the copy constructor
    # temp_code = ''
    # for i, constr_el in enumerate(class_constructors):

    #     # Identify arguments
    #     args = funcutils.getArgs(constr_el)
    #     factory_args = funcutils.constrWrapperArgs(args, add_ref=True)

    #     # Skip if this is a copy constructor. (Another copy constructor is added below.)
    #     if (len(args) == 1) and (args[0]['id'] == class_el.get('id')):
    #         continue

    #     # If default arguments are use, we need overloaded constructors to connect to the overloaded
    #     # factory function pointers
    #     n_overloads = funcutils.numberOfDefaultArgs(constr_el)

    #     # One constructor for each set of default arguments
    #     for remove_n_args in range(n_overloads+1):

    #         if remove_n_args == 0:
    #             use_args         = args
    #             factory_use_args = factory_args
    #         else:
    #             use_args         = args[:-remove_n_args]
    #             factory_use_args = factory_args[:-remove_n_args]

    #         args_bracket = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=True, include_namespace=True, use_wrapper_class=True)
    #         args_bracket_notypes = funcutils.constrArgsBracket(use_args, include_arg_name=True, include_arg_type=False, wrapper_to_pointer=True)
    #         factory_args_bracket = funcutils.constrArgsBracket(factory_use_args, include_arg_name=False, include_arg_type=True, include_namespace=True)

    #         # Factory pointer name
    #         factory_ptr_name = 'Factory_' + class_name['short'] + '_' + str(i)
    #         if remove_n_args > 0:
    #             factory_ptr_name += '_overload_' + str(remove_n_args)

    #         temp_code += 2*indent + short_wrapper_class_name + args_bracket + ' :\n'
    #         # temp_code += 3*indent + 'BEptr( Factory_' + class_name['short'] + '_' + str(i) + args_bracket_notypes + ' ),\n'  # FIXME: This is not general. Fix argument list.
    #         temp_code += 3*indent + wrapper_base_class_name + '( ' + 'nullptr_check< '+ abstr_class_name['long'] + '* (*)' + factory_args_bracket + ' >(' + factory_ptr_name +')'+ args_bracket_notypes + ', ' + 'false' + ' )'  # FIXME: This is not general. Fix argument list.
    #         if common_init_list_code != '':
    #             temp_code += ',\n' + common_init_list_code
    #         else:
    #             temp_code += '\n'
    #         temp_code += common_constructor_body

    #         temp_code += '\n'

    # if temp_code != '':
    #     code += '\n'
    #     code += 2*indent + '// Wrappers for original constructors: \n'    
    #     code += temp_code


    # Add special constructor based on abstract pointer (This is needed to allow return-by-value with the wrapper classes.)
    code += 2*indent + '// Special pointer-based constructor: \n'
    code += 2*indent + short_wrapper_class_name + '(' + abstr_class_name['long'] +'* in, bool memvar_in=false) :\n'
    # code += 3*indent + 'BEptr(in),\n'
    code += 3*indent + wrapper_base_class_name + '( in, memvar_in )'  # FIXME: This is not general. Fix argument list.
    if common_init_list_code != '':
        code += ',\n' + common_init_list_code
    else:
        code += '\n'
    code += common_constructor_body

    # Add copy constructor
    if has_copy_constructor:
        code += '\n'
        code += 2*indent + '// Copy constructor: \n'
        code += 2*indent + short_wrapper_class_name + '(const ' + short_wrapper_class_name +'& in) :\n'
        # code += 3*indent + 'BEptr(in),\n'
        code += 3*indent + wrapper_base_class_name + '(in)'  # FIXME: This is not general. Fix argument list.
        if common_init_list_code != '':
            code += ',\n' + common_init_list_code
        else:
            code += '\n'
        code += common_constructor_body


    # # Add copy constructor -- UPDATE: now placed in WrapperBase
    # if has_copy_constructor:
    #     code += '\n'
    #     code += 2*indent + '// Copy constructor: \n'
    #     code += 2*indent + short_wrapper_class_name + '(const ' + short_wrapper_class_name +'& in) :\n'
    #     # code += 3*indent + 'BEptr(in.BEptr->pointerCopy' + cfg.code_suffix + '()),\n'
    #     code += 3*indent + wrapper_base_class_name + '( in.BEptr->pointerCopy' + cfg.code_suffix + '(), false )'  # FIXME: This is not general. Fix argument list.
    #     if common_init_list_code != '':
    #         code += ',\n' + common_init_list_code
    #     else:
    #         code += '\n'
    #     code += common_constructor_body


    #
    # Add assignment operator -- UPDATE: now placed in WrapperBase
    #
    # if construct_assignment_operator:
    #     code += '\n'
    #     code += 2*indent + '// Assignment operator: \n'
    #     code += 2*indent + short_wrapper_class_name + '& ' + 'operator=(const ' + short_wrapper_class_name +'& in)\n'
    #     code += 2*indent + '{\n'
    #     # code += 3*indent + 'if (this != &in) { BEptr->pointerAssign' + cfg.code_suffix + '(in.BEptr.get()); }\n'
    #     code += 3*indent + 'if (this != &in) { BEptr->pointerAssign' + cfg.code_suffix + '(in.BEptr); }\n'
    #     code += 3*indent + 'return *this;'
    #     code += 2*indent + '}\n'

    # #
    # # Add destructor
    # #
    # code += '\n'
    # code += 2*indent + '// Destructor: \n'
    # code += 2*indent + '~' + short_wrapper_class_name + '()\n'
    # code += 2*indent + '{\n'
    # code += 3*indent + 'if(member_variable==false) { delete BEptr; }\n'
    # code += 2*indent + '}\n'

    # Close class body
    code += '};\n'

    
    # Add #include statements at the beginning
    include_statements = []

    # - Header where NULL is defined
    include_statements.append( '#include <cstddef>' )

    # # - Header where the function nullptr_check is defined
    # include_statements.append( '#include "nullptr_check' + cfg.header_extension + '"')

    # - Base class for all wrapper classes
    include_statements.append( '#include "' + cfg.wrapper_header_prefix + 'WrapperBase' + cfg.header_extension + '"')
    # include_statements.append( '#include "' + cfg.new_header_files['WrapperBase']['wrapper'] + '"')

    # - Abstract class for the original class
    include_statements.append( '#include "' + cfg.new_header_files[class_name['long']]['abstract'] + '"' )
    # include_statements.append( '#include "' + os.path.join( cfg.add_path_to_includes, cfg.new_header_files[class_name['long']]['abstract']) + '"' )

    # - Wrapper parent classes
    for parent_dict in loaded_parent_classes:
        include_statements.append('#include "' + cfg.new_header_files[ parent_dict['class_name']['long'] ]['wrapper'] + '"')

    # - Any other types (excluding the current wrapper class)
    # include_statements += classutils.getIncludeStatements(all_types_in_class, 'wrapper', add_extra_include_path=False)
    include_statements += utils.getIncludeStatements(class_el, convert_loaded_to='wrapper', add_extra_include_path=False, exclude_types=[class_name])


    # Remove duplicates and construct code
    include_statements = list( set(include_statements) )
    include_statements_code = '\n'.join(include_statements) + 2*'\n'
    code = include_statements_code + code


    # Add include guards
    # code = utils.addIncludeGuard(code, short_wrapper_class_fname)

    # Return code
    return code

#
# END: generateWrapperHeader
#   
