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
#wrapper_class_headers = {}

#
# Main function for parsing classes
#

def run():

    # Prepare returned dict
    return_code_dict = OrderedDict()

    # Prepare dict to keep track of include statements in source files
    src_includes = OrderedDict()

    #
    # Loop over all classes 
    #
    
    for class_name_full, class_el in cfg.class_dict.items():

        # Print current class
        print
        print '~~~~ Current class: ' + class_name_full + ' ~~~~'
        print

        # Check if we've done this class already
        if class_name_full in classes_done:
            print ' '*15 + '--> Class already done'
            continue


        # NOT NEEDED WHEN PARSING ONLY LISTED CLASSES:
        #
        # # Check if this class is native to the source code
        # if utils.isNative(class_el):
        #     if 'extern' in class_el.keys() and class_el.get('extern') == "1":
        #         print ' '*15 + '--> Class %s RECOGNIZED but IGNORED' % class_name_full
        #         continue
        #     else:
        #         print ' '*15 + '--> Class %s ACCEPTED' % class_name_full
        #         pass
        # else:
        #     print ' '*15 + '--> Class IGNORED'
        #     continue


        # Check if this is a template class
        if '<' in class_name_full:
            is_template = True
        else:
            is_template = False


        # Make list of all types in use in this class
        # Also, make a separate list of all std types
        check_member_elements = utils.getMemberElements(class_el)
        cfg.all_types_in_class = {};
        cfg.std_types_in_class = {};

        class_id = class_el.get('id')
        for mem_el in check_member_elements:

            if mem_el.tag in ['Constructor', 'Destructor', 'Method', 'OperatorMethod']:
                args_list = funcutils.getArgs(mem_el)
                for arg_dict in args_list:
                    arg_type_el = cfg.id_dict[arg_dict['id']]

                    arg_type = arg_dict['type']

                    if arg_type not in cfg.all_types_in_class.keys():
                        cfg.all_types_in_class[arg_type] = arg_type_el

                        if utils.isStdType(arg_type_el):
                            cfg.std_types_in_class[arg_type] = arg_type_el


            if ('type' in mem_el.keys()) or ('returns' in mem_el.keys()):
                mem_type, mem_type_kw, mem_type_id = utils.findType(mem_el)
                type_el = cfg.id_dict[mem_type_id]

                if mem_type not in cfg.all_types_in_class.keys():
                    cfg.all_types_in_class[mem_type] = type_el

                    if utils.isStdType(type_el):
                        cfg.std_types_in_class[mem_type] = type_el


                # if utils.isStdType(type_el):
                #     if mem_type not in cfg.std_types_in_class:
                #         cfg.std_types_in_class.append(mem_type)

        # while len(check_member_elements)>0:
        #     print 'lenght: ', len(check_member_elements)
        #     el = check_member_elements[-1]

        #     if 'type' in el.keys():
        #         mem_type, mem_type_kw, mem_type_id = utils.findType(el)
        #         if mem_type not in cfg.types_used:
        #             cfg.types_used.append(mem_type) 

        #     append_elements = utils.getMemberElements(el)
        #     check_member_elements = append_elements + check_member_elements
        #     check_member_elements.pop()


        
        # Set a bunch of generally useful variables 
        class_name             = class_name_full.split('<',1)[0]
        short_class_name       = class_el.get('name').split('<',1)[0]
        short_abstr_class_name = classutils.getAbstractClassName(class_name, prefix=cfg.abstr_class_prefix, short=True)

        src_file_el        = cfg.id_dict[class_el.get('file')]
        src_file_name      = src_file_el.get('name')
        src_dir            = os.path.split(src_file_name)[0]

        short_abstr_class_fname = cfg.abstr_header_prefix + short_class_name + cfg.header_extension
        abstr_class_fname       = os.path.join(cfg.extra_output_dir, short_abstr_class_fname)

        namespaces    = class_name.split('::')[:-1]
        has_namespace = bool(len(namespaces))

        has_copy_constructor    = classutils.checkCopyConstructor(class_el)
        has_assignment_operator = classutils.checkAssignmentOperator(class_el)


        # Prepare entries in return_code_dict and src_includes
        if abstr_class_fname not in return_code_dict.keys():
            return_code_dict[abstr_class_fname] = []
        if src_file_name not in return_code_dict.keys():
            return_code_dict[src_file_name] = []
        if src_file_name not in src_includes.keys():
            src_includes[src_file_name] = []


        # Register an include statement for this header file, to go in the file cfg.all_headers_fname
        include_line = '#include "' + short_abstr_class_fname + '"\n' 
        all_headers_file_path = os.path.join(cfg.extra_output_dir, cfg.all_headers_fname)

        if all_headers_file_path not in return_code_dict.keys():
            return_code_dict[all_headers_file_path] = []
        if include_line not in [ inc_tuple[1] for inc_tuple in return_code_dict[all_headers_file_path] ]:
            return_code_dict[all_headers_file_path].append( (0, include_line) )


        # Treat the first specialization of a template class differently
        if is_template and class_name not in template_done:
            template_bracket, template_types = utils.getTemplateBracket(class_el)
            
            empty_templ_class_decl = ''
            empty_templ_class_decl += classutils.constrEmptyTemplClassDecl(short_abstr_class_name, namespaces, template_bracket, indent=cfg.indent)
            empty_templ_class_decl += classutils.constrTemplForwDecl(short_class_name, namespaces, template_bracket, indent=cfg.indent)

            return_code_dict[abstr_class_fname].append( (0, empty_templ_class_decl) )


        # Get template arguments for specialization, 
        # and check that they are acceptable
        if is_template and class_name_full not in templ_spec_done:
            spec_template_types = utils.getSpecTemplateTypes(class_el)
            for template_type in spec_template_types:
                if (template_type not in cfg.accepted_types):
                    raise Exception("The template specialization type '" + template_type + "' for class " + class_name_full + " is not among accepted types.")


        # Construct code for the abstract class header file and register it

        class_decl = ''

        # - DEBUG: For debug purposes, include <iostream> in all abstract header files
        #          (To allow virtual member functions to print a warning if they are executed.)
        class_decl += '#include <iostream>  // FOR DEBUG: Allow virtual member functions to print a warning if executed.\n'
        class_decl += '\n'

        # - Add forward declarations at the top of the header file  
        #
        #   FIXME: For now, this includes forward declarations of *all* loaded native classes.
        #          Should limit this to only the forward declarations needed for the given abstract class.
        class_decl += '// Forward declarations:\n'
        class_decl += utils.constrForwardDecls()
        class_decl += '\n'


        # - Add the the code for the abstract class
        if (is_template == True) and (class_name_full in templ_spec_done):
            pass
        elif (is_template == True) and (class_name_full not in templ_spec_done):
            class_decl += classutils.constrAbstractClassDecl(class_el, short_class_name, short_abstr_class_name, namespaces, 
                                                                indent=cfg.indent, template_types=spec_template_types, 
                                                                has_copy_constructor=has_copy_constructor, has_assignment_operator=has_assignment_operator)
            class_decl += '\n'
        else:
            class_decl += classutils.constrAbstractClassDecl(class_el, short_class_name, short_abstr_class_name, namespaces, indent=cfg.indent, 
                                                                has_copy_constructor=has_copy_constructor, has_assignment_operator=has_assignment_operator)
            class_decl += '\n'

        return_code_dict[abstr_class_fname].append( (-1, class_decl) )


        # Add abstract class to inheritance list of original class

        line_number = int(class_el.get('line'))

        # Get file content, but replace all comments with whitespace
        f = open(src_file_name, 'r')
        file_content = f.read()
        f.close()
        file_content_nocomments = utils.removeComments(file_content, insert_blanks=True)

        # Find index of the \n in line number line_number
        newline_pos = utils.findNewLinePos(file_content_nocomments, line_number)

        # First, find position of class name
        search_limit = newline_pos
        while search_limit > -1:
            pos = file_content_nocomments[:search_limit].rfind(short_class_name)
            pre_char  = file_content_nocomments[pos-1]
            post_char = file_content_nocomments[pos+len(short_class_name)]
            if (pre_char in [' ','\n','\t']) and (post_char in [' ', ':', '\n', '<', '{']):
                break
            else:
                search_limit = pos

        class_name_pos = pos


        # Special preparations for template classes:
        if is_template:
        
            # - Determine whether this is the source for the general template 
            #   or for a specialization (look for '<' after class name)
            temp_pos = class_name_pos + len(short_class_name)
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
        if (class_el.get('bases') == "") and (class_name_full not in added_parent):

            # - Calculate insert position
            insert_pos = class_name_pos + len(short_class_name)
            if is_template and src_is_specialization:
                insert_pos += len(add_template_bracket)

            # - Generate code
            add_code = ' : public virtual ' + short_abstr_class_name
            if is_template == True:
                add_code += add_template_bracket

        # If there are previous parent classes
        else:

            # - Get colon position
            if is_template and src_is_specialization:
                temp_pos = class_name_pos + len(short_class_name) + len(add_template_bracket)
            else:
                temp_pos = class_name_pos + len(short_class_name)
            colon_pos = temp_pos + file_content_nocomments[temp_pos:newline_pos].find(':')

            # - Calculate insert position
            insert_pos = colon_pos + 1

            # - Generate code
            add_code = ' public virtual ' + short_abstr_class_name
            if is_template == True:
                add_code += add_template_bracket
            add_code += ','

        # - Register new code
        if src_file_name not in return_code_dict.keys():
            return_code_dict[src_file_name] = []
        return_code_dict[src_file_name].append( (insert_pos, add_code) )

        # - Update added_parent dict
        added_parent.append(class_name_full)


        # Generate code for #include statement in orginal header/source file 
        src_include_line = '#include "' + os.path.join(cfg.add_path_to_includes, short_abstr_class_fname) + '"'
        if src_include_line in src_includes[src_file_name]:
            pass
        else:
            # - Find position
            if is_template == True:
                insert_pos = file_content_nocomments[:class_name_pos].rfind('template')
            else:
                insert_pos = file_content_nocomments[:class_name_pos].rfind('class')
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
            include_code += use_indent + src_include_line + '\n'
            include_code += use_indent
            for ns in namespaces:
                include_code += 'namespace ' + ns + ' { '
            include_code += '\n'*has_namespace

            # - Register code
            return_code_dict[src_file_name].append( (insert_pos, include_code) )

            # - Register include line
            src_includes[src_file_name].append(src_include_line)


        # Generate wrappers for all member functions that make use of native types

        # - Create lists of all 'non-artificial' members of the class
        member_methods   = []
        member_variables = []
        if 'members' in class_el.keys():
            for mem_id in class_el.get('members').split():
                el = cfg.id_dict[mem_id]
                if not 'artificial' in el.keys():
                    if (el.tag == 'Method') and (not funcutils.ignoreFunction(el)):
                        if funcutils.usesNativeType(el):
                            member_methods.append(el)
                    elif (el.tag in ('Field', 'Variable')) and (el.get('access') == 'public'):
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

        # - Generate code for each member method
        wrapper_code  = '\n'
        current_access = None
        for method_el in member_methods:
            method_access = method_el.get('access')
            if method_access != current_access:
                wrapper_code += ' '*(len(namespaces)+1)*cfg.indent + method_access +':\n'
                current_access = method_access
            wrapper_code += classutils.constrWrapperFunction(method_el, indent=cfg.indent, n_indents=len(namespaces)+2)
            # wrapper_code += ' '*(len(namespaces)+2)*cfg.indent + 'WRAPPER CODE FOR ' + mem_el.get('name') + ' GOES HERE!\n'

        # - Register code
        return_code_dict[src_file_name].append( (insert_pos, wrapper_code) )            

        # - Generate a reference-returning method for each (public) member variable:
        ref_func_code = ''
        if len(member_variables) > 0:
            n_indents = len(namespaces)
            ref_func_code += '\n'
            ref_func_code += ' '*cfg.indent*(n_indents+1) + 'public:\n'
            for var_el in member_variables:
                ref_func_code += classutils.constrVariableRefFunction(var_el, virtual=False, indent=cfg.indent, n_indents=n_indents+2)
                ref_func_code += '\n'

        # - Register code
        if ref_func_code != '':
            return_code_dict[src_file_name].append( (insert_pos, ref_func_code) )            


        # Generate pointer-based copy and assignment functions
        n_indents = len(namespaces)
        ptr_code  = '\n'
        ptr_code += ' '*cfg.indent*(n_indents+1) + 'public:\n'

        ptr_code += classutils.constrPtrCopyFunc(short_abstr_class_name, short_class_name, virtual=False, indent=cfg.indent, n_indents=n_indents+2)
        ptr_code += classutils.constrPtrAssignFunc(short_abstr_class_name, short_class_name, virtual=False, indent=cfg.indent, n_indents=n_indents+2)
        
        # - Register code
        return_code_dict[src_file_name].append( (insert_pos, ptr_code) )


        # Generate info for factory 
        
        # fact_subdict = {}
        # fact_subdict['include']  = '#include "' + os.path.basename(class_file_name) + '"\n'
        # fact_subdict['func_def'] = constrFactoryFunction(class_el, indent=cfg.indent)
        # factories_dict[short_class_name] = fact_subdict

        factory_file_content  = ''
        if is_template and class_name in template_done:
            pass
        else:
            src_file_name_base = os.path.basename(src_file_name)
            factory_file_content += '#include "' + os.path.join(cfg.add_path_to_includes, src_file_name_base) + '"\n'
            factory_file_content += '\n'
        if is_template:
            factory_file_content += classutils.constrFactoryFunction(class_el, class_name_full, indent=cfg.indent, template_types=spec_template_types)
        else:
            factory_file_content += classutils.constrFactoryFunction(class_el, class_name_full, indent=cfg.indent)
        factory_file_content += '\n'

        # # Add include statements
        # include_list = []
        # for fact_subdict in factories_dict.values():
        #     inc_line = fact_subdict['include']
        #     if inc_line not in include_list:
        #         include_list.append(inc_line)
        # factory_file_content += ''.join(include_list)

        # factory_file_content += 2*'\n'

        # # Add function definitions
        # for fact_subdict in factories_dict.values():
        #     factory_file_content += fact_subdict['func_def']

        # - Generate factory file name
        dir_name = cfg.extra_output_dir
        # dir_name = os.path.split(src_file_el.get('name'))[0]
        factory_file_name = os.path.join(dir_name, cfg.factory_file_prefix + short_class_name + cfg.source_extension)

        # - Register code
        if factory_file_name not in return_code_dict.keys():
            return_code_dict[factory_file_name] = []
        return_code_dict[factory_file_name].append( (-1, factory_file_content) )


        #
        # Generate a header containing the GAMBIT wrapper class
        #
        short_wrapper_class_fname = cfg.wrapper_header_prefix + short_class_name + cfg.header_extension
        wrapper_class_fname = os.path.join(cfg.extra_output_dir, short_wrapper_class_fname)
        
        wrapper_class_file_content = generateWrapperHeader(class_el, short_class_name, short_abstr_class_name, namespaces, 
                                                           short_abstr_class_fname, has_copy_constructor, has_assignment_operator)

        # - Update dict of wrapper class header files
        #wrapper_class_headers[short_class_name] = wrapper_class_fname

        # - Register code
        if wrapper_class_fname not in return_code_dict.keys():
            return_code_dict[wrapper_class_fname] = []
        return_code_dict[wrapper_class_fname].append( (0, wrapper_class_file_content) )


        #
        # Keep track of classes done
        #
        classes_done.append(class_name_full)
        if is_template: 
            if class_name not in template_done:
                template_done.append(class_name)
            if class_name_full not in templ_spec_done:
                templ_spec_done.append(class_name_full)
        
        print '...Done'
        print
        print

    #
    # END: Loop over all classes in cfg.class_dict
    #    

    print
    print 'CLASSES DONE:    ', classes_done
    print 'TEMPLATE DONE:   ', template_done
    print 'TEMPL_SPEC DONE: ', templ_spec_done
    print


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
def generateWrapperHeader(class_el, short_class_name, short_abstr_class_name, namespaces, 
                          short_abstr_class_fname, has_copy_constructor, has_assignment_operator):

    # Useful variables
    abstr_class_name = '::'.join(namespaces) + '::'*bool(len(namespaces)) + short_abstr_class_name
    indent           = ' '*cfg.indent

    # Useful lists
    class_variables    = []
    class_functions    = []
    class_constructors = []
    class_members      = utils.getMemberElements(class_el, include_artificial=False)
    class_members_full = utils.getMemberElements(class_el, include_artificial=True)
    for mem_el in class_members:
        if (mem_el.tag == 'Method') and (mem_el.get('access') == 'public'):
            if funcutils.ignoreFunction(mem_el, limit_pointerness=True):
                warnings.warn('The member "%s" in class "%s" makes use of a non-accepted type and will be ignored.' % (mem_el.get('name'), short_class_name))
            else:
                class_functions.append(mem_el)

        elif (mem_el.tag in ('Field', 'Variable')) and (mem_el.get('access') == 'public'):
            if utils.isAcceptedType(mem_el):
                class_variables.append(mem_el)
            else:
                warnings.warn('The member "%s" in class "%s" makes use of a non-accepted type and will be ignored.' % (mem_el.get('name'), short_class_name))
        else:
            pass
    for mem_el in class_members_full:
        if (mem_el.tag == 'Constructor') and (mem_el.get('access') == 'public'):
            if funcutils.ignoreFunction(mem_el, limit_pointerness=True):
                warnings.warn('The member "%s" in class "%s" makes use of a non-accepted type and will be ignored.' % (mem_el.get('name'), short_class_name))
            else:
                class_constructors.append(mem_el)


    #
    # Start code generation
    #

    code = ''

    # FOR DEBUG: print all types used in this class
    code += '// FOR DEBUG:\n'
    code += '// Types originally used in this class:\n'
    for type_name in cfg.all_types_in_class.keys():
        code += '// - ' + type_name + '\n'

    # # Add include statement for abstract class header
    # code += '\n'
    # code += '#include "' + short_abstr_class_fname + '"\n'

    # # Figure out what types can be used and generate the necessary include statements
    # ignore_types = []
    # for type_name, type_el in cfg.all_types_in_class.items():
    #     basic_type_name = type_name.replace('*','').replace('&','')

    #     if utils.isFundamental(type_el):
    #         continue

    #     elif utils.isNative(type_el):
    #         if basic_type_name in cfg.loaded_classes:
    #             basic_type_name_short = basic_type_name.split('::')[-1]
    #             wrapper_header = cfg.wrapper_header_prefix + basic_type_name_short + cfg.header_extension
    #             code += '#include "' + wrapper_header + '"\n'            
    #         else:
    #             warnings.warn('The type "%s" is indetified as native to the source code, yet it is not registred as an accepted type. All class members using this type will be ignored.' % (basic_type_name))
    #             ignore_types.append(type_el)

    #     elif utils.isStdType(type_el):
    #             warnings.warn('The type "%s" is indetified as a std type. We need to identify the correct header to inlcude...' % (basic_type_name))

    #     else:
    #             warnings.warn('The type "%s" is unknown. All class members using this type will be ignored.' % (basic_type_name))
    #             ignore_types.append(type_el)


    # Create one factory function pointer for each constructor.
    # Initialize pointers to NULL - they will later be filled by dynamic loading
    # We append a number (e.g. '_1') to the function pointer name to make sure each overload of a constructor
    # gets a unique.
    temp_code = ''
    for i, constr_el in enumerate(class_constructors):

        # Identify arguments, translate argument type of loaded classes
        # and construct the argument bracket
        args = funcutils.getArgs(constr_el)
        w_args = funcutils.constrWrapperArgs(args, add_ref=True)
        args_bracket = funcutils.constrArgsBracket(w_args, include_arg_name=False, include_arg_type=True)

        # Construct factory pointer code
        temp_code += abstr_class_name + '* (*Factory_' + short_class_name + '_' + str(i) + ')' + args_bracket + ' = NULL;\n'

    if temp_code != '':
        code += '\n'
        code += "// Factory function pointers to be filled by dynamic loading\n"
        code += temp_code


    #
    # Generate code for wrapper class
    #

    short_wrapper_class_name = short_class_name + cfg.code_suffix

    # Class declaration line
    code += '\n'
    code += 'class ' + short_wrapper_class_name + '\n'

    # Class body
    code += '{\n'

    # Start private block
    code += indent + 'private:\n'
    
    # Private variable: bool member_element
    code += 2*indent + 'bool member_variable;\n'

    # Start public block
    code += indent + 'public:\n'

    # Variables:
    code += 2*indent + '// Member variables: \n'

    # Public variable: pointer to abstract class: Abstract_ClassX* BEptr
    code += 2*indent + abstr_class_name + '* BEptr;\n'

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

    # Add wrappers for all member functions
    for func_el in class_functions:

        # Function name
        func_name = func_el.get('name')

        # Determine return type
        return_type, return_type_kw, return_type_id = utils.findType(func_el)
        return_type_el = cfg.id_dict[return_type_id]

        # Return keywords
        return_kw_str = ' '.join(return_type_kw)        
        return_kw_str += ' '*bool(len(return_type_kw))

        # Arguments
        args = funcutils.getArgs(func_el)
        args_bracket = funcutils.constrArgsBracket(args, include_arg_name=True, include_arg_type=True, include_namespace=True, use_wrapper_class=True)

        # Convert return type if loaded class
        if utils.isLoadedClass(return_type, byname=True):
            use_return_type = classutils.toWrapperType(return_type)
        else:
            use_return_type = return_type

        # Write declaration line
        code += 2*indent + return_kw_str + use_return_type + ' ' + func_name + args_bracket + '\n'

        # Write function body
        code += 2*indent + '{\n'

        if return_type == 'void':
            code += 3*indent
        else:
            code += 3*indent + 'return '

        args_bracket_notypes = funcutils.constrArgsBracket(args, include_arg_name=True, include_arg_type=False, wrapper_to_pointer=True)
        code += 'BEptr->' + func_name + args_bracket_notypes + ';\n'

        code += 2*indent + '}\n'


    # Add special member function: _set_member(bool) - set the variable 'member_variable' to true/false
    code += '\n'
    code += 2*indent + '// Special member function to set member_variable: \n'
    code += 2*indent + 'void _set_member_variable(bool in) { member_variable = in; }\n'    


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
    common_init_list_code += 3*indent + 'member_variable(false)\n'

    common_constructor_body = ''
    common_constructor_body += 2*indent + '{\n'
    for var_el in class_variables:
        var_name = var_el.get('name')
        var_type, var_type_kw, var_type_id = utils.findType(var_el)
        if utils.isLoadedClass(var_type, byname=True):
            common_constructor_body += 3*indent + var_name + '._set_member_variable(true);\n'
    common_constructor_body += 2*indent + '}\n'


    # Add wrappers for all original constructors except the copy constructor
    temp_code = ''
    for i, constr_el in enumerate(class_constructors):

        # Identify arguments
        args = funcutils.getArgs(constr_el)

        # Skip if this is a copy constructor. (Another copy constructor is added below.)
        if (len(args) == 1) and (args[0]['id'] == class_el.get('id')):
            continue

        args_bracket = funcutils.constrArgsBracket(args, include_arg_name=True, include_arg_type=True, include_namespace=True, use_wrapper_class=True)
        args_bracket_notypes = funcutils.constrArgsBracket(args, include_arg_name=True, include_arg_type=False, wrapper_to_pointer=True)

        temp_code += 2*indent + short_wrapper_class_name + args_bracket + ' :\n'
        temp_code += 3*indent + 'BEptr( Factory_' + short_class_name + '_' + str(i) + args_bracket_notypes + ' ),\n'  # FIXME: This is not general. Fix argument list.
        temp_code += common_init_list_code
        temp_code += common_constructor_body

        temp_code += '\n'

    if temp_code != '':
        code += '\n'
        code += 2*indent + '// Wrappers for original constructors: \n'    
        code += temp_code


    # Add special constructor based on abstract pointer (This is needed to allow return-by-value with the wrapper classes.)
    code += 2*indent + '// Special pointer-based constructor: \n'
    code += 2*indent + short_wrapper_class_name + '(' + abstr_class_name +'* in) :\n'
    code += 3*indent + 'BEptr(in),\n'
    code += common_init_list_code
    code += common_constructor_body


    # Add copy constructor
    if has_copy_constructor:
        code += '\n'
        code += 2*indent + '// Copy constructor: \n'
        code += 2*indent + short_wrapper_class_name + '(const ' + short_wrapper_class_name +'& in) :\n'
        code += 3*indent + 'BEptr(in.BEptr->pointerCopy' + cfg.code_suffix + '()),\n'
        code += common_init_list_code
        code += common_constructor_body

    #
    # Add assignment operator
    #
    if has_assignment_operator:
        code += '\n'
        code += 2*indent + '// Assignment operator: \n'
        code += 2*indent + short_wrapper_class_name + '& ' + short_wrapper_class_name + '::operator=(const ' + short_wrapper_class_name +'& in)\n'
        code += 2*indent + '{\n'
        code += 3*indent + 'if (this != &in) { BEptr->pointerAssign' + cfg.code_suffix + '(in.BEptr); }\n'
        code += 2*indent + '}\n'

    #
    # Add destructor
    #
    code += '\n'
    code += 2*indent + '// Destructor: \n'
    code += 2*indent + '~' + short_wrapper_class_name + '()\n'
    code += 2*indent + '{\n'
    code += 3*indent + 'if(member_variable==false) { delete BEptr; }\n'
    code += 2*indent + '}\n'

    # Close class body
    code += '};\n'

    # Add include guards
    guard_var = '__' + cfg.wrapper_header_prefix + short_class_name + '_' + cfg.header_extension.replace('.','') + '__'
    guard_var = guard_var.upper()
    guard_code_top = '#ifndef ' + guard_var + '\n' + '#define ' + guard_var + '\n'
    guard_code_bottom = '#endif /* ' + guard_var + ' */\n'
    
    code = guard_code_top + '\n' + code + '\n' + guard_code_bottom

    # Return code
    return code

#
# END: generateWrapperHeader
#   
