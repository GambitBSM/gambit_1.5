############################
#                          #
#  Class parsing for BOSS  #
#                          #
############################

import xml.etree.ElementTree as ET
from collections import OrderedDict
import os

import modules.cfg as cfg
import modules.utils as utils
import modules.classutils as classutils
import modules.funcutils as funcutils

#
# Module-level globals
#

classes_done    = []
template_done   = []
templ_spec_done = []
added_parent    = []

class_name             = None
short_class_name       = None
short_abstr_class_name = None

src_file_el   = None
src_file_name = None
src_dir       = None

short_abstr_class_fname = None
abstr_class_fname       = None

namespaces    = None
has_namespace = None


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

        class_id = class_el.get('id')
        for mem_el in check_member_elements:

            if mem_el.tag in ['Constructor', 'Destructor', 'Method', 'OperatorMethod']:
                args_list = funcutils.getArgs(mem_el)
                for arg_dict in args_list:
                    arg_type_el = cfg.id_dict[arg_dict['id']]

                    arg_type = arg_dict['type']

                    if arg_type not in cfg.all_types_in_class:
                        cfg.all_types_in_class.append(arg_type)

                        if utils.isStdType(arg_type_el):
                            cfg.std_types_used.append(arg_type) 


            if ('type' in mem_el.keys()) or ('returns' in mem_el.keys()):
                mem_type, mem_type_kw, mem_type_id = utils.findType(mem_el)
                type_el = cfg.id_dict[mem_type_id]

                if mem_type not in cfg.all_types_in_class:
                    cfg.all_types_in_class.append(mem_type)

                    if utils.isStdType(type_el):
                        cfg.std_types_used.append(mem_type) 


                # if utils.isStdType(type_el):
                #     if mem_type not in cfg.std_types_used:
                #         cfg.std_types_used.append(mem_type)

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


        
        # Set the global module-level variables 
        class_name             = class_name_full.split('<',1)[0]
        short_class_name       = class_el.get('name').split('<',1)[0]
        short_abstr_class_name = classutils.getAbstractClassName(class_name, prefix=cfg.abstr_class_prefix, short=True)

        src_file_el        = cfg.id_dict[class_el.get('file')]
        src_file_name      = src_file_el.get('name')
        src_dir            = os.path.split(src_file_name)[0]
        new_src_file_name  = os.path.join(cfg.extra_output_dir, os.path.basename(src_file_name))

        short_abstr_class_fname = cfg.abstr_header_prefix + short_class_name + cfg.header_extension
        abstr_class_fname       = os.path.join(cfg.extra_output_dir, short_abstr_class_fname)

        namespaces    = class_name.split('::')[:-1]
        has_namespace = bool(len(namespaces))


        # Prepare entries in return_code_dict and src_includes
        if abstr_class_fname not in return_code_dict.keys():
            return_code_dict[abstr_class_fname] = []
        if new_src_file_name not in return_code_dict.keys():
            return_code_dict[new_src_file_name] = []
        if new_src_file_name not in src_includes.keys():
            src_includes[new_src_file_name] = []


        # Register an include statement for this header file, to go in the file cfg.all_headers_fname
        include_line = '#include "' + short_abstr_class_fname.lower() + '"\n' 
        all_headers_file_path = os.path.join(cfg.extra_output_dir, cfg.all_headers_fname)

        if all_headers_file_path not in return_code_dict.keys():
            return_code_dict[all_headers_file_path] = []
        if include_line not in [ inc_tuple[1] for inc_tuple in return_code_dict[all_headers_file_path] ]:
            return_code_dict[all_headers_file_path].append( (0, include_line) )


        # Treat the first specialization of a template class differently
        if is_template and class_name not in template_done:
            template_bracket, template_types = utils.getTemplateBracket(class_el)
            
            empty_templ_class_decl = ''
            empty_templ_class_decl += classutils.constructEmptyTemplClassDecl(short_abstr_class_name, namespaces, template_bracket, indent=cfg.indent)
            empty_templ_class_decl += classutils.constructTemplForwDecl(short_class_name, namespaces, template_bracket, indent=cfg.indent)

            return_code_dict[abstr_class_fname].append( (0, empty_templ_class_decl) )


        # Get template arguments for specialization, 
        # and check that they are acceptable
        if is_template and class_name_full not in templ_spec_done:
            spec_template_types = utils.getSpecTemplateTypes(class_el)
            for template_type in spec_template_types:
                if template_type not in cfg.accepted_types:
                    raise Exception("The template specialization type '" + template_type + "' for class " + class_name_full + " is not among accepted types.")


        # Construct code for the abstract class header file and register it

        class_decl = ''

        # - Add forward declarations at the top of the header file  
        #
        #   FIXME: For now, this includes forward declarations of *all* accepted native classes.
        #          Should limit this to only the forward declarations needed for the given abstract class.
        class_decl += '// Forward declarations:\n'
        class_decl += utils.constrForwardDecls()
        class_decl += '\n'


        # - Add the the code for the abstract class
        #
        # if (is_template == True) and (class_name_full in templ_spec_done):
        #     pass
        # elif (is_template == True) and (class_name_full not in templ_spec_done):
        #     class_decl += classutils.constructAbstractClassDecl(class_el, short_class_name, short_abstr_class_name, namespaces, indent=cfg.indent, template_bracket=add_template_bracket)
        #     class_decl += '\n'
        # else:
        #     class_decl += classutils.constructAbstractClassDecl(class_el, short_class_name, short_abstr_class_name, namespaces, indent=cfg.indent)
        #     class_decl += '\n'

        if (is_template == True) and (class_name_full in templ_spec_done):
            pass
        elif (is_template == True) and (class_name_full not in templ_spec_done):
            class_decl += classutils.constructAbstractClassDecl(class_el, short_class_name, short_abstr_class_name, namespaces, indent=cfg.indent, template_types=spec_template_types)
            class_decl += '\n'
        else:
            class_decl += classutils.constructAbstractClassDecl(class_el, short_class_name, short_abstr_class_name, namespaces, indent=cfg.indent)
            class_decl += '\n'

        return_code_dict[abstr_class_fname].append( (-1, class_decl) )


        # Add abstract class to inheritance list of original class

        # if (is_template == True) and (class_name_full in templ_spec_done):
        #     pass
        #     print
        #     print 'IT DID HAPPEN!'
        #     print '---> THIS IS A TEMPLATE SPECIALIZATION THAT HAS BEEN DONE: ', class_name_full
        #     print
        # else:

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
        if new_src_file_name not in return_code_dict.keys():
            return_code_dict[new_src_file_name] = []
        return_code_dict[new_src_file_name].append( (insert_pos, add_code) )

        # - Update added_parent dict
        added_parent.append(class_name_full)


        # Generate code for #include statement in orginal header/source file 
        src_include_line = '#include "' + os.path.join(cfg.add_path_to_includes, short_abstr_class_fname) + '"'
        if src_include_line in src_includes[new_src_file_name]:
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
            return_code_dict[new_src_file_name].append( (insert_pos, include_code) )

            # - Register include line
            src_includes[new_src_file_name].append(src_include_line)


        # # Generate converter functions ('upcast', 'abstractify') for the original class

        # # - Determine insert position  -- DUPLICATED BELOW
        # rel_pos_start, rel_pos_end = utils.getBracketPositions(file_content_nocomments[class_name_pos:], delims=['{','}'])
        # class_body_start = class_name_pos + rel_pos_start
        # class_body_end   = class_name_pos + rel_pos_end

        # # - Generate code
        # converter_code  = '\n'
        # converter_code += ' '*cfg.indent*(len(namespaces)+1) + 'public:\n'
        # if is_template:
        #     converter_code += classutils.constructUpcastFunction(short_abstr_class_name, indent=cfg.indent, n_indents=len(namespaces)+2, template_bracket=add_template_bracket)
        #     converter_code += '\n'
        #     converter_code += classutils.constructAbstractifyFunction(short_abstr_class_name, indent=cfg.indent, n_indents=len(namespaces)+2, template_bracket=add_template_bracket)
        # else:
        #     converter_code += classutils.constructUpcastFunction(short_abstr_class_name, indent=cfg.indent, n_indents=len(namespaces)+2)
        #     converter_code += '\n'
        #     converter_code += classutils.constructAbstractifyFunction(short_abstr_class_name, indent=cfg.indent, n_indents=len(namespaces)+2)
        # converter_code += ' '*cfg.indent*len(namespaces)

        # # - Register code
        # insert_pos = class_body_end
        # return_code_dict[src_file_name].append( (insert_pos, converter_code) )


        # Generate wrappers for all member functions

        # - Determine insert position
        rel_pos_start, rel_pos_end = utils.getBracketPositions(file_content_nocomments[class_name_pos:], delims=['{','}'])
        class_body_start = class_name_pos + rel_pos_start
        class_body_end   = class_name_pos + rel_pos_end

        # - Create list of all 'non-artificial' members of the class
        member_methods = []
        if 'members' in class_el.keys():
            for mem_id in class_el.get('members').split():
                el = cfg.id_dict[mem_id]
                if (el.tag == 'Method') and (not 'artificial' in el.keys()) and (not funcutils.ignoreFunction(el)):
                    member_methods.append(el)

        # - Determine insert position
        insert_pos = class_body_end

        # - Generate code for each member
        wrapper_code  = '\n'
        current_access = None
        for method_el in member_methods:
            method_access = method_el.get('access')
            if method_access != current_access:
                wrapper_code += ' '*(len(namespaces)+1)*cfg.indent + method_access +':\n'
                current_access = method_access
            wrapper_code += classutils.constructWrapperFunction(method_el, indent=cfg.indent, n_indents=len(namespaces)+2)
            # wrapper_code += ' '*(len(namespaces)+2)*cfg.indent + 'WRAPPER CODE FOR ' + mem_el.get('name') + ' GOES HERE!\n'

        # - Register code
        return_code_dict[new_src_file_name].append( (insert_pos, wrapper_code) )            


        # Generate info for factory 
        
        # fact_subdict = {}
        # fact_subdict['include']  = '#include "' + os.path.basename(class_file_name) + '"\n'
        # fact_subdict['func_def'] = constructFactoryFunction(class_el, indent=cfg.indent)
        # factories_dict[short_class_name] = fact_subdict

        factory_file_content  = ''
        if is_template and class_name in template_done:
            pass
        else:
            src_file_name_base = os.path.basename(src_file_name)
            factory_file_content += '#include "' + os.path.join(cfg.add_path_to_includes, src_file_name_base) + '"\n'
            factory_file_content += '\n'
        if is_template:
            factory_file_content += classutils.constructFactoryFunction(class_el, class_name_full, indent=cfg.indent, template_types=spec_template_types)
        else:
            factory_file_content += classutils.constructFactoryFunction(class_el, class_name_full, indent=cfg.indent)
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


        # Keep track of classes done
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







# #
# # Function for parsing template classes
# #

# def runTemplate(class_el):

#     # Figure out template variables
#     if class_name not in template_done:
#         template_bracket, template_types = utils.getTemplateBracket(class_el)
#     spec_template_types = utils.getSpecTemplateTypes(class_el)
#     print ' --> template bracket, types, spec. types: ', template_bracket, template_types, spec_template_types


#     # Add class_name to list of registred template classes
#     if class_name not in template_done:
#         template_done.append(class_name)



#     # Return result
#     return return_code_dict









