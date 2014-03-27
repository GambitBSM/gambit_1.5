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

#
# Module-level globals
#

classes_done  = []
template_done = []
added_parent  = []

#
# Main function for parsing classes
#

def run():

    # Prepare returned dict
    return_code_dict = OrderedDict()

    #
    # Loop over all classes 
    #
    
    for class_name_full, class_el in cfg.class_dict.items():


        # # Print current class
        # print 'Current class: ' + class_name_full


        # Check if we've done this class already
        if class_name_full in classes_done:
            print ' '*15 + '--> Class already done'
            continue


        # Check if this class is native to the source code
        if utils.isNative(class_el):
            if 'extern' in class_el.keys() and class_el.get('extern') == "1":
                # print ' '*15 + '--> Class %s RECOGNIZED but IGNORED' % class_name_full
                continue
            else:
                # print ' '*15 + '--> Class %s ACCEPTED' % class_name_full
                pass
        else:
            # print ' '*15 + '--> Class IGNORED'
            continue

        # Print current class
        print 'Current class: ' + class_name_full


        # Check if this is a template class
        if '<' in class_name_full:
            is_template = True

            # HACK
            # runTemplate(class_el, class_name_full)

            class_name  = class_name_full.split('<',1)[0]
        else:
            is_template = False
            class_name = class_name_full

        
        # Determine a set of useful variables
        src_file_el   = cfg.id_dict[class_el.get('file')]
        src_file_name = src_file_el.get('name')
        src_dir       = os.path.split(src_file_name)[0]
        short_class_name = class_el.get('name').split('<',1)[0]


        # If template class, figure out template variables
        if is_template == True:
            if not class_name in template_done:
                template_bracket, template_types = utils.getTemplateBracket(class_el)
            spec_template_types = utils.getSpecTemplateTypes(class_el)
            print ' --> template bracket, types, spec. types: ', template_bracket, template_types, spec_template_types


        # Define name and path for new header file 'abstract_CLASS.hpp'
        
        abstr_class_fname = utils.generateAbstractFilePath(class_el, short_class_name, prefix=cfg.abstr_header_prefix, suffix='.hpp')
        short_abstr_class_fname = os.path.split(abstr_class_fname)[-1]

        # short_abstr_class_fname  = cfg.abstr_header_prefix + short_class_name.lower() + '.hpp'
        # src_file_el = cfg.id_dict[class_el.get('file')]
        # dir_name = os.path.split(src_file_el.get('name'))[0]
        # abstr_class_fname = os.path.join(dir_name,short_abstr_class_fname)


        # Register an include statement for this header file, to go in the file cfg.all_headers_fname
        include_line = '#include "' + short_abstr_class_fname + '"\n' 
        if cfg.all_headers_fname not in return_code_dict.keys():
            return_code_dict[cfg.all_headers_fname] = []
        return_code_dict[cfg.all_headers_fname].append( (0, include_line) )


        # Construct abstract class declaration and register the new code
        class_decl = ''
        if (is_template == True) and (class_name not in template_done):
            class_decl += classutils.constructEmptyTemplClassDecl(class_name_full, template_bracket, indent=cfg.indent)
            class_decl += classutils.constructAbstractClassDecl(class_el, indent=cfg.indent, template_bracket=template_bracket)
        else:
            class_decl += classutils.constructAbstractClassDecl(class_el, indent=cfg.indent)
            class_decl += '\n'


        if abstr_class_fname not in return_code_dict.keys():
            return_code_dict[abstr_class_fname] = []
        return_code_dict[abstr_class_fname].append( (0, class_decl) )


        # Add abstract class to inheritance list of original class
        if (is_template == True) and (class_name in template_done):
            pass
        else:
            line_number   = int(class_el.get('line'))
            short_abstr_class_name = classutils.getAbstractClassName(class_name, short=True)

            # Get file content, but replace all comments with whitespace
            f = open(src_file_name, 'r+')
            file_content = f.read()
            f.close()
            file_content_nocomments = utils.removeComments(file_content, insert_blanks=True)

            # Find index of the \n in line number line_number
            newline_pos = utils.findNewLinePos(file_content_nocomments, line_number)

            # Determine position for code insertion

            # - First, find position of class name
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

            # - If no previous parent classes:
            if (class_el.get('bases') == "") and (class_name not in added_parent):

                # - Calculate insert position
                insert_pos = class_name_pos + len(short_class_name)

                # - Generate code
                add_code = ' : public virtual ' + short_abstr_class_name
                if is_template == True:
                    add_code += '<' + ','.join(template_types) + '>'

            # - If there are previous parent classes
            else:

                # - Get colon position
                colon_pos = class_name_pos + file_content_nocomments[class_name_pos:newline_pos].find(':')

                # - Calculate insert position
                insert_pos = colon_pos + 1

                # - Generate code
                add_code = ' public virtual ' + short_abstr_class_name
                if is_template == True:
                    add_code += '<' + ','.join(template_types) + '>'
                add_code += ','

            # - Register new code
            if src_file_name not in return_code_dict.keys():
                return_code_dict[src_file_name] = []
            return_code_dict[src_file_name].append( (insert_pos, add_code) )

            # - Update added_parent dict
            added_parent.append(class_name)


            # Generate code for #include statement in orginal header/source file 
            if src_file_name not in return_code_dict.keys():
                return_code_dict[src_file_name] = []

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
            namespaces = class_name.split('::')[:-1]
            has_namespace = bool(len(namespaces))

            include_code += use_indent
            for ns in namespaces:
                include_code += '} '
            include_code += '\n'*has_namespace
            include_code += use_indent + '#include "' + short_abstr_class_fname + '"\n'
            include_code += use_indent
            for ns in namespaces:
                include_code += 'namespace ' + ns + ' { '
            include_code += '\n'*has_namespace

            # - Register code
            if src_file_name not in return_code_dict.keys():
                return_code_dict[src_file_name] = []
            return_code_dict[src_file_name].append( (insert_pos, include_code) )


        # Generate info for factory 
        
        # fact_subdict = {}
        # fact_subdict['include']  = '#include "' + os.path.basename(class_file_name) + '"\n'
        # fact_subdict['func_def'] = constructFactoryFunction(class_el, indent=cfg.indent)
        # factories_dict[short_class_name] = fact_subdict

        factory_file_content  = ''
        factory_file_content += '#include "' + os.path.basename(src_file_name) + '"\n'
        factory_file_content += '\n'
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
        dir_name = os.path.split(src_file_el.get('name'))[0]
        factory_file_name = os.path.join(dir_name, cfg.factory_file_prefix + short_class_name + '.cpp')

        # - Register code
        if factory_file_name not in return_code_dict.keys():
            return_code_dict[factory_file_name] = []
        return_code_dict[factory_file_name].append( (0, factory_file_content) )


        # Keep track of classes done
        classes_done.append(class_name_full)
        if (is_template == True) and (class_name not in template_done):
            template_done.append(class_name)

        
        print '...Done'
        print
        print

    print
    print 'CLASSES DONE:  ', classes_done
    print 'TEMPLATE DONE: ', template_done
    print

    #
    # Return result
    #

    return return_code_dict







#
# Function for parsing template classes
#

def runTemplate(class_el):

    # Get class name without template bracket
    class_name  = class_name_full.split('<',1)[0]


    # Figure out template variables
    if class_name not in template_done:
        template_bracket, template_types = utils.getTemplateBracket(class_el)
    spec_template_types = utils.getSpecTemplateTypes(class_el)
    print ' --> template bracket, types, spec. types: ', template_bracket, template_types, spec_template_types


    # Add class_name to list of registred template classes
    if class_name not in template_done:
        template_done.append(class_name)



    # Return result
    return return_code_dict









