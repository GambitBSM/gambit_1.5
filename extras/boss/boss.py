#######################################################
#                                                     #
#  First version of BOSS - Backend-On-a-Stick Script  #
#                                                     #
#######################################################

#
# This is a Python package for making the classes of a C++ library 
# available for dynamic use through the 'dlopen' system.
#
# To run BOSS, the library source code must first be parsed and analyzed
# using 'gccxml'. The resulting XML files is the starting point for BOSS.
# 

import xml.etree.ElementTree as ET
import os
import sys
import warnings
from collections import OrderedDict
from optparse import OptionParser

import modules.cfg as cfg
import modules.classutils as classutils
import modules.classparse as classparse
import modules.funcparse as funcparse
import modules.funcutils as funcutils
import modules.utils as utils


# ====== main ========

def main():

    #
    # Initialization
    #

    # Prepare container variables
    new_code = OrderedDict()

    # Create the extra output directory if it does not exist
    try:
      os.mkdir(cfg.extra_output_dir)
    except OSError, e:
      if e.errno == 17:
        pass
      else:
        raise

    # Parse command line arguments and options
    parser = OptionParser(usage="usage: %prog [options] xmlfiles",
                          version="%prog 0.1")
    parser.add_option("-l", "--list",
                      action="store_true",
                      dest="list_classes_flag",
                      default=False,
                      help="output a list of the available classes and functions")
    parser.add_option("-d", "--debug-mode",
                      action="store_true",
                      dest="debug_mode_flag",
                      default=False,
                      help="Test run. No source files are changed.")
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("No xml file(s) specified")

    # Get the xml file names from command line input
    xml_files = args


    #
    # Loop over all the xml files
    #

    for current_xml_file in xml_files:

        # Output xml file name
        print
        print '=====:  Current XML file: %s  :=====' % current_xml_file
        print

        # Set global xml file name and parse it using ElementTree
        cfg.xml_file_name = current_xml_file
        tree = ET.parse(cfg.xml_file_name)
        root = tree.getroot()


        # Update global dict: id --> xml element (all elements)
        cfg.id_dict = OrderedDict([ (el.get('id'), el) for el in root.getchildren() ]) 


        # Update global dict: file name --> file xml element
        cfg.file_dict = OrderedDict([ (el.get('name'), el) for el in root.findall('File') ])


        # Update global dict: std type --> type xml element
        for el in root.getchildren():
            is_std_type = False
            try:
                is_std_type = utils.isStdType(el)
            except Exception:
                pass
            if is_std_type:
                if 'demangled' in el.keys():
                    std_type_name = el.get('demangled')    
                else:
                    std_type_name = el.get('name')

                cfg.std_types_dict[std_type_name] = el
                
                # print 'GOT STD TYPE: ', el.get('name')

        # for el in root.findall('Class'):
        #     if utils.isStdType(el):
        #         cfg.std_types_dict[el.get('name')] = el
        # for el in root.findall('Struct'):
        #     if utils.isStdType(el):
        #         cfg.std_types_dict[el.get('name')] = el


        # Update global dict: typedef name --> typedef xml element
        cfg.typedef_dict = OrderedDict() 
        for el in root.findall('Typedef'):

            # Only accept native typedefs:
            if utils.isNative(el):

                typedef_name = el.get('name')

                try:
                    type_name, type_kw, type_id = utils.findType(el)
                except:
                    warnings.warn("Could not determine type of xml element with 'id'=%s" % el.get('id')) 
                    continue
                type_el = cfg.id_dict[type_id]


                # If underlying type is a fundamental or standard type, accept it right away
                if utils.isFundamental(type_el) or utils.isStdType(type_el):
                    cfg.typedef_dict[typedef_name] = el
                
                # If underlying type is a class/struct, check if it's acceptable
                elif type_el.tag in ['Class', 'Struct']:
                    if 'demangled' in type_el.keys():
                        complete_type_name = type_el.get('demangled').replace(' ','')
                    else:
                        complete_type_name = type_el.get('name')
                    
                    if complete_type_name in cfg.accepted_classes:
                        cfg.typedef_dict[typedef_name] = el
                
                # If neither fundamental or class/struct, ignore it.
                else:
                    pass

        # print
        # print '****'
        # for td in cfg.typedef_dict.keys():
        #     print td
        # print '****'
        # print

        # Update global dict: class name --> class xml element
        for el in root.findall('Class'):
            demangled_class_name = el.get('demangled').replace(' ','')
            if demangled_class_name in cfg.accepted_classes:
                cfg.class_dict[demangled_class_name] = el
            elif (options.list_classes_flag) and (utils.isNative(el)):
                cfg.class_dict[demangled_class_name] = el


        # Update global dict: function name --> function xml element
        for el in root.findall('Function'):
            if 'demangled' in el.keys():
                demangled_name = el.get('demangled')
                # Template functions have more complicated 'demangled' entries...
                if '<' in demangled_name:
                    func_name_full = demangled_name.split(' ',1)[1].split('(',1)[0]
                else:
                    func_name_full = demangled_name.split('(',1)[0]
                func_name_full = func_name_full.replace(' ','')
                if func_name_full in cfg.accepted_functions:
                    cfg.func_dict[func_name_full] = el
                elif (options.list_classes_flag) and (utils.isNative(el)):
                    cfg.func_dict[func_name_full] = el


        # Update global list: accepted types
        fundamental_types  = [ el.get('name') for el in root.findall('FundamentalType')]
        cfg.accepted_types = fundamental_types + cfg.std_types_dict.keys() + cfg.accepted_classes + cfg.typedef_dict.keys()


        #
        # If -l option is given, print a list of all avilable classes and functions, then exit.
        #

        if options.list_classes_flag:

            print 'Classes:'
            print '--------'
            for demangled_class_name, class_el in cfg.class_dict.items():
                print ' - ' + demangled_class_name
            print

            print 'Functions:'
            print '----------'
            for func_name_full, func_el in cfg.func_dict.items():
                extended_name = func_el.get('demangled').replace(' ','')
                print ' - ' + extended_name
            print

            sys.exit()


        #
        # Parse classes
        #

        temp_new_code = classparse.run()

        for src_file_name in temp_new_code.keys():
            if src_file_name not in new_code.keys():
                new_code[src_file_name] = []
            for code_tuple in temp_new_code[src_file_name]:
                new_code[src_file_name].append(code_tuple)

        #
        # Parse functions
        #

        temp_new_code = funcparse.run()

        for src_file_name in temp_new_code.keys():
            if src_file_name not in new_code.keys():
                new_code[src_file_name] = []
            for code_tuple in temp_new_code[src_file_name]:
                new_code[src_file_name].append(code_tuple)

        #
        # Create header with typedefs
        #
        code_tuple = utils.constrTypedefHeader()


        # typedef_code = ''
        # insert_pos   = 0

        # # - Forward declarations
        # typedef_code += '// Forward declarations:\n'
        # for class_name_full, class_el in cfg.class_dict.items():

        #     class_name_short       = class_name_full.split('<',1)[0]
        #     abstr_class_name       = classutils.getAbstractClassName(class_name_full)
        #     abstr_class_name_short = abstr_class_name.split('<',1)[0]

        #     # print 'class_name_full        ', class_name_full
        #     # print 'class_name_short       ', class_name_short
        #     # print 'abstr_class_name       ', abstr_class_name
        #     # print 'abstr_class_name_short ', abstr_class_name_short


        #     if '<' in class_name_full:
        #         is_template = True
        #     else:
        #         is_template = False

        #     if is_template:
        #         template_bracket = utils.getTemplateBracket(class_el)[0]
        #         # spec_template_types = utils.getSpecTemplateTypes(class_el)
        #         # spec_template_bracket = '<' + ','.join(spec_template_types) + '>'

        #     if is_template:
        #         typedef_code += 'template ' + template_bracket + '\n'
        #         typedef_code += 'class ' + abstr_class_name_short + ';\n'
                
        #         typedef_code += 'template ' + template_bracket + '\n'
        #         typedef_code += 'class ' + class_name_short + ';\n'
                
        #         typedef_code += 'class ' + abstr_class_name + ';\n'
        #         typedef_code += 'class ' + class_name_full + ';\n'
        #     else:
        #         typedef_code += 'class ' + abstr_class_name + ';\n'
        #         typedef_code += 'class ' + class_name_full + ';\n'
        # typedef_code += '\n'

        # # - Typedefs
        # typedef_code += '// Typedefs:\n'
        # for typedef_name, typedef_el in cfg.typedef_dict.items():
        #     type_name, type_kw, type_id = utils.findType(typedef_el)
        #     typedef_code += 'typedef ' + type_name + ' ' + typedef_name + ';\n'

        # code_tuple = (insert_pos, typedef_code)
        
        typedef_file_path = os.path.join(cfg.extra_output_dir, cfg.all_typedefs_fname)
        if typedef_file_path not in new_code.keys():
            new_code[typedef_file_path] = []
        new_code[typedef_file_path].append(code_tuple)



    #
    # END: loop over xml files
    #

    #
    # (Write new code to source files)
    #

    for src_file_name in new_code.keys():

        code_tuples = new_code[src_file_name]
        code_tuples.sort( key=lambda x : x[0], reverse=True )

        new_src_file_name  = os.path.join(cfg.extra_output_dir, os.path.basename(src_file_name))

        if code_tuples == []:
            continue

        print 
        print
        print 'FILE : ', new_src_file_name
        print '========================================='


        if os.path.isfile(src_file_name):
            if options.debug_mode_flag == True:
                f = open(src_file_name, 'a+')
            else:
                f = open(src_file_name, 'r')
            f.seek(0)
            file_content = f.read()
            f.close()
            new_file_content = file_content

        else:
            new_file_content = ''

        for pos,code in code_tuples:
            if pos == -1:
                new_file_content = new_file_content + code    
            else:
                new_file_content = new_file_content[:pos] + code + new_file_content[pos:]

        # Do the writing!
        if options.debug_mode_flag == True:
            pass
        else:
            f = open(new_src_file_name, 'w')
            f.write(new_file_content)
            f.close()

        print new_file_content
        print


    print 'STANDARD TYPES USED:'
    print '===================='
    for std_type in cfg.std_types_used:
        print std_type
    print

    print 'ALL TYPES USED:'
    print '==============='
    for type_name in cfg.all_types_in_class:
        print type_name
    print


# ====== END: main ========

if  __name__ =='__main__':main()

