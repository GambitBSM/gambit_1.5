#!/usr/bin/python

########################################
#                                      #
#   BOSS - Backend-On-a-Stick Script   #
#                                      #
########################################

#
# This is a Python package for making the classes of a C++ library 
# available for dynamic use through the 'dlopen' system.
# BOSS relies on 'gccxml' being installed and callable from the command line.
#
# Default usage:
# ./boss [list of class header files]
# 

import xml.etree.ElementTree as ET
import os
import sys
import warnings
import shutil
import glob
import subprocess
from collections import OrderedDict
from optparse import OptionParser

import modules.cfg as cfg
import modules.gb as gb
import modules.classutils as classutils
import modules.classparse as classparse
import modules.funcparse as funcparse
import modules.funcutils as funcutils
import modules.utils as utils
import modules.filehandling as filehandling


# ====== main ========

def main():

    print
    print
    print '    ========================================'
    print '    ||                                    ||'
    print '    ||  BOSS - Backend-On-a-Stick-Script  ||'
    print '    ||                                    ||'
    print '    ========================================'
    print 
    print 


    # Parse command line arguments and options
    parser = OptionParser(usage="usage: %prog [options] <header files>",
                          version="%prog 0.1")
    parser.add_option("-l", "--list",
                      action="store_true",
                      dest="list_flag",
                      default=False,
                      help="Output a list of the available classes and functions.")
    parser.add_option("-b", "--backup-sources",
                      action="store_true",
                      dest="backup_sources_flag",
                      default=False,
                      help="Create backup source files (blah.cpp.boss) before BOSS mangles them.")
    parser.add_option("-d", "--debug-mode",
                      action="store_true",
                      dest="debug_mode_flag",
                      default=False,
                      help="Test run. No source files are changed.")
    (options, args) = parser.parse_args()


    # Check that arguments list is not empty
    if len(args) == 0:

        print 
        print 'Missing input arguments. For instructions, run boss.py --help'
        print 

        # Exit
        sys.exit()



    # Get the input file names from command line input. 
    input_files = args

    # Sort them to make sure screen output is identical regardless of ordering of input files.
    input_files.sort()


    #
    # Run gccxml for all input header/source files
    #

    print
    print 'Parsing the input files:'
    print '------------------------'
    print 

    xml_files = []
    for input_file_path in input_files:

        # Get path and filename for the input file
        input_file_dir, input_file_short_name = os.path.split(input_file_path)

        # Construct file name for xml file produced by gccxml
        xml_output_path = os.path.join('temp', input_file_path.replace('/','_').replace('.','_') + '.xml' )

        # List all include paths
        include_paths_list = [cfg.include_path] + cfg.additional_include_paths

        # Timeout limit and process poll interval [seconds]
        timeout = 20.
        poll = 0.2

        # Run gccxml
        try:
            utils.gccxmlRunner(input_file_path, include_paths_list, xml_output_path, timeout_limit=timeout, poll_interval=poll)
        except:
            raise

        # Append xml file to list of xml files
        xml_files.append(xml_output_path)

    #
    # END: Run gccxml on input files
    #
    print


    #
    # If -l option is given, print a list of all classes and functions, then exit.
    #

    if options.list_flag:

        all_class_names    = []
        all_function_names = []

        for xml_file in xml_files:

            tree = ET.parse(xml_file)
            root = tree.getroot()

            # Set the global xml id dict. (Needed by the functions called from utils.)
            gb.id_dict = OrderedDict([ (el.get('id'), el) for el in root.getchildren() ]) 
           
            # Find all available classes
            for el in (root.findall('Class') + root.findall('Struct')):
                
                # Skip classes that are not loadable (incomplete, abstract, ...)
                try:
                    is_loadable = utils.isLoadable(el)
                except KeyError:
                    continue

                if not is_loadable:
                    continue

                demangled_class_name = el.get('demangled')
                if utils.isNative(el):
                    all_class_names.append(demangled_class_name)

            # Find all available functions
            for el in root.findall('Function'):
                if 'demangled' in el.keys():
                    demangled_name = el.get('demangled')

                    # Template functions have more complicated 'demangled' entries...
                    if '<' in demangled_name:
                        func_name_full = demangled_name.split(' ',1)[1].split('(',1)[0]
                    else:
                        func_name_full = demangled_name.split('(',1)[0]

                    if utils.isNative(el):
                        all_function_names.append(func_name_full)

        # END: Loop over xml files

        # Remove duplicates
        all_class_names    = list(OrderedDict.fromkeys(all_class_names))
        all_function_names = list(OrderedDict.fromkeys(all_function_names))
        
        # Output lists
        print 'Classes:'
        print '--------'
        for demangled_class_name in all_class_names:
            print ' - ' + demangled_class_name
        print
        print 'Functions:'
        print '----------'
        for demangled_func_name in all_function_names:
            print ' - ' + demangled_func_name
        print

        # Exit
        sys.exit()


    #
    # Initialization
    #

    # Create the output directories if they do not exist.
    filehandling.createOutputDirectories()


    #
    # Remove from cfg.loaded_classes all classes that are not loadabe (not found, incomplete, abstract, ...)
    #

    # Remove duplicates from cfg.loaded_classes
    cfg.loaded_classes = list(OrderedDict.fromkeys(cfg.loaded_classes))

    # Create list of names for all loadable classes
    loadable_classes_names_list = []

    for xml_file in xml_files:

        tree = ET.parse(xml_file)
        root = tree.getroot()

        # Set the global xml id dict. (Needed by the functions called from utils.)
        gb.id_dict = OrderedDict([ (el.get('id'), el) for el in root.getchildren() ]) 

        # Loop all class xml elements found in the current xml file
        for class_el in (root.findall('Class') + root.findall('Struct')):

            # Store name of any loadable classes
            try:
                is_loadable = utils.isLoadable(class_el)
            except KeyError:
                pass

            if is_loadable:
                if 'demangled' in class_el.keys():
                    loadable_classes_names_list.append(class_el.get('demangled'))

    # Remove unloadable classes from cfg.loaded_classes
    for class_name_long_templ in cfg.loaded_classes[:]:
        if class_name_long_templ not in loadable_classes_names_list:
            cfg.loaded_classes.remove(class_name_long_templ)
            print 'INFO: Class %s is not loadable. Possible explanations: not found, incomplete, abstract, ...' % (class_name_long_templ)



    #
    # TODO: Remove from cfg.loaded_functions all functions that are not loadable
    #



    #
    # Main loop over all xml files
    #

    for current_xml_file in xml_files:

        # Reset some global variables for each new xml file
        gb.xml_file_name = ''
        gb.id_dict.clear()
        gb.file_dict.clear()
        gb.all_classes_dict.clear()
        gb.std_types_dict.clear()
        gb.typedef_dict.clear()
        gb.class_dict.clear()
        gb.loaded_classes_in_xml.clear()
        gb.func_dict.clear()

        # Output xml file name
        print
        print ':::::  Current XML file: %s  :::::' % current_xml_file
        print

        # Set global xml file name and parse it using ElementTree
        gb.xml_file_name = current_xml_file
        tree = ET.parse(gb.xml_file_name)
        root = tree.getroot()

        # Update global dict: id --> xml element (all elements)
        gb.id_dict = OrderedDict([ (el.get('id'), el) for el in root.getchildren() ]) 


        # Update global dict: file name --> file xml element
        gb.file_dict = OrderedDict([ (el.get('name'), el) for el in root.findall('File') ])


        # Update global dict: class name --> class xml element (all classes)
        for class_el in (   root.findall('Class') 
                          + root.findall('Struct') 
                          + root.findall('FundamentalType') 
                          + root.findall('Typedef') ):

            if 'demangled' in class_el.keys():
                class_name_long_templ = class_el.get('demangled')
                gb.all_classes_dict[class_name_long_templ] = class_el
            elif 'name' in class_el.keys():
                class_name_long_templ = class_el.get('name')
                gb.all_classes_dict[class_name_long_templ] = class_el


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

                gb.std_types_dict[std_type_name] = el
                

        # Update global dict: class name --> class xml element
        for el in (root.findall('Class') + root.findall('Struct')):

            try:
                class_name = classutils.getClassNameDict(el)
            except KeyError:
                continue

            # Check if we have done this class already
            if class_name in gb.classes_done:
                print 'INFO: (Class %s already done)' % (class_name['long_templ'])
                continue

            # Check that class is requested
            if (class_name['long_templ'] in cfg.loaded_classes):

                # Skip classes that are not loadable (incomplete, abstract, ...) 
                # (This should not be needed as we already check all classes listed in cfg.loaded_classes...)
                if not utils.isLoadable(el, print_warning=True):
                    continue

                # Store class xml element
                gb.loaded_classes_in_xml[class_name['long_templ']] = el



        # Update global dict: typedef name --> typedef xml element
        gb.typedef_dict = OrderedDict() 
        for el in root.findall('Typedef'):

            # Only accept native typedefs:
            if utils.isNative(el):

                typedef_name = el.get('name')

                type_dict = utils.findType(el)
                type_el = type_dict['el']

                # If underlying type is a fundamental or standard type, accept it right away
                if utils.isFundamental(type_el) or utils.isStdType(type_el):
                    gb.typedef_dict[typedef_name] = el
                
                # If underlying type is a class/struct, check if it's acceptable
                elif type_el.tag in ['Class', 'Struct']:
                    if 'demangled' in type_el.keys():
                        complete_type_name = type_el.get('demangled')
                    else:
                        complete_type_name = type_el.get('name')
                    
                    if complete_type_name in cfg.loaded_classes:
                        gb.typedef_dict[typedef_name] = el
                
                # If neither fundamental or class/struct, ignore it.
                else:
                    pass


        # Update global dict: function name --> function xml element
        for el in root.findall('Function'):
            if 'demangled' in el.keys():
                demangled_name = el.get('demangled')

                # Template functions have more complicated 'demangled' entries...
                if '<' in demangled_name:
                    func_name_full = demangled_name.split(' ',1)[1].split('(',1)[0]
                else:
                    func_name_full = demangled_name.split('(',1)[0]

                if func_name_full in cfg.loaded_functions:
                    gb.func_dict[func_name_full] = el


        # Update global dict: function name --> function xml element
        for el in root.findall('Function'):
            if 'demangled' in el.keys():
                demangled_func_name = el.get('demangled')
                # Template functions have more complicated 'demangled' entries...
                if '<' in demangled_func_name:
                    func_name_full = demangled_func_name.split(' ',1)[1].split('(',1)[0]
                else:
                    func_name_full = demangled_func_name.split('(',1)[0]
                if func_name_full in cfg.loaded_functions:
                    gb.func_dict[func_name_full] = el



        # If requested, append any (native) parent classes to the cfg.loaded_classes list
        if cfg.load_parent_classes:

            for class_name, class_el in gb.loaded_classes_in_xml.items():

                parents_el_list = utils.getAllParentClasses(class_el, only_native_classes=True)

                for el in parents_el_list:

                    # Skip classes that are not loadable (incomplete, abstract, ...)
                    if not utils.isLoadable(el, print_warning=True):
                        continue

                    if 'demangled' in el.keys():
                        demangled_name = el.get('demangled')
                        
                        # - Update cfg.loaded_classes
                        if demangled_name not in cfg.loaded_classes:
                            cfg.loaded_classes.append(demangled_name)
                        
                        # - Update gb.class_dict
                        if demangled_name not in gb.loaded_classes_in_xml.keys():
                            gb.loaded_classes_in_xml[demangled_name] = el


        # Update global list: accepted types
        fundamental_types  = [ el.get('name') for el in root.findall('FundamentalType')]
        enumeration_types  = [ '::'.join( utils.getNamespaces(el, include_self=True) ) for el in root.findall('Enumeration')]

        gb.accepted_types  = fundamental_types + enumeration_types + gb.std_types_dict.keys() + cfg.loaded_classes + gb.typedef_dict.keys()

        # Remove from gb.accepted_types all classes that use of native types as template arguments
        # (BOSS cannot deal with this yet...)
        for i in range(len(gb.accepted_types))[::-1]:

            type_name = gb.accepted_types[i]

            # Get list of all template arguments (unpack any nested template arguments)
            unpacked_template_args = []
            utils.unpackAllSpecTemplateTypes(type_name, unpacked_template_args)

            # If no template arguments, continue
            if unpacked_template_args == []:
                continue

            else:
                # print 'TEMPL_ARGS:', unpacked_template_args
                for templ_arg in unpacked_template_args:

                    # Remove asterix and/or ampersand
                    base_templ_arg = utils.getBasicTypeName(templ_arg)

                    # Check that this type is listed in gb.all_classes_dict (all native types should be)
                    if base_templ_arg in gb.all_classes_dict.keys():

                        # Get xml entry for the type
                        class_el = gb.all_classes_dict[base_templ_arg]

                        # If this is a native type, remove the current entry (i) in gb.accepted_types
                        if utils.isNative(class_el):
                            
                            # Remove entry i from gb.accepted_types
                            gb.accepted_types.pop(i)
                            break


        # print
        # print 'ALL CLASSES:'
        # print '------------'
        # for entry in gb.all_classes_dict.keys():
        #     print 'ALL CLASSES:', entry
        # print

        # print
        # print 'STD CLASSES:'
        # print '------------'
        # for entry in gb.std_types_dict.keys():
        #     print 'STD CLASSES:', entry
        # print

        # print
        # print 'ACCEPTED TYPES:'
        # print '---------------'
        # for entry in gb.accepted_types:
        #     print 'ACCEPTED:', entry
        # print


        # Update global dict: new header files
        for class_name in cfg.loaded_classes:
            
            class_name_short = class_name.split('<',1)[0].split('::')[-1]
            class_name_long  = class_name.split('<',1)[0]

            if class_name_long not in gb.new_header_files.keys():
              
                abstract_header_name     = cfg.abstr_header_prefix + class_name_short + cfg.header_extension
                wrapper_header_name      = cfg.wrapper_header_prefix + class_name_short + cfg.header_extension
                wrapper_decl_header_name = cfg.wrapper_header_prefix + class_name_short + '_decl' + cfg.header_extension
                wrapper_def_header_name  = cfg.wrapper_header_prefix + class_name_short + '_def'  + cfg.header_extension

                abstract_header_fullpath     = os.path.join(gb.gambit_backend_types_basedir, gb.gambit_backend_name_full, cfg.abstr_header_prefix + class_name_short + cfg.header_extension )
                wrapper_header_fullpath      = os.path.join(gb.gambit_backend_types_basedir, gb.gambit_backend_name_full, cfg.wrapper_header_prefix + class_name_short + cfg.header_extension )
                wrapper_decl_header_fullpath = os.path.join(gb.gambit_backend_types_basedir, gb.gambit_backend_name_full, cfg.wrapper_header_prefix + class_name_short + '_decl' + cfg.header_extension )
                wrapper_def_header_fullpath  = os.path.join(gb.gambit_backend_types_basedir, gb.gambit_backend_name_full, cfg.wrapper_header_prefix + class_name_short + '_def'  + cfg.header_extension )
                
                gb.new_header_files[class_name_long] = {    'abstract': abstract_header_name, 
                                                            'wrapper': wrapper_header_name, 
                                                            'wrapper_decl': wrapper_decl_header_name,
                                                            'wrapper_def': wrapper_def_header_name,
                                                            'abstract_fullpath': abstract_header_fullpath, 
                                                            'wrapper_fullpath': wrapper_header_fullpath, 
                                                            'wrapper_decl_fullpath': wrapper_decl_header_fullpath,
                                                            'wrapper_def_fullpath': wrapper_def_header_fullpath    }


        #
        # Parse classes
        #

        classparse.run()


        #
        # Parse functions
        #

        funcparse.run()


        #
        # Create header with forward declarations of all abstract classes
        #

        abs_frwd_decls_header_path = os.path.join(cfg.extra_output_dir, gb.frwd_decls_abs_fname + cfg.header_extension)
        utils.constrAbsForwardDeclHeader(abs_frwd_decls_header_path)


        #
        # Create header with forward declarations of all wrapper classes
        #

        wrp_frwd_decls_header_path = os.path.join(cfg.extra_output_dir, gb.frwd_decls_wrp_fname + cfg.header_extension)
        utils.constrWrpForwardDeclHeader(wrp_frwd_decls_header_path)
        

        #
        # Create header with declarations of all enum types
        #

        enum_decls_header_path = os.path.join(cfg.extra_output_dir, gb.enum_decls_wrp_fname + cfg.header_extension)
        utils.constrEnumDeclHeader(root.findall('Enumeration'), enum_decls_header_path)


    #
    # END: loop over xml files
    #


    #
    # Write new files
    #

    for src_file_name, code_dict in gb.new_code.iteritems():

        add_include_guard = code_dict['add_include_guard']
        code_tuples = code_dict['code_tuples']

        code_tuples.sort( key=lambda x : x[0], reverse=True )

        new_src_file_name  = os.path.join(cfg.extra_output_dir, os.path.basename(src_file_name))

        if code_tuples == []:
            continue

        boss_backup_exists = False
        if os.path.isfile(src_file_name):
            try:
                f = open(src_file_name + '.boss', 'r')
                boss_backup_exists = True
            except IOError, e:
                if e.errno != 2:
                    raise e
                f = open(src_file_name, 'r')

            f.seek(0)
            file_content = f.read()
            f.close()
            new_file_content = file_content
        else:
            new_file_content = ''

        if options.backup_sources_flag and not boss_backup_exists and new_file_content:
            f = open(src_file_name + '.boss', 'w')
            f.write(new_file_content)
            f.close()

        for pos,code in code_tuples:

            if pos == -1:
                new_file_content = new_file_content + code    
            else:
                new_file_content = new_file_content[:pos] + code + new_file_content[pos:]

        # Add include guard where requested
        if add_include_guard:

            short_new_src_file_name = os.path.basename(new_src_file_name)

            if 'include_guard_prefix' in code_dict.keys():
                prefix = code_dict['include_guard_prefix']
            else:
                prefix = ''

            new_file_content = utils.addIncludeGuard(new_file_content, short_new_src_file_name, prefix=prefix ,suffix=gb.gambit_backend_name_full)

        # Do the writing! (If debug_mode, only print file content to screen)
        if options.debug_mode_flag:
            print 
            print
            print 'FILE : ', new_src_file_name
            print '========================================='
            print new_file_content
            print
        else:
            f = open(new_src_file_name, 'w')
            f.write(new_file_content)
            f.close()



    # 
    # Copy files from common_headers/ and replace any code template tags
    # 

    filehandling.createCommonHeaders()


    #
    # Move files to correct directories
    #

    filehandling.moveFilesAround()


    #
    # Run through all the generated files and use the code tags __START_GAMBIT_NAMESPACE__ and __END_GAMBIT_NAMESPACE__ to construct
    # the correct namespace.
    #

    construct_namespace_in_files = glob.glob( os.path.join(gb.gambit_backend_dir_complete, '*') )

    filehandling.replaceNamespaceTags(construct_namespace_in_files, gb.gambit_backend_namespace, '__START_GAMBIT_NAMESPACE__', '__END_GAMBIT_NAMESPACE__')


    #
    # Run through all the generated files and remove tags that are no longer needed
    #

    all_generated_files = glob.glob( os.path.join(cfg.extra_output_dir, '*') ) + glob.glob( os.path.join(gb.gambit_backend_dir_complete, '*') )
    remove_tags_list = [ '__START_GAMBIT_NAMESPACE__', 
                         '__END_GAMBIT_NAMESPACE__', 
                         '__INSERT_CODE_HERE__' ]

    filehandling.removeCodeTagsFromFiles(all_generated_files, remove_tags_list)


    #
    # Copy files to the correct locations within the source tree of the original code
    #

    print
    print
    print 'Copying generated files to original source tree:'
    print '------------------------------------------------'
    print 

    filehandling.copyFilesToSourceTree(verbose=True)


    #
    # Parse all factory function source files using gccxml
    #

    print
    print 
    print 
    print 'Parsing the generated factory function source files:'
    print '----------------------------------------------------'
    print 

    factory_xml_files = OrderedDict()
    for class_name in gb.classes_done:

        # Construct factory file name
        factory_source_fname_short = cfg.factory_file_prefix + class_name['short'] + cfg.source_extension
        factory_source_path        = os.path.join(cfg.source_path, factory_source_fname_short)

        # Construct file name for xml file produced by gccxml
        xml_output_path = os.path.join('temp', factory_source_path.replace('/','_').replace('.','_') + '.xml' )

        # List all include paths
        include_paths_list = [cfg.include_path] + cfg.additional_include_paths

        # Timeout limit and process poll interval [seconds]
        timeout = 20.
        poll = 0.2

        # Run gccxml
        try:
            utils.gccxmlRunner(factory_source_path, include_paths_list, xml_output_path, timeout_limit=timeout, poll_interval=poll)
        except:
            raise

        # Add factory xml file to dict
        factory_xml_files[class_name['long']] = xml_output_path

    #
    # END: Parse all factory function source files
    #


    #
    # Generate header file 'loaded_types.hpp'
    #

    print
    print
    print 'Generating file loaded_types.hpp:'
    print '---------------------------------'
    print 

    filehandling.createLoadedTypesHeader(factory_xml_files)



    #
    # Done!
    #

    print
    print 'Done!'
    print '-----' 
    print

# ====== END: main ========

if  __name__ =='__main__':main()

