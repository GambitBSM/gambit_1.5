#!/usr/bin/python
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
import shutil
import glob
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

    print
    print

    #
    # Initialization
    #

    # Prepare container variables
    new_code = OrderedDict()

    # Create the extra output directory if it does not exist
    try:
      os.mkdir(cfg.extra_output_dir)
    except OSError, e:
      if e.errno != 17:
        raise e

    # Set GAMBIT variables
    cfg.gambit_backend_safeversion = cfg.gambit_backend_version.replace('.','_')

    cfg.all_wrapper_fname = 'loadedtypes_' + cfg.gambit_backend_name.lower() + '_' + cfg.gambit_backend_safeversion.lower()

    # cfg.gambit_full_namespace = cfg.gambit_base_namespace.rstrip('::') + '::' + cfg.gambit_backend_name + '_' + cfg.gambit_backend_safeversion

    cfg.gambit_backend_name_full = cfg.gambit_backend_name + '_' + cfg.gambit_backend_safeversion


    # Initialize variables
    cfg.n_classes_done = 0


    # Parse command line arguments and options
    parser = OptionParser(usage="usage: %prog [options] xmlfiles",
                          version="%prog 0.1")
    parser.add_option("-l", "--list",
                      action="store_true",
                      dest="list_flag",
                      default=False,
                      help="Output a list of the available classes and functions.")
    parser.add_option("-c", "--choose",
                      action="store_true",
                      dest="choose_flag",
                      default=False,
                      help="Choose from a list of the available classes and functions at runtime.")
    parser.add_option("-p", "--set-path-id",
                      action="store_true",
                      dest="set_path_ids_flag",
                      default=False,
                      help="Set the path IDs used to consider a class or function part of a package.")
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

    if len(args) < 1:
        parser.error("No xml file(s) specified")

    # Get the xml file names from command line input
    xml_files = args

    # Set up a few more things before the file loop
    if options.choose_flag:
        print 'Overwriting any values in cfg.loaded_classes and cfg.loaded_functions...'
        cfg.loaded_classes   = []
        cfg.loaded_functions = []
    if options.set_path_ids_flag:
        print 'Overwriting any values in cfg.accepted_paths...'
        cfg.accepted_paths = []
        print 'Path ID example: "pythia" can ID any file within \n'
        print '    "/anywhere/pythiaxxxx/..." as being part of pythia.\n'
        print 'Enter all the path IDs to use:' 
        while(True):
            # Also allows for comma separated lists
            paths = raw_input('(blank to stop)\n').partition(',')
            if all([path.strip() == '' for path in paths]):
                break
            while(paths[0].strip()):
                cfg.accepted_paths.append(paths[0].strip())
                paths = paths[2].partition(',')


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
                
        # for std_type in cfg.std_types_dict.keys():
        #     print [std_type]


        # for el in root.findall(('Class') + root.findall('Struct')):
        #     if utils.isStdType(el):
        #         cfg.std_types_dict[el.get('name')] = el
        # for el in root.findall('Struct'):
        #     if utils.isStdType(el):
        #         cfg.std_types_dict[el.get('name')] = el


        # Update global dict: class name --> class xml element
        # Classes before typedefs, so the user can choose the classes they want
        if options.choose_flag:
            # ultimately:   terminal_size = (rows - 8, columns - 8)
            terminal_size = os.popen('stty size', 'r').read().split()
            terminal_size = (int(terminal_size[0]) - 8, int(terminal_size[1]) - 8)
            qtr = terminal_size[1] / 4
            coords = [0, 0]
            count = 0
            mini_class_dict = {}
            for el in (root.findall('Class') + root.findall('Struct')):
                demangled_class_name = el.get('demangled')
                if coords[0] == 0:
                    # then, new screen of classes to choose from
                    coords[0] += 1
                    count = 0
                    mini_class_dict.clear()
                    print '\n\n Choose what classes you need from this list: \n\n   ',
                if utils.isNative(el):
                    count += 1
                    mini_class_dict[str(count)] = demangled_class_name
                    item = str(count) + '. ' + demangled_class_name
                    item += ' '*(qtr - (len(item)+1)%qtr)
                    if coords[1] and len(item) + coords[1] + 1 >= terminal_size[1]:
                        coords[0] += 1
                        coords[1] = 0
                        print '\n   ',
                    coords[1] += len(item) + 1
                    print item,
                    if coords[0] >= terminal_size[0]:
                        coords[0] = 0
                        choices = raw_input('\n\n Input a comma separated list of the numbers indicating the classes you want:\n').partition(',')
                        while(choices[0].strip()):
                            cfg.loaded_classes.append(mini_class_dict[choices[0].strip()])
                            choices = choices[2].partition(',')
            else:
                choices = raw_input('\n\n Input a comma separated list of the numbers indicating the classes you want:\n').partition(',')
                while(choices[0].strip()):
                    cfg.loaded_classes.append(mini_class_dict[choices[0].strip()])
                    choices = choices[2].partition(',')
        for el in (root.findall('Class') + root.findall('Struct')):
            demangled_class_name = el.get('demangled')
            if demangled_class_name in cfg.loaded_classes:
                cfg.class_dict[demangled_class_name] = el
            elif options.list_flag and utils.isNative(el):
                cfg.class_dict[demangled_class_name] = el


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
                        complete_type_name = type_el.get('demangled')
                    else:
                        complete_type_name = type_el.get('name')
                    
                    if complete_type_name in cfg.loaded_classes:
                        cfg.typedef_dict[typedef_name] = el
                
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
                    cfg.func_dict[func_name_full] = el
                elif (options.list_flag) and (utils.isNative(el)):
                    cfg.func_dict[func_name_full] = el


        # Update global dict: function name --> function xml element
        if options.choose_flag:
            # # ultimately:   terminal_size = (rows - 8, columns - 8)
            # terminal_size = os.popen('stty size', 'r').read().split()
            # terminal_size = (int(terminal_size[0]) - 8, int(terminal_size[1]) - 8)
            # qtr = terminal_size[1] / 4
            coords = [0, 0]
            count = 0
            mini_func_dict = {}
            for el in root.findall('Function'):
                if 'demangled' in el.keys():
                    demangled_func_name = el.get('demangled')
                    # Template functions have more complicated 'demangled' entries...
                    if '<' in demangled_func_name:
                        func_name_full = demangled_func_name.split(' ',1)[1].split('(',1)[0]
                    else:
                        func_name_full = demangled_func_name.split('(',1)[0]

                    if coords[0] == 0:
                        # then, new screen of functions to choose from
                        coords[0] += 1
                        count = 0
                        mini_func_dict.clear()
                        print '\n\n Choose what functions you need from this list: \n\n   ',
                    if utils.isNative(el):
                        count += 1
                        mini_func_dict[str(count)] = func_name_full
                        item = str(count) + '. ' + func_name_full
                        item += ' '*(qtr - (len(item)+1)%qtr)
                        if coords[1] and len(item) + coords[1] + 1 >= terminal_size[1]:
                            coords[0] += 1
                            coords[1] = 0
                            print '\n   ',
                        coords[1] += len(item) + 1
                        print item,
                        if coords[0] >= terminal_size[0]:
                            coords[0] = 0
                            choices = raw_input('\n\n Input a comma separated list of the numbers indicating the functions you want:\n').partition(',')
                            while(choices[0].strip()):
                                cfg.loaded_functions.append(mini_func_dict[choices[0].strip()])
                                choices = choices[2].partition(',')
            else:
                choices = raw_input('\n\n Input a comma separated list of the numbers indicating the functions you want:\n').partition(',')
                while(choices[0].strip()):
                    cfg.loaded_functions.append(mini_func_dict[choices[0].strip()])
                    choices = choices[2].partition(',')
        for el in root.findall('Function'):
            if 'demangled' in el.keys():
                demangled_func_name = el.get('demangled')
                # Template functions have more complicated 'demangled' entries...
                if '<' in demangled_func_name:
                    func_name_full = demangled_func_name.split(' ',1)[1].split('(',1)[0]
                else:
                    func_name_full = demangled_func_name.split('(',1)[0]
                if func_name_full in cfg.loaded_functions:
                    cfg.func_dict[func_name_full] = el
                elif (options.list_flag) and (utils.isNative(el)):
                    cfg.func_dict[func_name_full] = el



        # If requested, append any (native) parent classes to the cfg.loaded_classes list
        if cfg.load_parent_classes:

            for class_name in cfg.loaded_classes:
                class_el = cfg.class_dict[class_name]
                parents_el_list = utils.getAllParentClasses(class_el, only_native_classes=True)

                for el in parents_el_list:
                    if 'demangled' in el.keys():
                        demangled_name = el.get('demangled')
                        
                        # - Update cfg.loaded_classes
                        if demangled_name not in cfg.loaded_classes:
                            cfg.loaded_classes.append(demangled_name)
                        
                        # - Update cfg.class_dict
                        if demangled_name not in cfg.class_dict.keys():
                            cfg.class_dict[demangled_name] = el


        # Update global list: accepted types
        fundamental_types  = [ el.get('name') for el in root.findall('FundamentalType')]
        cfg.accepted_types = fundamental_types + cfg.std_types_dict.keys() + cfg.loaded_classes + cfg.typedef_dict.keys()


        # Update global dict: new header files
        for class_name in cfg.loaded_classes:
            
            class_name_short = class_name.split('<',1)[0].split('::')[-1]
            class_name_long  = class_name.split('<',1)[0]

            if class_name_long not in cfg.new_header_files.keys():
              
                abstract_header_name     = cfg.abstr_header_prefix + class_name_short + cfg.header_extension
                wrapper_header_name      = cfg.wrapper_header_prefix + class_name_short + cfg.header_extension
                wrapper_decl_header_name = cfg.wrapper_header_prefix + class_name_short + '_decl' + cfg.header_extension
                wrapper_def_header_name  = cfg.wrapper_header_prefix + class_name_short + '_def'  + cfg.header_extension

                abstract_header_fullpath     = os.path.join(cfg.gambit_backend_basedir, cfg.gambit_backend_name_full, cfg.abstr_header_prefix + class_name_short + cfg.header_extension )
                wrapper_header_fullpath      = os.path.join(cfg.gambit_backend_basedir, cfg.gambit_backend_name_full, cfg.wrapper_header_prefix + class_name_short + cfg.header_extension )
                wrapper_decl_header_fullpath = os.path.join(cfg.gambit_backend_basedir, cfg.gambit_backend_name_full, cfg.wrapper_header_prefix + class_name_short + '_decl' + cfg.header_extension )
                wrapper_def_header_fullpath  = os.path.join(cfg.gambit_backend_basedir, cfg.gambit_backend_name_full, cfg.wrapper_header_prefix + class_name_short + '_def'  + cfg.header_extension )
                
                cfg.new_header_files[class_name_long] = {   'abstract': abstract_header_name, 
                                                            'wrapper': wrapper_header_name, 
                                                            'wrapper_decl': wrapper_decl_header_name,
                                                            'wrapper_def': wrapper_def_header_name,
                                                            'abstract_fullpath': abstract_header_fullpath, 
                                                            'wrapper_fullpath': wrapper_header_fullpath, 
                                                            'wrapper_decl_fullpath': wrapper_decl_header_fullpath,
                                                            'wrapper_def_fullpath': wrapper_def_header_fullpath    }



        #
        # If -l option is given, print a list of all avilable classes and functions, then exit.
        #

        if options.list_flag:

            if not options.choose_flag:
                print 'Classes:'
                print '--------'
                for demangled_class_name in cfg.class_dict:
                    print ' - ' + demangled_class_name
                print

                print 'Functions:'
                print '----------'
                for func_el in cfg.func_dict.values():
                    extended_name = func_el.get('demangled')
                    print ' - ' + extended_name
                print

            sys.exit()


        # #
        # # Write header with base wrapper class
        # #

        # wrapper_base_class_name = 'WrapperBase'

        # code_tuple = classutils.generateWrapperBaseHeader(wrapper_base_class_name)

        # wrapper_base_header_fname = cfg.wrapper_header_prefix + wrapper_base_class_name + cfg.header_extension
        # wrapper_base_header_path  = os.path.join(cfg.extra_output_dir, wrapper_base_header_fname)

        # # - Update cfg.new_header_files
        # if wrapper_base_class_name not in cfg.new_header_files:
        #     cfg.new_header_files[wrapper_base_class_name] = {'abstract': '', 'wrapper': wrapper_base_header_fname}

        # # - Register header file in the 'new_code' dict
        # if wrapper_base_header_path not in new_code.keys():
        #     new_code[wrapper_base_header_path] = {'code_tuples':[], 'add_include_guard':True}
        # new_code[wrapper_base_header_path]['code_tuples'].append(code_tuple)


        #
        # Parse classes
        #

        temp_new_code = classparse.run()

        for src_file_name in temp_new_code.keys():

            if src_file_name not in new_code:
                new_code[src_file_name] = {'code_tuples':[], 'add_include_guard':False}

            new_code[src_file_name]['code_tuples']      += temp_new_code[src_file_name]['code_tuples']
            new_code[src_file_name]['add_include_guard'] = temp_new_code[src_file_name]['add_include_guard']


        # for src_file_name in temp_new_code:
        #     if src_file_name not in new_code:
        #         new_code[src_file_name] = []
        #     for code_tuple in temp_new_code[src_file_name]:
        #         new_code[src_file_name].append(code_tuple)
        #     new_code[s]


        #
        # Parse functions
        #

        temp_new_code = funcparse.run()

        for src_file_name in temp_new_code.keys():

            if src_file_name not in new_code:
                new_code[src_file_name] = {'code_tuples':[], 'add_include_guard':False}

            new_code[src_file_name]['code_tuples']      += temp_new_code[src_file_name]['code_tuples']
            new_code[src_file_name]['add_include_guard'] = temp_new_code[src_file_name]['add_include_guard']

        # for src_file_name in temp_new_code:
        #     if src_file_name not in new_code:
        #         new_code[src_file_name] = []
        #     for code_tuple in temp_new_code[src_file_name]:
        #         new_code[src_file_name].append(code_tuple)


        #
        # Create header with typedefs
        #

        code_tuple = utils.constrTypedefHeader()

        typedef_header_path = os.path.join(cfg.extra_output_dir, cfg.all_typedefs_fname + cfg.header_extension)
        
        if typedef_header_path not in new_code.keys():
            new_code[typedef_header_path] = {'code_tuples':[], 'add_include_guard':True}
        new_code[typedef_header_path]['code_tuples'].append(code_tuple)



        #
        # Create header with forward declarations of all abstract classes
        #

        code_tuple = utils.constrAbsForwardDeclHeader()

        abs_frwd_decls_header_path = os.path.join(cfg.extra_output_dir, cfg.frwd_decls_abs_fname + cfg.header_extension)
        
        if abs_frwd_decls_header_path not in new_code.keys():
            new_code[abs_frwd_decls_header_path] = {'code_tuples':[], 'add_include_guard':True}
        new_code[abs_frwd_decls_header_path]['code_tuples'].append(code_tuple)


        #
        # Create header with forward declarations of all wrapper classes
        #

        code_tuple = utils.constrWrpForwardDeclHeader()

        wrp_frwd_decls_header_path = os.path.join(cfg.extra_output_dir, cfg.frwd_decls_wrp_fname + cfg.header_extension)
        
        if wrp_frwd_decls_header_path not in new_code.keys():
            new_code[wrp_frwd_decls_header_path] = {'code_tuples':[], 'add_include_guard':True}
        new_code[wrp_frwd_decls_header_path]['code_tuples'].append(code_tuple)


    #
    # END: loop over xml files
    #


    # #
    # # If cfg.gambit_mode = True, add namespace tags to thw wrapper typedefs header
    # #

    # if cfg.gambit_mode:
    #     wrapper_typedefs_path = cfg.wrapper_typedefs_fname + cfg.header_extension
    #     if wrapper_typedefs_path not in new_code.keys():
    #         new_code[wrapper_typedefs_path] = {'code_tuples':[], 'add_include_guard':True}
    #     new_code[wrapper_typedefs_path]['code_tuples'].append( (0, '__START_GAMBIT_NAMESPACE__\n') )
    #     new_code[wrapper_typedefs_path]['code_tuples'].append( (-1,'__END_GAMBIT_NAMESPACE__\n') )


    #
    # Add include statement at beginning of wrapper deleter source file
    #    
    w_deleter_include = '#include "' + os.path.join(cfg.wrapper_deleter_fname + cfg.header_extension) + '"\n'
    w_deleter_source_path = os.path.join(cfg.extra_output_dir, cfg.wrapper_deleter_fname + cfg.source_extension)
    if w_deleter_source_path not in new_code.keys():
        new_code[w_deleter_source_path] = {'code_tuples':[], 'add_include_guard':False}
    new_code[w_deleter_source_path]['code_tuples'].append( (0, w_deleter_include) )


    # #
    # # Add include statement at beginning of wrapper deleter header file
    # #    
    # w_deleter_include = '#include "' + os.path.join(cfg.add_path_to_includes, cfg.all_wrapper_fname + '_decl' + cfg.header_extension) + '"\n'
    # w_deleter_header_path = os.path.join(cfg.extra_output_dir, cfg.wrapper_deleter_fname + cfg.header_extension)
    # if w_deleter_header_path not in new_code.keys():
    #     new_code[w_deleter_header_path] = {'code_tuples':[], 'add_include_guard':True}
    # new_code[w_deleter_header_path]['code_tuples'].append( (0, w_deleter_include) )



    #
    # Write new code to source files
    #

    for src_file_name, code_dict in new_code.iteritems():

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
                if not code in new_file_content:
                    new_file_content = new_file_content + code    
                else:
                    # Ignore duplicate code. (Generate warning is duplicated code is not blank.)
                    if code.strip() != '':
                        warnings.warn('Ignored duplicate code: \n' + code + '\n\n')
                    else:
                        pass

            elif new_file_content[pos:].find(code) != 0:
                new_file_content = new_file_content[:pos] + code + new_file_content[pos:]
            else:
                # Ignore duplicate code. (Generate warning is duplicated code is not blank.)
                if code.strip() != '':
                    warnings.warn('Ignored duplicate code: \n' + code + '\n\n')
                else:
                    pass

        # Add include guard where requested
        if add_include_guard:
            short_new_src_file_name = os.path.basename(new_src_file_name)
            new_file_content = utils.addIncludeGuard(new_file_content, short_new_src_file_name, extra_string=cfg.gambit_backend_name_full)

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

    # #
    # # Generate 'abstracts' header
    # #

    # abstracts_includes_header_path = os.path.join(cfg.extra_output_dir, 'abstracts_includes_TEMP.hpp')
    # abstracts_typedefs_header_path = os.path.join(cfg.extra_output_dir, 'abstracts_typedefs.hpp')
    # abstracts_header_path          = os.path.join(cfg.extra_output_dir, 'abstracts' + cfg.header_extension)

    # f = open(abstracts_includes_header_path, 'r')
    # abstracts_includes_header_content = f.read()
    # f.close()

    # f = open(abstracts_typedefs_header_path, 'r')
    # abstracts_typedefs_header_content = f.read()
    # f.close()

    # abstracts_header_content  = '\n' 
    # abstracts_header_content += '// Abstract classes:\n'
    # abstracts_header_content += '// -----------------\n'
    # abstracts_header_content += '\n' 
    # abstracts_header_content += abstracts_includes_header_content
    # abstracts_header_content += '\n' 
    # abstracts_header_content += '// Typedefs for abstract classes:\n'
    # abstracts_header_content += '// ------------------------------\n'
    # abstracts_header_content += '\n' 
    # abstracts_header_content += '#include "' + os.path.join(cfg.gambit_backend_basedir, cfg.gambit_backend_name_full, 'identification.hpp') + '"\n'
    # abstracts_header_content += abstracts_typedefs_header_content
    # abstracts_header_content += '#include "backend_undefs.hpp"\n'
    # abstracts_header_content += '\n' 


    # # Add include guards
    # abstracts_header_content = utils.addIncludeGuard(abstracts_header_content, os.path.basename(abstracts_header_path))

    # # Write the 'abstracts' header (if not in debug_mode)
    # if not options.debug_mode_flag:
    #     f = open(abstracts_header_path, 'w')
    #     f.write(abstracts_header_content)
    #     f.close()

    # # Remove the temporary headers 'abstracts_includes_TEMP.hpp' and 'abstracts_typedefs_TEMP.hpp' we used to construct the master header
    # os.remove(abstracts_includes_header_path)
    # # os.remove(abstracts_typedefs_header_path)



    # #
    # # Generate 'wrappers' header
    # #

    # wrappers_decl_includes_header_path = os.path.join(cfg.extra_output_dir, 'wrappers_decl_includes_TEMP.hpp')
    # wrappers_def_includes_header_path  = os.path.join(cfg.extra_output_dir, 'wrappers_def_includes_TEMP.hpp')
    # wrappers_typedefs_header_path      = os.path.join(cfg.extra_output_dir, 'wrappers_typedefs_TEMP.hpp')
    # wrappers_header_path               = os.path.join(cfg.extra_output_dir, 'wrappers' + cfg.header_extension)

    # f = open(wrappers_decl_includes_header_path, 'r')
    # wrappers_decl_includes_header_content = f.read()
    # f.close()

    # f = open(wrappers_def_includes_header_path, 'r')
    # wrappers_def_includes_header_content = f.read()
    # f.close()

    # f = open(wrappers_typedefs_header_path, 'r')
    # wrappers_typedefs_header_content = f.read()
    # f.close()

    # wrappers_header_content  = '\n' 
    # wrappers_header_content += '// Wrapper classes:\n'
    # wrappers_header_content += '// ----------------\n'
    # wrappers_header_content += '\n' 
    # wrappers_header_content += wrappers_decl_includes_header_content
    # wrappers_header_content += '\n' 
    # wrappers_header_content += wrappers_def_includes_header_content
    # wrappers_header_content += '\n' 
    # wrappers_header_content += '// Typedefs for wrapper classes:\n'
    # wrappers_header_content += '// -----------------------------\n'
    # wrappers_header_content += '\n' 
    # wrappers_header_content += '#include "' + os.path.join(cfg.gambit_backend_basedir, cfg.gambit_backend_name_full, 'identification.hpp') + '"\n'
    # wrappers_header_content += wrappers_typedefs_header_content
    # wrappers_header_content += '#include "backend_undefs.hpp"\n'
    # wrappers_header_content += '\n' 


    # # Add include guards
    # wrappers_header_content = utils.addIncludeGuard(wrappers_header_content, os.path.basename(wrappers_header_path))

    # # Write the 'wrappers' header (if not in debug_mode)
    # if not options.debug_mode_flag:
    #     f = open(wrappers_header_path, 'w')
    #     f.write(wrappers_header_content)
    #     f.close()

    # # Remove the temporary headers 'wrappers_decl_includes_TEMP.hpp', 'wrappers_def_includes_TEMP.hpp' and 'wrappers_typedefs_TEMP.hpp' 
    # # we used to construct the master header
    # os.remove(wrappers_decl_includes_header_path)
    # os.remove(wrappers_def_includes_header_path)
    # os.remove(wrappers_typedefs_header_path)




    # #
    # # Generate the master header
    # #

    # all_wrapper_decl_header_path = os.path.join(cfg.extra_output_dir, cfg.all_wrapper_fname + '_decl' + cfg.header_extension)
    # all_wrapper_def_header_path  = os.path.join(cfg.extra_output_dir, cfg.all_wrapper_fname + '_def' + cfg.header_extension)
    # master_header_path           = os.path.join(cfg.extra_output_dir, cfg.all_wrapper_fname + cfg.header_extension)

    # f = open(all_wrapper_decl_header_path, 'r')
    # all_wrapper_decl_header_content = f.read()
    # f.close()

    # f = open(all_wrapper_def_header_path, 'r')
    # all_wrapper_def_header_content = f.read()
    # f.close()

    # master_header_content  = '\n' 
    # master_header_content += '// Loaded types declarations:\n'
    # master_header_content += '// --------------------------\n'
    # master_header_content += '\n' 
    # master_header_content += all_wrapper_decl_header_content
    # master_header_content += '\n' 
    # master_header_content += '// Loaded types implementations:\n'
    # master_header_content += '// -----------------------------\n'
    # master_header_content += '\n' 
    # master_header_content += all_wrapper_def_header_content
    # master_header_content += '\n' 
    # master_header_content += '// Typedefs: wrapper class name --> original class name:\n'
    # master_header_content += '// -----------------------------------------------------\n'
    # master_header_content += '\n' 
    # master_header_content += '#include "' + cfg.wrapper_typedefs_fname + cfg.header_extension + '"\n' 


    # # Add include guards
    # master_header_content = utils.addIncludeGuard(master_header_content, cfg.all_wrapper_fname)

    # # Write the master header (if not in debug_mode)
    # if not options.debug_mode_flag:
    #     f = open(master_header_path, 'w')
    #     f.write(master_header_content)
    #     f.close()

    # # # Remove the temporary '_decl' and '_def' header files we used to construct the master header
    # # os.remove(all_wrapper_decl_header_path)
    # # os.remove(all_wrapper_def_header_path)




    # 
    # Copy and rename files from ./common_headers to output directory
    # 

    # - abstractbase.hpp

    source_file_name = 'common_headers/abstractbase.hpp'
    target_file_name = os.path.join(cfg.extra_output_dir, 'abstractbase.hpp')
    shutil.copyfile(source_file_name, target_file_name)


    # - wrapperbase.hpp

    source_file_name = 'common_headers/wrapperbase.hpp'
    target_file_name = os.path.join(cfg.extra_output_dir, 'wrapperbase.hpp')
    
    f_src = open(source_file_name, 'r')
    content = f_src.read()
    f_src.close()

    new_content = content.replace('__CODE_SUFFIX__', cfg.code_suffix)

    f_target = open(target_file_name, 'w')
    f_target.write(new_content)
    f_target.close()


    # # - boss_loaded_classes.hpp

    # source_file_name = 'common_headers/boss_loaded_classes.hpp'
    # target_file_name = os.path.join(cfg.extra_output_dir, 'boss_loaded_classes' + cfg.header_extension)

    # f_src = open(source_file_name, 'r')
    # content = f_src.read()
    # f_src.close()

    # new_content = content.replace('__HEADER_EXTENSION__', cfg.header_extension)

    # f_target = open(target_file_name, 'w')
    # f_target.write(new_content)
    # f_target.close()


    # - nullptr_check.hpp

    # source_file_name = 'common_headers/nullptr_check.hpp'
    # target_file_name = os.path.join(cfg.extra_output_dir, 'nullptr_check' + cfg.header_extension)
    # shutil.copy (source_file_name, target_file_name)


    # 
    # Organize output files: A base folder with files needed only by the external code, and a subfolder
    # with the files that should be copied to GAMBIT
    #

    gambit_backend_dir_complete = os.path.join(cfg.extra_output_dir, cfg.gambit_backend_basedir, cfg.gambit_backend_name_full) 

    # - Create directories
    try:
        os.mkdir(cfg.extra_output_dir)
    except OSError, e:
        if e.errno == 17:
            pass
        else:
            raise

    if cfg.gambit_backend_basedir != '':
        try:
            os.mkdir( os.path.join(cfg.extra_output_dir, cfg.gambit_backend_basedir) )
        except OSError, e:
            if e.errno == 17:
                pass
            else:
                raise

    try:
        os.mkdir( gambit_backend_dir_complete )
    except OSError, e:
        if e.errno == 17:
            pass
        else:
            raise


    # NEW: move files to correct directories

    # - To gambit_backend_dir_complete
    move_files_list  = []

    # -- abstract class headers    
    move_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.abstr_header_prefix + '*') )

    # -- wrapper class headers
    move_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.wrapper_header_prefix + '*') )

    # -- header with forward declarations for all abstract classes
    move_files_list += [ os.path.join(cfg.extra_output_dir, cfg.frwd_decls_abs_fname + cfg.header_extension) ]

    # -- header with forward declarations for all wrapper classes
    move_files_list += [ os.path.join(cfg.extra_output_dir, cfg.frwd_decls_wrp_fname + cfg.header_extension) ]

    for mv_file in move_files_list:
        shutil.move(mv_file, gambit_backend_dir_complete)


    
    # NEW: copy files to correct directories

    # - To gambit_backend_dir_complete
    copy_files_list  = []
    
    copy_files_list += ['headers_by_hand/loaded_types.hpp' ]
    copy_files_list += ['headers_by_hand/identification.hpp' ]

    for cp_file in copy_files_list:
        shutil.copy(cp_file, gambit_backend_dir_complete)


    # - To cfg.extra_output_dir
    copy_files_list  = []

    copy_files_list += ['common_headers/cats.hpp' ]
    copy_files_list += ['common_headers/backend_undefs.hpp' ]

    for cp_file in copy_files_list:
        shutil.copy(cp_file, cfg.extra_output_dir)



    # # - First, move files that should *only* go into gambit_backend_dir_complete
    # move_files_list  = []
    # move_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.wrapper_typedefs_fname + cfg.header_extension) )   # + wrapper typedefs header
    # move_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.all_wrapper_fname + cfg.header_extension) )        # + master header that includes all wrapper headers + wrapper typedefs

    # for mv_file in move_files_list:
    #     shutil.move(mv_file, gambit_backend_dir_complete)

    # # - Next, copy to gambit_backend_dir_complete the files that should go in both gambit_backend_dir_complete and cfg.extra_output_dir
    # copy_files_list  = []
    # copy_files_list += glob.glob( os.path.join(cfg.extra_output_dir, 'AbstractBase' + cfg.header_extension) )             # + abstract base class header
    # copy_files_list += glob.glob( os.path.join(cfg.extra_output_dir, 'WrapperBase' + cfg.header_extension) )              # + wrapper base class header
    # copy_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.abstr_header_prefix + '*') )                     # + abstract class headers
    # copy_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.wrapper_header_prefix + '*') )                   # + wrapper class headers
    # copy_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.all_wrapper_fname + cfg.header_extension) )      # + header including all wrapper headers
    # copy_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.frwd_decls_abs_fname + cfg.header_extension) )   # + header with forward declarations for all abstract classes
    # copy_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.frwd_decls_wrp_fname + cfg.header_extension) )   # + header with forward declarations for all wrapper classes
    
    # for cp_file in copy_files_list:
    #     shutil.copy(cp_file, gambit_backend_dir_complete)

    # # - Finally move all (remaining) files to cfg.extra_output_dir
    # move_files_list  = []
    # move_files_list += glob.glob( os.path.join(cfg.extra_output_dir, '*' + cfg.header_extension) )   # + header files
    # move_files_list += glob.glob( os.path.join(cfg.extra_output_dir, '*' + cfg.source_extension) )   # + source files

    # for mv_file in move_files_list:
    #     shutil.move(mv_file, cfg.extra_output_dir)


    # # #
    # # # If cfg.gambit_mode = True, in the directory gambit_backend_dir_complete the 'master header' cfg.all_wrapper_fname should
    # # # be renamed to loadedtypes_BACKENDNAME_VERSION.hpp
    # # #
    # # if cfg.gambit_mode:
    # #     move_from = os.path.join(gambit_backend_dir_complete, cfg.all_wrapper_fname + cfg.header_extension)
    # #     move_to   = os.path.join(gambit_backend_dir_complete, 'loadedtypes_' + cfg.gambit_backend_name.lower() + '_' + cfg.gambit_backend_safeversion.lower() + cfg.header_extension)
    # #     shutil.move(move_from, move_to)


    #
    # Run through all the generated files and use the code tags __START_GAMBIT_NAMESPACE__ and __END_GAMBIT_NAMESPACE__ to construct
    # the correct namespace. If cfg.gambit_mode = False, then simply remove these tags
    #

    open_namespace_tag  = '__START_GAMBIT_NAMESPACE__'
    close_namespace_tag = '__END_GAMBIT_NAMESPACE__'

    remove_tags_from_files = glob.glob( os.path.join(cfg.extra_output_dir, '*') )
    construct_namespace_in_files = []
    if cfg.gambit_mode:
        construct_namespace_in_files += glob.glob( os.path.join(gambit_backend_dir_complete, '*') )
    else:
        remove_tags_from_files += glob.glob( os.path.join(gambit_backend_dir_complete, '*') )


    # - Remove namespace tags from given files
    for file_path in remove_tags_from_files:

        if not os.path.isfile(file_path):
            continue

        f = open(file_path, 'r')
        content = f.read()
        f.close()

        new_content = content.replace(open_namespace_tag, '').replace(close_namespace_tag, '')

        f = open(file_path, 'w')
        f.write(new_content)
        f.close()


    # - Construct namespace in given files
    for file_path in construct_namespace_in_files:

        if not os.path.isfile(file_path):
            continue

        f = open(file_path, 'r')
        content = f.read()
        f.close()

        new_content = utils.constrNamespaceFromTags(content, cfg.gambit_backend_namespace, open_namespace_tag, close_namespace_tag)

        f = open(file_path, 'w')
        f.write(new_content)
        f.close()


    print 
    print



# ====== END: main ========

if  __name__ =='__main__':main()

