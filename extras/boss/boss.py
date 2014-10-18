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


# ====== main ========

def main():

    print
    print '    ================================='
    print '    ||                             ||'
    print '    ||  BOSS - Backend-On-a-Stick  ||'
    print '    ||                             ||'
    print '    ================================='
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


    # Get the input file names from command line input
    input_files = args

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
    # Run gccxml for all input header/source files
    #
    xml_files = []
    for input_file_path in input_files:

        # Get path and filename for the input file
        input_file_dir, input_file_short_name = os.path.split(input_file_path)

        # Construct file name for xml file produced by gccxml
        # xml_file_path = os.path.join('temp', os.path.splitext(input_file_short_name)[0] + '.xml' )
        xml_file_path = os.path.join('temp', input_file_path.replace('/','_').replace('.','_') + '.xml' )


        # Construct gccxml command to run
        gccxml_cmd = 'gccxml '

        # - Add include paths
        if cfg.include_path != '':
            gccxml_cmd += '-I' + cfg.include_path + ' '
        for add_incl_path in cfg.additional_include_paths:
            gccxml_cmd += '-I' + add_incl_path + ' '

        # - Add the input file (full path)
        gccxml_cmd += input_file_path

        # - Add gccxml option that specifies the xml output file: input_file_short_name.xml
        gccxml_cmd += ' -fxml=' + xml_file_path


        # Run gccxml
        print 'Runing command: ' + gccxml_cmd
        proc = subprocess.Popen(gccxml_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait()

        if proc.returncode != 0:        
            print 'ERROR: gccxml failed. Printing output:'
            output = proc.stderr.read()
            print 
            print '==== START GCCXML OUTPUT ===='
            print
            print output
            print '==== END GCCXML OUTPUT ===='
            print
            sys.exit()
        
        else:
            print 'Command finished successfully.'


        # Append xml file to list of xml files
        xml_files.append(xml_file_path)

    #
    # END: Run gccxml on input files
    #



    #
    # Loop over all the xml files
    #

    for current_xml_file in xml_files:

        # Output xml file name
        print
        print '=====:  Current XML file: %s  :=====' % current_xml_file
        print

        # Set global xml file name and parse it using ElementTree
        gb.xml_file_name = current_xml_file
        tree = ET.parse(gb.xml_file_name)
        root = tree.getroot()


        # Update global dict: id --> xml element (all elements)
        gb.id_dict = OrderedDict([ (el.get('id'), el) for el in root.getchildren() ]) 


        # Update global dict: file name --> file xml element
        gb.file_dict = OrderedDict([ (el.get('name'), el) for el in root.findall('File') ])


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
                
        # for std_type in gb.std_types_dict.keys():
        #     print [std_type]


        # for el in root.findall(('Class') + root.findall('Struct')):
        #     if utils.isStdType(el):
        #         gb.std_types_dict[el.get('name')] = el
        # for el in root.findall('Struct'):
        #     if utils.isStdType(el):
        #         gb.std_types_dict[el.get('name')] = el


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
                gb.class_dict[demangled_class_name] = el
            elif options.list_flag and utils.isNative(el):
                gb.class_dict[demangled_class_name] = el


        # Update global dict: typedef name --> typedef xml element
        gb.typedef_dict = OrderedDict() 
        for el in root.findall('Typedef'):

            # Only accept native typedefs:
            if utils.isNative(el):

                typedef_name = el.get('name')

                try:
                    type_name, type_kw, type_id = utils.findType(el)
                except:
                    warnings.warn("Could not determine type of xml element with 'id'=%s" % el.get('id')) 
                    continue
                type_el = gb.id_dict[type_id]


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
                elif (options.list_flag) and (utils.isNative(el)):
                    gb.func_dict[func_name_full] = el


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
                    gb.func_dict[func_name_full] = el
                elif (options.list_flag) and (utils.isNative(el)):
                    gb.func_dict[func_name_full] = el



        # If requested, append any (native) parent classes to the cfg.loaded_classes list
        if cfg.load_parent_classes:

            for class_name in cfg.loaded_classes:
                class_el = gb.class_dict[class_name]
                parents_el_list = utils.getAllParentClasses(class_el, only_native_classes=True)

                for el in parents_el_list:
                    if 'demangled' in el.keys():
                        demangled_name = el.get('demangled')
                        
                        # - Update cfg.loaded_classes
                        if demangled_name not in cfg.loaded_classes:
                            cfg.loaded_classes.append(demangled_name)
                        
                        # - Update gb.class_dict
                        if demangled_name not in gb.class_dict.keys():
                            gb.class_dict[demangled_name] = el


        # Update global list: accepted types
        fundamental_types  = [ el.get('name') for el in root.findall('FundamentalType')]
        gb.accepted_types  = fundamental_types + gb.std_types_dict.keys() + cfg.loaded_classes + gb.typedef_dict.keys()


        # Update global dict: new header files
        for class_name in cfg.loaded_classes:
            
            class_name_short = class_name.split('<',1)[0].split('::')[-1]
            class_name_long  = class_name.split('<',1)[0]

            if class_name_long not in gb.new_header_files.keys():
              
                abstract_header_name     = cfg.abstr_header_prefix + class_name_short + cfg.header_extension
                wrapper_header_name      = cfg.wrapper_header_prefix + class_name_short + cfg.header_extension
                wrapper_decl_header_name = cfg.wrapper_header_prefix + class_name_short + '_decl' + cfg.header_extension
                wrapper_def_header_name  = cfg.wrapper_header_prefix + class_name_short + '_def'  + cfg.header_extension

                abstract_header_fullpath     = os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, cfg.abstr_header_prefix + class_name_short + cfg.header_extension )
                wrapper_header_fullpath      = os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, cfg.wrapper_header_prefix + class_name_short + cfg.header_extension )
                wrapper_decl_header_fullpath = os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, cfg.wrapper_header_prefix + class_name_short + '_decl' + cfg.header_extension )
                wrapper_def_header_fullpath  = os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, cfg.wrapper_header_prefix + class_name_short + '_def'  + cfg.header_extension )
                
                gb.new_header_files[class_name_long] = {   'abstract': abstract_header_name, 
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
                for demangled_class_name in gb.class_dict:
                    print ' - ' + demangled_class_name
                print

                print 'Functions:'
                print '----------'
                for func_el in gb.func_dict.values():
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

        # # - Update gb.new_header_files
        # if wrapper_base_class_name not in gb.new_header_files:
        #     gb.new_header_files[wrapper_base_class_name] = {'abstract': '', 'wrapper': wrapper_base_header_fname}

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


        # #
        # # Create header with typedefs
        # #

        # code_tuple = utils.constrTypedefHeader()

        # typedef_header_path = os.path.join(cfg.extra_output_dir, gb.all_typedefs_fname + cfg.header_extension)
        
        # if typedef_header_path not in new_code.keys():
        #     new_code[typedef_header_path] = {'code_tuples':[], 'add_include_guard':True}
        # new_code[typedef_header_path]['code_tuples'].append(code_tuple)



        #
        # Create header with forward declarations of all abstract classes
        #

        code_tuple = utils.constrAbsForwardDeclHeader()

        abs_frwd_decls_header_path = os.path.join(cfg.extra_output_dir, gb.frwd_decls_abs_fname + cfg.header_extension)
        
        if abs_frwd_decls_header_path not in new_code.keys():
            new_code[abs_frwd_decls_header_path] = {'code_tuples':[], 'add_include_guard':True}
        new_code[abs_frwd_decls_header_path]['code_tuples'].append(code_tuple)


        #
        # Create header with forward declarations of all wrapper classes
        #

        code_tuple = utils.constrWrpForwardDeclHeader()

        wrp_frwd_decls_header_path = os.path.join(cfg.extra_output_dir, gb.frwd_decls_wrp_fname + cfg.header_extension)
        
        if wrp_frwd_decls_header_path not in new_code.keys():
            new_code[wrp_frwd_decls_header_path] = {'code_tuples':[], 'add_include_guard':True}
        new_code[wrp_frwd_decls_header_path]['code_tuples'].append(code_tuple)


    #
    # END: loop over xml files
    #

    #
    # Add include statement at beginning of wrapper deleter source file
    #    
    w_deleter_include = '#include "' + os.path.join(gb.wrapper_deleter_fname + cfg.header_extension) + '"\n'
    w_deleter_source_path = os.path.join(cfg.extra_output_dir, gb.wrapper_deleter_fname + cfg.source_extension)
    if w_deleter_source_path not in new_code.keys():
        new_code[w_deleter_source_path] = {'code_tuples':[], 'add_include_guard':False}
    new_code[w_deleter_source_path]['code_tuples'].append( (0, w_deleter_include) )


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
            new_file_content = utils.addIncludeGuard(new_file_content, short_new_src_file_name, extra_string=gb.gambit_backend_name_full)

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
    # Copy and edit files from common_headers/
    # 

    # - abstractbase.hpp

    source_file_name = 'common_headers/abstractbase.hpp'
    target_file_name = os.path.join(cfg.extra_output_dir, 'abstractbase.hpp')
    shutil.copyfile(source_file_name, target_file_name)


    # - wrapperbase.hpp

    source_file_name = 'common_headers/wrapperbase.hpp'
    target_file_name = os.path.join(cfg.extra_output_dir, 'wrapperbase.hpp')

    new_content = utils.replaceCodeTags(source_file_name, file_input=True)

    f_target = open(target_file_name, 'w')
    f_target.write(new_content)
    f_target.close()


    # - identification.hpp

    source_file_name = 'common_headers/identification.hpp'
    target_file_name = os.path.join(cfg.extra_output_dir, 'identification.hpp')

    new_content = utils.replaceCodeTags(source_file_name, file_input=True)

    f_target = open(target_file_name, 'w')
    f_target.write(new_content)
    f_target.close()



    # 
    # Organize output files: A base folder with files needed only by the external code, and a subfolder
    # with the files that should be copied to GAMBIT
    #

    gambit_backend_dir_complete = os.path.join(cfg.extra_output_dir, cfg.gambit_backend_basedir, gb.gambit_backend_name_full) 

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


    # Move files to correct directories

    # - To gambit_backend_dir_complete
    move_files_list  = []

    # -- abstract class headers    
    move_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.abstr_header_prefix + '*') )

    # -- wrapper class headers
    move_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.wrapper_header_prefix + '*') )

    # -- header with forward declarations for all abstract classes
    move_files_list += [ os.path.join(cfg.extra_output_dir, gb.frwd_decls_abs_fname + cfg.header_extension) ]

    # -- header with forward declarations for all wrapper classes
    move_files_list += [ os.path.join(cfg.extra_output_dir, gb.frwd_decls_wrp_fname + cfg.header_extension) ]

    # -- identification.hpp
    move_files_list += [ os.path.join(cfg.extra_output_dir, 'identification.hpp') ]    

    for mv_file in move_files_list:
        shutil.move(mv_file, gambit_backend_dir_complete)


    
    # Copy files to correct directories

    # - To gambit_backend_dir_complete
    copy_files_list  = []
    
    copy_files_list += ['headers_by_hand/loaded_types.hpp' ]

    for cp_file in copy_files_list:
        shutil.copy(cp_file, gambit_backend_dir_complete)


    # - To cfg.extra_output_dir
    copy_files_list  = []

    copy_files_list += ['common_headers/cats.hpp' ]
    copy_files_list += ['common_headers/backend_undefs.hpp' ]

    for cp_file in copy_files_list:
        shutil.copy(cp_file, cfg.extra_output_dir)


    #
    # Run through all the generated files and use the code tags __START_GAMBIT_NAMESPACE__ and __END_GAMBIT_NAMESPACE__ to construct
    # the correct namespace.
    #

    open_namespace_tag  = '__START_GAMBIT_NAMESPACE__'
    close_namespace_tag = '__END_GAMBIT_NAMESPACE__'

    remove_tags_from_files       = glob.glob( os.path.join(cfg.extra_output_dir, '*') )
    construct_namespace_in_files = glob.glob( os.path.join(gambit_backend_dir_complete, '*') )


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

        new_content = utils.constrNamespaceFromTags(content, gb.gambit_backend_namespace, open_namespace_tag, close_namespace_tag)

        f = open(file_path, 'w')
        f.write(new_content)
        f.close()



    #
    # Copy files to the correct locations within the source tree of the original code
    #

    print
    print 'Copying generated files to original source tree:'
    print '------------------------------------------------'
    print 


    # - Copy all manipulated original files
    for original_file_short_name, original_file_full_path in gb.original_file_paths.items():

        cp_source = os.path.join(cfg.extra_output_dir,original_file_short_name)
        cp_target = original_file_full_path

        if os.path.isfile(cp_source):
            shutil.copyfile(cp_source, cp_target)
            print '   ' + cp_source + '   --->   ' + cp_target
    print


    # - Copy factory source files
    for class_name_long in gb.classes_done:

        class_name_short = utils.removeNamespace(class_name_long, return_namespace=False)
        class_name_short = utils.removeTemplateBracket(class_name_short)

        factory_source_fname_short = cfg.factory_file_prefix + class_name_short + cfg.source_extension

        cp_source = os.path.join(cfg.extra_output_dir, factory_source_fname_short)
        cp_target = os.path.join(cfg.source_path, factory_source_fname_short)

        if os.path.isfile(cp_source):
            shutil.copyfile(cp_source, cp_target)
            print '   ' + cp_source + '   --->   ' + cp_target
    print


    # - Copy 'extras' source files (containing implementations for the helper functions that BOSS adds to the original classes)
    for class_name_long in gb.classes_done:

        class_name_short = utils.removeNamespace(class_name_long, return_namespace=False)
        class_name_short = utils.removeTemplateBracket(class_name_short)

        extra_source_fname_short = class_name_short + '_extras' + gb.code_suffix + cfg.source_extension

        cp_source = os.path.join(cfg.extra_output_dir, extra_source_fname_short)
        cp_target = os.path.join(cfg.source_path, extra_source_fname_short)

        if os.path.isfile(cp_source):
            shutil.copyfile(cp_source, cp_target)
            print '   ' + cp_source + '   --->   ' + cp_target
    print


    # - Copy files that are common for all BOSS runs
    files_to_incl_path = []
    files_to_incl_path.append( 'abstractbase.hpp' )                                      # abstractbase.hpp
    files_to_incl_path.append( gb.abstract_typedefs_fname + cfg.header_extension )       # abstracttypedefs.hpp
    files_to_incl_path.append( 'backend_undefs.hpp' )                                    # backend_undefs.hpp
    files_to_incl_path.append( 'cats.hpp' )                                              # cats.hpp
    files_to_incl_path.append( 'wrapperbase.hpp' )                                       # wrapperbase.hpp
    files_to_incl_path.append( gb.wrapper_deleter_fname + cfg.header_extension )         # wrapperdeleter.hpp
    files_to_incl_path.append( gb.wrapper_typedefs_fname + cfg.header_extension )        # wrappertypedefs.hpp

    for file_name_short in files_to_incl_path:
        cp_source = os.path.join(cfg.extra_output_dir, file_name_short)
        cp_target = os.path.join(cfg.include_path, file_name_short)

        if os.path.isfile(cp_source):
            shutil.copyfile(cp_source, cp_target)
            print '   ' + cp_source + '   --->   ' + cp_target
    print

    files_to_src_path = []
    files_to_src_path.append( gb.wrapper_deleter_fname + cfg.source_extension )         # wrapperdeleter.cpp

    for file_name_short in files_to_src_path:
        cp_source = os.path.join(cfg.extra_output_dir, file_name_short)
        cp_target = os.path.join(cfg.source_path, file_name_short)

        if os.path.isfile(cp_source):
            shutil.copyfile(cp_source, cp_target)
            print '   ' + cp_source + '   --->   ' + cp_target
    print


    # - Copy the entire backend_types/ directory to the include directory
    cp_source = os.path.join(cfg.extra_output_dir, cfg.gambit_backend_basedir)
    cp_target = os.path.join(cfg.include_path, cfg.gambit_backend_basedir)

    if os.path.isdir(cp_source):
        print '   ' + cp_source + '   --->   ' + cp_target
        shutil.copytree(cp_source, cp_target)




    print 
    print

# ====== END: main ========

if  __name__ =='__main__':main()

