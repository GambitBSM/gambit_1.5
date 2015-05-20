######################################
#                                    #
#  File handling functions for BOSS  #
#                                    #
######################################

import xml.etree.ElementTree as ET
from collections import OrderedDict
import os
import warnings
import subprocess
import shutil
import glob

import modules.cfg as cfg
import modules.gb as gb
import modules.utils as utils



# ====== createOutputDirectories ========

# Creates the directory structure for the files generate by BOSS.
# Base directory is cfg.extra_output_dir.

def createOutputDirectories():

    try:
        os.mkdir(cfg.extra_output_dir)
    except OSError, e:
        if e.errno == 17:
            pass
        else:
            raise

    if gb.gambit_backend_types_basedir != '':
        try:
            os.mkdir( os.path.join(cfg.extra_output_dir, gb.gambit_backend_types_basedir) )
        except OSError, e:
            if e.errno == 17:
                pass
            else:
                raise

    try:
        os.mkdir( gb.gambit_backend_dir_complete )
    except OSError, e:
        if e.errno == 17:
            pass
        else:
            raise

# ====== END: createOutputDirectories ========



# ====== moveFilesAround ========

# Move files to their correct subdirectories within cfg.extra_output_dir.

def moveFilesAround():

    # - To gb.gambit_backend_dir_complete
    move_files_list  = []

    # -- abstract class headers    
    move_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.abstr_header_prefix + '*') )

    # -- wrapper class headers
    move_files_list += glob.glob( os.path.join(cfg.extra_output_dir, cfg.wrapper_header_prefix + '*') )

    # -- header with forward declarations for all abstract classes
    move_files_list += [ os.path.join(cfg.extra_output_dir, gb.frwd_decls_abs_fname + cfg.header_extension) ]

    # -- header with forward declarations for all wrapper classes
    move_files_list += [ os.path.join(cfg.extra_output_dir, gb.frwd_decls_wrp_fname + cfg.header_extension) ]

    # -- header with copies of all enum type declarations
    move_files_list += [ os.path.join(cfg.extra_output_dir, gb.enum_decls_wrp_fname + cfg.header_extension) ]

    # -- identification.hpp
    move_files_list += [ os.path.join(cfg.extra_output_dir, 'identification.hpp') ]    

    for mv_file in move_files_list:
        shutil.move(mv_file, gb.gambit_backend_dir_complete)

# ====== END: moveFilesAround ========



# ====== createCommonHeaders ========

# Copy header files from common_headers/ to cfg.extra_output_dir
# and replace any code template tags with proper code

def createCommonHeaders():

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


    # - cats.hpp

    source_file_name = 'common_headers/cats.hpp'
    target_file_name = os.path.join(cfg.extra_output_dir, 'cats.hpp')    
    shutil.copyfile(source_file_name, target_file_name)    
    

    # - backend_undefs.hpp

    source_file_name = 'common_headers/backend_undefs.hpp'
    target_file_name = os.path.join(cfg.extra_output_dir, 'backend_undefs.hpp')    
    shutil.copyfile(source_file_name, target_file_name)    


# ====== END: createCommonHeaders ========



# ====== replaceNamespaceTags ========

# Run through a list of files and use code tags open_tag and close_tag
# to construct the namespace constr_namespace.

def replaceNamespaceTags(files_list, constr_namespace, open_tag, close_tag):

    for file_path in files_list:

        if not os.path.isfile(file_path):
            continue

        f = open(file_path, 'r')
        content = f.read()
        f.close()

        new_content = utils.constrNamespaceFromTags(content, constr_namespace, open_tag, close_tag)

        f = open(file_path, 'w')
        f.write(new_content)
        f.close()

# ====== END: replaceNamespaceTags ========



# ====== removeCodeTagsFromFiles ========

# Run through a list of files and remove code tags.

def removeCodeTagsFromFiles(files_list, remove_tags_list):

    for file_path in files_list:

        if not os.path.isfile(file_path):
            continue

        # Read file.
        f = open(file_path, 'r')
        content = f.read()
        f.close()

        # Generate new file content by removing code tags.
        new_content = utils.removeCodeTags(content, remove_tags_list)

        # Write new content to file.
        f = open(file_path,'w')
        f.write(new_content)
        f.close()

# ====== END: removeCodeTagsFromFiles ========



# ====== copyFilesToSourceTree ========

# Copy generated files to original source tree.

def copyFilesToSourceTree(verbose=False):

    # Construct a list of (source,target) touples for all copy operations.
    source_target_touples = []

    # - Add all manipulated original files
    for original_file_short_name, original_file_full_path in gb.original_file_paths.items():

        cp_source = os.path.join(cfg.extra_output_dir,original_file_short_name)
        cp_target = original_file_full_path
        source_target_touples.append( (cp_source, cp_target) )

    # - Add factory source files
    for class_name in gb.classes_done:

        factory_source_fname_short = cfg.factory_file_prefix + class_name['short'] + cfg.source_extension

        cp_source = os.path.join(cfg.extra_output_dir, factory_source_fname_short)
        cp_target = os.path.join(cfg.source_path, factory_source_fname_short)
        source_target_touples.append( (cp_source, cp_target) )

    # - Add 'extras' source files (containing implementations for the helper functions that BOSS adds to the original classes)
    for class_name in gb.classes_done:

        extra_source_fname_short = class_name['short'] + '_extras' + gb.code_suffix + cfg.source_extension

        cp_source = os.path.join(cfg.extra_output_dir, extra_source_fname_short)
        cp_target = os.path.join(cfg.source_path, extra_source_fname_short)
        source_target_touples.append( (cp_source, cp_target) )


    # - Add standard BOSS files that are to be copied to the original source tree

    # -- abstractbase.hpp
    cp_source = os.path.join(cfg.extra_output_dir, 'abstractbase.hpp')
    cp_target = os.path.join(cfg.include_path, gb.gambit_backend_incl_dir, 'abstractbase.hpp')
    source_target_touples.append( (cp_source, cp_target) )

    # -- abstracttypedefs.hpp
    cp_source = os.path.join(cfg.extra_output_dir, gb.abstract_typedefs_fname + cfg.header_extension)
    cp_target = os.path.join(cfg.include_path, gb.gambit_backend_incl_dir, gb.abstract_typedefs_fname + cfg.header_extension)
    source_target_touples.append( (cp_source, cp_target) )

    # -- backend_undefs.hpp
    cp_source = os.path.join(cfg.extra_output_dir, 'backend_undefs.hpp')
    cp_target = os.path.join(cfg.include_path, gb.gambit_backend_incl_dir, 'backend_undefs.hpp')
    source_target_touples.append( (cp_source, cp_target) )

    # -- cats.hpp
    cp_source = os.path.join(cfg.extra_output_dir, 'cats.hpp')
    cp_target = os.path.join(cfg.include_path, gb.gambit_utils_incl_dir, 'cats.hpp')
    source_target_touples.append( (cp_source, cp_target) )

    # -- wrapperbase.hpp
    cp_source = os.path.join(cfg.extra_output_dir, 'wrapperbase.hpp')
    cp_target = os.path.join(cfg.include_path, gb.gambit_backend_incl_dir, 'wrapperbase.hpp')
    source_target_touples.append( (cp_source, cp_target) )

    # -- wrapperdeleter.hpp
    cp_source = os.path.join(cfg.extra_output_dir, gb.wrapper_deleter_fname + cfg.header_extension)
    cp_target = os.path.join(cfg.include_path, gb.gambit_backend_incl_dir, gb.wrapper_deleter_fname + cfg.header_extension)
    source_target_touples.append( (cp_source, cp_target) )

    # -- wrappertypedefs.hpp
    cp_source = os.path.join(cfg.extra_output_dir, gb.wrapper_typedefs_fname + cfg.header_extension)
    cp_target = os.path.join(cfg.include_path, gb.gambit_backend_incl_dir, gb.wrapper_typedefs_fname + cfg.header_extension)
    source_target_touples.append( (cp_source, cp_target) )


    # -- wrapperdeleter.cpp
    cp_source = os.path.join(cfg.extra_output_dir, gb.wrapper_deleter_fname + cfg.source_extension)
    cp_target = os.path.join(cfg.source_path, gb.wrapper_deleter_fname + cfg.source_extension)
    source_target_touples.append( (cp_source, cp_target) )


    # - Add the entire backend_types/ directory to the include directory of the original source tree
    cp_source = os.path.join(cfg.extra_output_dir, gb.gambit_backend_types_basedir)
    cp_target = os.path.join(cfg.include_path, gb.gambit_backend_types_basedir)
    source_target_touples.append( (cp_source, cp_target) )


    # Perform copy operations
    for cp_source, cp_target in source_target_touples:

        if os.path.isfile(cp_source):
            target_dir_name = os.path.dirname(cp_target)
            if not os.path.exists(target_dir_name): 
                os.makedirs(target_dir_name)
            shutil.copyfile(cp_source, cp_target)

        elif os.path.isdir(cp_source):
            shutil.copytree(cp_source, cp_target)

        else:
            continue

        if verbose: 
            n_spaces = max(50-len(cp_source), 2)
            print '   ' + cp_source + ' '*n_spaces  + '--->   ' + cp_target

# ====== END: copyFilesToSourceTree ========



# ====== createLoadedTypesHeader ========

# Generate the header file loaded_types.hpp. This header will 
# contain the symbol names and function signatures for all 
# the generated factory functions.

def createLoadedTypesHeader(factory_xml_files_dict):

    # First update the 'symbol' entry in the dictionaries containing the factory function info
    for class_name in gb.classes_done:

        # Set useful variables
        xml_file = factory_xml_files_dict[class_name['long']]
        info_dicts_list = gb.factory_info[class_name['long']]

        # Get all function elements in the xml file
        factory_func_elements = OrderedDict()
        tree = ET.parse(xml_file)
        root = tree.getroot()
        for func_el in root.findall('Function'):
            func_name = func_el.get('name')
            factory_func_elements[func_name] = func_el

        for info_dict in info_dicts_list:

            factory_el = factory_func_elements[ info_dict['name'] ]
            info_dict['symbol'] = factory_el.get('mangled')
    
    # Generate the code for loaded_types.hpp
    loaded_types_header_content = utils.constrLoadedTypesHeaderContent()

    # Write to file
    loaded_types_output_path = os.path.join(gb.gambit_backend_dir_complete, 'loaded_types.hpp')
    f = open(loaded_types_output_path, 'w')
    f.write(loaded_types_header_content)
    f.close()

# ====== END: createLoadedTypesHeader ========


