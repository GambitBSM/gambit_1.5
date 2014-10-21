

## =========================================================

##
## Added a command line option for creating the loaded_types.hpp header based on the shared library 
##

    parser.add_option("-s", "--symbols-header",
                      action="store",
                      type="string",
                      nargs=1,
                      dest="shared_lib_path",
                      default=None,
                      help="Generate a header 'loaded_types.hpp' containing the symbol names GAMBIT needs from the shared library.")

    # If the -s is given, produce the header file 'loaded_types.hpp' and exit
    if options.shared_lib_path:

        if len(args) > 0:
            parser.error("Unrecognized argument(s) -- The -s option only accepts a single file name.")

        header_code = utils.constrLoadedTypesHeader(options.shared_lib_path)
        header_path = os.path.join(cfg.extra_output_dir, cfg.gambit_backend_basedir, gb.gambit_backend_name_full, 'loaded_types.hpp')

        f = open(header_path,'w')
        f.write(header_code)
        f.close()

        print 'INFO: Header file written to: %s' % header_path
        print
        return


    # If not the -s option is given, proceed as normal with .xml files input
    if len(args) < 1:
        parser.error("No xml file(s) specified.")


## =========================================================

## 
## Add GAMBIT namespace tags to the wrapper typedefs header
## 

    wrapper_typedefs_path = gb.wrapper_typedefs_fname + cfg.header_extension
    if wrapper_typedefs_path not in new_code.keys():
        new_code[wrapper_typedefs_path] = {'code_tuples':[], 'add_include_guard':True}
    new_code[wrapper_typedefs_path]['code_tuples'].append( (0, '__START_GAMBIT_NAMESPACE__\n') )
    new_code[wrapper_typedefs_path]['code_tuples'].append( (-1,'__END_GAMBIT_NAMESPACE__\n') )


## =========================================================

## 
## Add include statement at beginning of wrapper deleter header file
## 

    w_deleter_include = '#include "' + os.path.join(cfg.add_path_to_includes, gb.all_wrapper_fname + '_decl' + cfg.header_extension) + '"\n'
    w_deleter_header_path = os.path.join(cfg.extra_output_dir, gb.wrapper_deleter_fname + cfg.header_extension)
    if w_deleter_header_path not in new_code.keys():
        new_code[w_deleter_header_path] = {'code_tuples':[], 'add_include_guard':True}
    new_code[w_deleter_header_path]['code_tuples'].append( (0, w_deleter_include) )

## =========================================================

##
## Generate 'abstracts' header
##

    abstracts_includes_header_path = os.path.join(cfg.extra_output_dir, 'abstracts_includes_TEMP.hpp')
    abstracts_typedefs_header_path = os.path.join(cfg.extra_output_dir, abstract_typedefs_fname + cfg.header_extension)
    abstracts_header_path          = os.path.join(cfg.extra_output_dir, 'abstracts' + cfg.header_extension)

    f = open(abstracts_includes_header_path, 'r')
    abstracts_includes_header_content = f.read()
    f.close()

    f = open(abstracts_typedefs_header_path, 'r')
    abstracts_typedefs_header_content = f.read()
    f.close()

    abstracts_header_content  = '\n' 
    abstracts_header_content += '// Abstract classes:\n'
    abstracts_header_content += '// -----------------\n'
    abstracts_header_content += '\n' 
    abstracts_header_content += abstracts_includes_header_content
    abstracts_header_content += '\n' 
    abstracts_header_content += '// Typedefs for abstract classes:\n'
    abstracts_header_content += '// ------------------------------\n'
    abstracts_header_content += '\n' 
    abstracts_header_content += '#include "' + os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, 'identification.hpp') + '"\n'
    abstracts_header_content += abstracts_typedefs_header_content
    abstracts_header_content += '#include "backend_undefs.hpp"\n'
    abstracts_header_content += '\n' 


    # Add include guards
    abstracts_header_content = utils.addIncludeGuard(abstracts_header_content, os.path.basename(abstracts_header_path))

    # Write the 'abstracts' header (if not in debug_mode)
    if not options.debug_mode_flag:
        f = open(abstracts_header_path, 'w')
        f.write(abstracts_header_content)
        f.close()

    # Remove the temporary headers 'abstracts_includes_TEMP.hpp' and 'abstracttypedefs.hpp' we used to construct the master header
    os.remove(abstracts_includes_header_path)
    # os.remove(abstracts_typedefs_header_path)

## =========================================================

##
## Generate 'wrappers' header
##

    wrappers_decl_includes_header_path = os.path.join(cfg.extra_output_dir, 'wrappers_decl_includes_TEMP.hpp')
    wrappers_def_includes_header_path  = os.path.join(cfg.extra_output_dir, 'wrappers_def_includes_TEMP.hpp')
    wrappers_typedefs_header_path      = os.path.join(cfg.extra_output_dir, 'wrappers_typedefs_TEMP.hpp')
    wrappers_header_path               = os.path.join(cfg.extra_output_dir, 'wrappers' + cfg.header_extension)

    f = open(wrappers_decl_includes_header_path, 'r')
    wrappers_decl_includes_header_content = f.read()
    f.close()

    f = open(wrappers_def_includes_header_path, 'r')
    wrappers_def_includes_header_content = f.read()
    f.close()

    f = open(wrappers_typedefs_header_path, 'r')
    wrappers_typedefs_header_content = f.read()
    f.close()

    wrappers_header_content  = '\n' 
    wrappers_header_content += '// Wrapper classes:\n'
    wrappers_header_content += '// ----------------\n'
    wrappers_header_content += '\n' 
    wrappers_header_content += wrappers_decl_includes_header_content
    wrappers_header_content += '\n' 
    wrappers_header_content += wrappers_def_includes_header_content
    wrappers_header_content += '\n' 
    wrappers_header_content += '// Typedefs for wrapper classes:\n'
    wrappers_header_content += '// -----------------------------\n'
    wrappers_header_content += '\n' 
    wrappers_header_content += '#include "' + os.path.join(cfg.gambit_backend_basedir, gb.gambit_backend_name_full, 'identification.hpp') + '"\n'
    wrappers_header_content += wrappers_typedefs_header_content
    wrappers_header_content += '#include "backend_undefs.hpp"\n'
    wrappers_header_content += '\n' 


    # Add include guards
    wrappers_header_content = utils.addIncludeGuard(wrappers_header_content, os.path.basename(wrappers_header_path))

    # Write the 'wrappers' header (if not in debug_mode)
    if not options.debug_mode_flag:
        f = open(wrappers_header_path, 'w')
        f.write(wrappers_header_content)
        f.close()

    # Remove the temporary headers 'wrappers_decl_includes_TEMP.hpp', 'wrappers_def_includes_TEMP.hpp' and 'wrappers_typedefs_TEMP.hpp' 
    # we used to construct the master header
    os.remove(wrappers_decl_includes_header_path)
    os.remove(wrappers_def_includes_header_path)
    os.remove(wrappers_typedefs_header_path)

## =========================================================

##
## Generate the master header
##

    all_wrapper_decl_header_path = os.path.join(cfg.extra_output_dir, gb.all_wrapper_fname + '_decl' + cfg.header_extension)
    all_wrapper_def_header_path  = os.path.join(cfg.extra_output_dir, gb.all_wrapper_fname + '_def' + cfg.header_extension)
    master_header_path           = os.path.join(cfg.extra_output_dir, gb.all_wrapper_fname + cfg.header_extension)

    f = open(all_wrapper_decl_header_path, 'r')
    all_wrapper_decl_header_content = f.read()
    f.close()

    f = open(all_wrapper_def_header_path, 'r')
    all_wrapper_def_header_content = f.read()
    f.close()

    master_header_content  = '\n' 
    master_header_content += '// Loaded types declarations:\n'
    master_header_content += '// --------------------------\n'
    master_header_content += '\n' 
    master_header_content += all_wrapper_decl_header_content
    master_header_content += '\n' 
    master_header_content += '// Loaded types implementations:\n'
    master_header_content += '// -----------------------------\n'
    master_header_content += '\n' 
    master_header_content += all_wrapper_def_header_content
    master_header_content += '\n' 
    master_header_content += '// Typedefs: wrapper class name --> original class name:\n'
    master_header_content += '// -----------------------------------------------------\n'
    master_header_content += '\n' 
    master_header_content += '#include "' + gb.wrapper_typedefs_fname + cfg.header_extension + '"\n' 


    # Add include guards
    master_header_content = utils.addIncludeGuard(master_header_content, gb.all_wrapper_fname)

    # Write the master header (if not in debug_mode)
    if not options.debug_mode_flag:
        f = open(master_header_path, 'w')
        f.write(master_header_content)
        f.close()

    # # Remove the temporary '_decl' and '_def' header files we used to construct the master header
    # os.remove(all_wrapper_decl_header_path)
    # os.remove(all_wrapper_def_header_path)

## =========================================================

## 
## 
## 

## =========================================================

## 
## 
## 

## =========================================================

## 
## 
## 

## =========================================================

## 
## 
## 



