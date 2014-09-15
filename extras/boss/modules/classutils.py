####################################
#                                  #
#  Utility functions for handling  # 
#  C++ classes with BOSS           #
#                                  #
####################################

import xml.etree.ElementTree as ET
from collections import OrderedDict
import warnings
import os

import modules.cfg as cfg
import modules.funcutils as funcutils
import modules.utils as utils


# ====== getAbstractClassName ========

def getAbstractClassName(input_name, prefix=cfg.abstr_class_prefix, short=False):

    if '::' in input_name:
        namespaces, short_class_name = input_name.rsplit('::',1)
        abstract_class_name = namespaces + '::' + cfg.abstr_class_prefix + short_class_name
    else:
        abstract_class_name = cfg.abstr_class_prefix + input_name

    if short == True:
        return abstract_class_name.rsplit('::',1)[-1]
    else:
        return abstract_class_name

# ====== END: getAbstractClassName ========


# ====== constrEmptyTemplClassDecl ========

def constrEmptyTemplClassDecl(abstr_class_name_short, namespaces, template_bracket, indent=4):

    n_indents  = len(namespaces)
    class_decl = ''

    class_decl += utils.constrNamespace(namespaces, 'open')

    class_decl += ' '*n_indents*indent + 'template ' + template_bracket + '\n'
    class_decl += ' '*n_indents*indent + 'class ' + abstr_class_name_short + ' {};\n'

    class_decl += utils.constrNamespace(namespaces, 'close')
    
    class_decl += '\n'

    return class_decl

# ====== END: constrEmptyTemplClassDecl ========



# ====== constrTemplForwDecl ========

def constrTemplForwDecl(class_name_short, namespaces, template_bracket, indent=4):

    n_indents = len(namespaces)
    forw_decl = ''

    forw_decl += utils.constrNamespace(namespaces, 'open')

    forw_decl += ' '*n_indents*indent + 'template ' + template_bracket + '\n'
    forw_decl += ' '*n_indents*indent + 'class ' + class_name_short + ';\n'

    forw_decl += utils.constrNamespace(namespaces, 'close')
    
    forw_decl += '\n'

    return forw_decl

# ====== END: constrTemplForwDecl ========



# ====== constrAbstractClassDecl ========

def constrAbstractClassDecl(class_el, class_name_short, abstr_class_name_short, namespaces, indent=4, template_types=[], 
                            has_copy_constructor=True,  construct_assignment_operator=True):

    n_indents = len(namespaces)

    # Check template_types argument:
    if len(template_types) > 0:
        is_template = True
        print '--> This is a template class'

    else:
        is_template = False


    # Create list of all 'non-artificial' members of the class
    member_elements = []
    if 'members' in class_el.keys():
        for mem_id in class_el.get('members').split():
            el = cfg.id_dict[mem_id]
            if not 'artificial' in el.keys():
                member_elements.append(el)

    # # Create list of member functions that needs overloading due to default argument values
    # overload_members = []
    # for i in range(len(member_elements)):
    #     if ( member_elements[i].tag == 'Method' ) and ( funcutils.numberOfDefaultArgs(member_elements[i])>0 ):
    #         overload_members.append(member_elements[i])


    # Get list of dicts with info on parent classes
    parent_classes = utils.getParentClasses(class_el)

    # Get wrapper class name (short)
    w_class_name = toWrapperType(class_name_short)


    #
    # Construct the abstract class declaration
    #
    
    class_decl = ''

    # class_decl = '#pragma GCC diagnostic push\n'
    # class_decl += '#pragma GCC diagnostic ignored "-Wunused-parameter"\n'
    # class_decl += '#pragma GCC diagnostic ignored "-Wreturn-type"\n'


    # - Add forward declarations needed by the 'destructor pattern'
    class_decl += '// Forward declarations needed by the destructor pattern.\n'
    class_decl += 'class ' + w_class_name + ';\n'
    class_decl += 'void wrapper_deleter(' + w_class_name + '*);\n'
    class_decl += '\n'

    # - Construct the beginning of the namespaces
    class_decl += utils.constrNamespace(namespaces, 'open')

    # - If this class is a template specialization, add 'template <>' at the top
    if is_template == True:
        class_decl += ' '*n_indents*indent + 'template <>\n'

    # - Construct the declaration line, with inheritance of abstract classes
    inheritance_line = ' : virtual public AbstractBase, '
    for parent_dict in parent_classes:

        if parent_dict['loaded']:
            inheritance_line += 'virtual ' + parent_dict['access'] + ' ' + parent_dict['abstr_class_name']['long_templ'] + ', '

        elif parent_dict['fundamental'] or parent_dict['std']:
            # inheritance_line += 'virtual ' + parent_dict['access'] + ' ' + parent_dict['class_name']['long_templ'] + ', '
            warnings.warn("In class '%s', the parent class '%s' is ignored (to avoid inheritance ambiguity)." % (abstr_class_name_short, parent_dict['class_name']['long_templ']))

        else:
            warnings.warn("In class '%s', the parent class '%s' is ignored (not loaded class or std/fundamental type)." % (abstr_class_name_short, parent_dict['class_name']['long_templ']))
            continue

    inheritance_line = inheritance_line.rstrip(', ')

    class_decl += ' '*n_indents*indent
    if is_template:
        class_decl += 'class ' + abstr_class_name_short + '<' + ','.join(template_types) + '>' + inheritance_line + '\n'
    else:
        class_decl += 'class ' + abstr_class_name_short + inheritance_line + '\n'

    # - Construct body of class declaration
    current_access = ''
    class_decl += ' '*n_indents*indent
    class_decl += '{' + '\n'
    done_members = []
    for el in member_elements:

        # Check access
        element_access = el.get('access')
        if current_access != element_access:
            class_decl += ' '*(n_indents+1)*indent
            class_decl += element_access + ':' +'\n'
            current_access = element_access

        #
        # Add code based on what element type this is
        #
        if el.tag in ['Constructor', 'Destructor']:
            pass   # (An empty virtual destructor will be added later)

        elif el.tag in ['Method', 'OperatorMethod']:

            # Check if this is an operator function
            is_operator = False
            if el.tag == 'OperatorMethod':
                is_operator = True

            # Check if this member function should be ignored based on unloaded types.
            if funcutils.ignoreFunction(el):
                print 'INFO: The member function "' + is_operator*'operator' + el.get('name') + '" is ignored due to non-accepted type(s).'
                continue

            # Check if this member function makes use of loaded types
            uses_loaded_type = funcutils.usesLoadedType(el)


            return_type, return_kw, return_id = utils.findType( cfg.id_dict[el.get('returns')] )
            return_type_base = return_type.replace('*','').replace('&','')
            
            return_kw_str = ' '.join(return_kw)        
            return_kw_str += ' '*bool(len(return_kw))

            return_el = cfg.id_dict[return_id]
            return_is_loaded = utils.isLoadedClass(return_el)
            args = funcutils.getArgs(el)
            w_args = funcutils.constrWrapperArgs(args, add_ref=True)


            # Check constness
            if ('const' in el.keys()) and (el.get('const')=='1'):
                is_const = True
            else:
                is_const = False

            # # Check if we're generating an overload of a previous member function.
            # overload_n = done_members.count(el)
            # is_overload = bool(overload_n)

            # # If so, remove some arguments with default values
            # if is_overload:
            #     args = args[:-overload_n]

            # if (return_type == 'void') or (return_type.count('*') > 0):
            #     w_return_type = return_type
            # elif return_is_loaded:
            #     w_return_type = return_type.replace(return_type_base, return_type_base+'*')
            # else:
            #     w_return_type = return_type

            # If default arguments are used, we need to generate overloaded versions
            n_overloads = funcutils.numberOfDefaultArgs(el)

            # One overloaded version for each set of default arguments
            for remove_n_args in range(n_overloads+1):

                if remove_n_args == 0:
                    use_args   = args
                    use_w_args = w_args
                else:
                    use_args   = args[:-remove_n_args]
                    use_w_args = w_args[:-remove_n_args]


                w_args_bracket = funcutils.constrArgsBracket(use_w_args, include_arg_name=True, include_arg_type=True, include_namespace=True)
                w_args_bracket_notypes = funcutils.constrArgsBracket(use_w_args, include_arg_name=True, include_arg_type=False)
                w_args_bracket_nonames = funcutils.constrArgsBracket(use_w_args, include_arg_name=False, include_arg_type=True, include_namespace=True)

                
                if is_operator:
                    if uses_loaded_type:
                        w_func_name = 'operator_' + cfg.operator_names[el.get('name')] + cfg.code_suffix
                    else:
                        w_func_name = 'operator' + el.get('name')    
                else:
                    # w_func_name = el.get('name') + cfg.code_suffix 
                    if uses_loaded_type or (remove_n_args>0):
                        w_func_name = el.get('name') + cfg.code_suffix 
                    else:
                        w_func_name = el.get('name')



                    # if is_operator:
                    #     w_func_name = 'operator_' + cfg.operator_names[el.get('name')] + cfg.code_suffix
                    # else:
                    #     w_func_name = el.get('name') + cfg.code_suffix 

                    # if remove_n_args > 0:
                    #     w_func_name += '_overload_' + str(remove_n_args)

                #
                # If the method makes use of a loaded class, construct a pair of wrapper methods.
                #
                if uses_loaded_type:

                    # Construct the virtual member function that is overridden, e.g.:  
                    #
                    #   virtual X* getX_GAMBIT(arguments) {}
                    #

                    if return_is_loaded:
                        if return_type.count('*') > 0:
                            w_return_type = toAbstractType(return_type, include_namespace=True)
                        else:
                            w_return_type = toAbstractType(return_type, include_namespace=True, add_pointer=True, remove_reference=True)
                    else:
                        w_return_type = return_type

                    class_decl += '\n'
                    class_decl += ' '*(n_indents+2)*indent
                    # class_decl += 'virtual ' + return_kw_str + w_return_type + ' ' + w_func_name + w_args_bracket_nonames + is_const*' const' + ' {std::cout << "Called virtual function" << std::endl;};' + '\n'
                    class_decl += 'virtual ' + return_kw_str + w_return_type + ' ' + w_func_name + w_args_bracket_nonames + is_const*' const' + ' =0;' + '\n'


                    # # Construct the member function with the original name, wrapping the overridden one, e.g.: 
                    # #
                    # #   Abstract__X* getX(arguments)
                    # #   {
                    # #     return reinterpret_cast<Abstract__X*>(getX_GAMBIT(arguments));
                    # #   }
                    # #

                    # if is_operator:
                    #     w2_func_name = 'operator' + el.get('name')
                    # else:
                    #     w2_func_name = el.get('name')
                    # w2_args_bracket = w_args_bracket

                    # class_decl += ' '*(n_indents+2)*indent 
                    # class_decl += return_kw_str + w_return_type + ' ' + w2_func_name + w2_args_bracket + is_const*' const' +'\n'
                    # class_decl += ' '*(n_indents+2)*indent + '{\n'
                    # if (return_type == 'void'):
                    #     class_decl += ' '*(n_indents+3)*indent + w_func_name + w_args_bracket_notypes + ';\n'
                    # else:
                    #     # OLD:
                    #     # if return_is_loaded:
                    #     #     class_decl += ' '*(n_indents+3)*indent + 'return reinterpret_cast<' + return_kw_str + w_return_type + '>(' + w_func_name + w_args_bracket_notypes + ');\n'
                    #     # else:
                    #     #     class_decl += ' '*(n_indents+3)*indent + 'return ' + w_func_name + w_args_bracket_notypes + ';\n'
                    #     # NEW:
                    #     class_decl += ' '*(n_indents+3)*indent + 'return ' + w_func_name + w_args_bracket_notypes + ';\n'

                    # class_decl += ' '*(n_indents+2)*indent + '}\n'

                #
                # If the method does not make use of any loaded class, construct a single virtual method
                #
                else:
                    class_decl += '\n'
                    class_decl += ' '*(n_indents+2)*indent
                    # class_decl += 'virtual ' + return_kw_str + return_type + ' ' + w_func_name + w_args_bracket_nonames + is_const*' const' + ' {std::cout << "Called virtual function" << std::endl;};' + '\n'
                    class_decl += 'virtual ' + return_kw_str + return_type + ' ' + w_func_name + w_args_bracket_nonames + is_const*' const' + ' =0;' + '\n'


        #
        # If element is a public member variable of accepted type, construct virtual method that returns a reference to this variable
        #
        elif (el.tag in ('Field', 'Variable')) and (el.get('access') == 'public') and utils.isAcceptedType(el):

            class_decl += '\n' 
            class_decl += constrVariableRefFunction(el, virtual=True, indent=indent, n_indents=n_indents+2)

        else:
            class_decl += ' '*(n_indents+2)*indent
            if 'name' in el.keys():
                class_decl += '// IGNORED: ' + el.tag + '  -- Name: ' + el.get('name') + '  -- XML id: ' + el.get('id') + '\n'
            else:
                class_decl += '// IGNORED: ' + el.tag + '  -- XML id: ' + el.get('id') + '\n'

    

    # - Construct 'pointerAssign' and 'pointerCopy' functions -- UPDATE: Always construct these functions
    # if has_copy_constructor or construct_assignment_operator:
        # class_decl += '\n'
        # class_decl += ' '*(n_indents+1)*indent + 'public:\n'
        # if construct_assignment_operator:
        #     class_decl += constrPtrAssignFunc(abstr_class_name_short, class_name_short, virtual=True, indent=indent, n_indents=n_indents+2)
        # if has_copy_constructor:
        #     class_decl += constrPtrCopyFunc(abstr_class_name_short, class_name_short, virtual=True, indent=indent, n_indents=n_indents+2)

    class_decl += '\n'
    class_decl += ' '*(n_indents+1)*indent + 'public:\n'
    class_decl += constrPtrAssignFunc(class_el, abstr_class_name_short, class_name_short, virtual=True, indent=indent, n_indents=n_indents+2)
    class_decl += constrPtrCopyFunc(class_el, abstr_class_name_short, class_name_short, virtual=True, indent=indent, n_indents=n_indents+2)



    # # - Construct the 'downcast' converter function
    # class_decl += '\n'
    # # class_decl += ' '*(n_indents+1)*indent + 'NEW CODE HERE!\n'
    # if is_template:
    #     template_bracket = '<' + ','.join(template_types) + '>'
    # else:
    #     template_bracket = ''
    # class_decl += ' '*(n_indents+1)*indent + 'public:\n'
    # class_decl += constrDowncastFunction(class_name_short, indent=indent, n_indents=n_indents+2, template_bracket=template_bracket)

    # # - Construct an empty virtual destructor
    # class_decl += ' '*(n_indents+2)*indent
    # class_decl += 'virtual ~' + abstr_class_name_short + '() {};\n'


    # - Construct code needed for 'destructor pattern' (abstract class and wrapper class must can delete each other)
    class_decl += '\n'
    class_decl += ' '*(n_indents+1)*indent + 'private:\n'
    class_decl += ' '*(n_indents+2)*indent + w_class_name + '* wptr;\n'

    class_decl += '\n'
    class_decl += ' '*(n_indents+1)*indent + 'public:\n'
    class_decl += ' '*(n_indents+2)*indent + 'void wrapper' + cfg.code_suffix + '(' + w_class_name + '* wptr_in)\n'
    class_decl += ' '*(n_indents+2)*indent + '{\n'
    class_decl += ' '*(n_indents+3)*indent + 'wptr = wptr_in;\n'
    class_decl += ' '*(n_indents+3)*indent + 'is_wrapped(true);\n'
    class_decl += ' '*(n_indents+2)*indent + '}\n'
    class_decl += '\n'
    class_decl += ' '*(n_indents+2)*indent + w_class_name + '* wrapper' + cfg.code_suffix + '()\n'
    class_decl += ' '*(n_indents+2)*indent + '{\n'
    class_decl += ' '*(n_indents+3)*indent + 'return wptr;\n'
    class_decl += ' '*(n_indents+2)*indent + '}\n'
    class_decl += '\n'
    class_decl += ' '*(n_indents+2)*indent + 'virtual ~' + abstr_class_name_short + '()\n' 
    class_decl += ' '*(n_indents+2)*indent + '{\n'
    class_decl += ' '*(n_indents+3)*indent + 'if (can_delete_wrapper())\n' 
    class_decl += ' '*(n_indents+3)*indent + '{\n'
    class_decl += ' '*(n_indents+4)*indent + 'can_delete_me(false);\n'
    class_decl += ' '*(n_indents+4)*indent + 'wrapper_deleter(wptr);\n'
    class_decl += ' '*(n_indents+3)*indent + '}\n'
    class_decl += ' '*(n_indents+2)*indent + '}\n'

    # - Close the class body
    class_decl += ' '*n_indents*indent + '};' + '\n'

    # - Construct the closing of the namespaces
    class_decl += utils.constrNamespace(namespaces, 'close')

    # class_decl += '#pragma GCC diagnostic pop\n'
    return class_decl

# ====== END: constrAbstractClassDecl ========



# ====== constrFactoryFunction ========

def constrFactoryFunction(class_el, class_name, indent=4, template_types=[], skip_copy_constructors=False, use_wrapper_class=False, add_include_statements=True):

    # Replace '*' and '&' in list of template types
    template_types = [e.replace('*','P').replace('&','R') for e in template_types]

    # Check for copy constructor. (We do not generate factory functions for copy constructors)
    if skip_copy_constructors:
        has_copy_constructor, copy_constr_id = checkCopyConstructor(class_el, return_id=True)

    # Create list of all constructors of the class
    constructor_elements = []
    if 'members' in class_el.keys():
        for mem_id in class_el.get('members').split():
            el = cfg.id_dict[mem_id]
            if (el.tag == 'Constructor') and (el.get('access') == 'public'): #and ('artificial' not in el.keys()):  #(el.get('explicit') == "1"):
                if skip_copy_constructors and (el.get('id') == copy_constr_id):
                    pass
                else:
                    constructor_elements.append(el)

    # List to hold include statements that are generated based on the types used
    # in the constructors
    if add_include_statements:
        include_statements = []

    # Construct factory function definition(s)
    func_def = ''
    for el in constructor_elements:

        if add_include_statements:
            # - Generate include statements based on the types used in the constructor
            include_statements += utils.getIncludeStatements(el, convert_loaded_to='none', input_element='function', add_extra_include_path=True)
            include_statements += utils.getIncludeStatements(el, convert_loaded_to='wrapper', input_element='function', add_extra_include_path=True)

        # Useful variables
        factory_name = 'Factory_' + class_name['short']
        if len(template_types) > 0:
            factory_name += '_' + '_'.join(template_types)

        # We need to generate as many overloaded versions as there are arguments with default values
        n_overloads = funcutils.numberOfDefaultArgs(el)

        # Identify arguments
        args = funcutils.getArgs(el)

        # Translate argument type of loaded classes
        if use_wrapper_class:
            w_args = funcutils.constrWrapperArgs(args, add_ref=True, convert_loaded_to_abstract=False)
        else:
            w_args = funcutils.constrWrapperArgs(args, add_ref=True, convert_loaded_to_abstract=True)

        # Invent argument names if missing
        argc = 1
        for i in range(len(args)):
            if args[i]['name'] == '':
                args[i]['name'] = 'arg_' + str(argc)
                argc += 1

        # Generate one factory function for each set of default arguments
        for remove_n_args in range(n_overloads+1):

            if remove_n_args == 0:
                use_args   = args
                use_w_args = w_args
            else:
                use_args   = args[:-remove_n_args]
                use_w_args = w_args[:-remove_n_args]

            # Construct bracket with input arguments
            if use_wrapper_class:
                args_bracket         = funcutils.constrArgsBracket(use_w_args, include_namespace=True, use_wrapper_class=True)
                args_bracket_notypes = funcutils.constrArgsBracket(use_args, include_arg_type=False, cast_to_original=True, wrapper_to_pointer=True)
            else:
                args_bracket         = funcutils.constrArgsBracket(use_w_args, include_namespace=True)
                args_bracket_notypes = funcutils.constrArgsBracket(use_args, include_arg_type=False, cast_to_original=True)
            
            # Generate declaration line:
            if use_wrapper_class:
                return_type = toWrapperType(class_name['long'])
            else:
                return_type = toAbstractType(class_name['long'], add_pointer=True)
            func_def += return_type + ' ' + factory_name + args_bracket + '\n'

            # Generate body
            func_def += '{' + '\n'
            if use_wrapper_class:
                func_def += indent*' ' + 'return ' + return_type + '( new ' + class_name['long'] + args_bracket_notypes + ' );' + '\n'
            else:
                func_def += indent*' ' + 'return new ' + class_name['long'] + args_bracket_notypes + ';' + '\n'
            func_def += '}' + 2*'\n'


    if add_include_statements:
        include_statements = list( set(include_statements) )
        include_statements_code = '\n'.join(include_statements) + 2*'\n'
        func_def = include_statements_code + func_def

    return func_def

# ====== END: constrFactoryFunction ========



# # ====== constrDowncastFunction ========

# def constrDowncastFunction(cast_to_class, indent=4, n_indents=0, template_bracket=''):

#     # Example of output C++ code:
#     #
#     #     public: A* downcast() 
#     #     {
#     #         return reinterpret_cast<A*>(this); 
#     #     }
#     #

#     base_indent = n_indents*indent
#     use_class_name = cast_to_class + template_bracket

#     func_def  = ''
#     func_def += base_indent*' ' + use_class_name + '* downcast()\n'
#     func_def += base_indent*' ' + '{\n'
#     func_def += (base_indent+indent)*' ' + 'return reinterpret_cast<' + use_class_name +'*>(this);\n'
#     func_def += base_indent*' ' + '}\n'

#     return func_def

# # ====== END: constrDowncastFunction ========



# ====== constrUpcastFunction ========

# def constrUpcastFunction(cast_to_class, indent=4, n_indents=0, template_bracket=''):

#     # Example of output C++ code:
#     #
#     #     public: Abstract__A* upcast() 
#     #     { 
#     #         return dynamic_cast<Abstract__A*>(this); 
#     #     }
#     #

#     base_indent = n_indents*indent
#     use_class_name = cast_to_class + template_bracket

#     func_def  = ''
#     func_def += base_indent*' ' + use_class_name + '* upcast()\n'
#     func_def += base_indent*' ' + '{\n'
#     func_def += (base_indent+indent)*' ' + 'return dynamic_cast<' + use_class_name +'*>(this);\n'
#     func_def += base_indent*' ' + '}\n'

#     return func_def

# ====== END: constrUpcastFunction ========



# ====== constrAbstractifyFunction ========

# def constrAbstractifyFunction(cast_to_class, indent=4, n_indents=0, template_bracket=''):

#     # Example of output C++ code:
#     #
#     #     public: Abstract__A abstractify() 
#     #     {
#     #         Abstract__A* ptr = reinterpret_cast<Abstract__A*>(this);
#     #         return *ptr;
#     #     }
#     #

#     base_indent = n_indents*indent
#     use_class_name = cast_to_class + template_bracket

#     func_def  = ''
#     func_def += base_indent*' ' + use_class_name + ' abstractify()\n'
#     func_def += base_indent*' ' + '{\n'
#     func_def += (base_indent+indent)*' ' + use_class_name +'* ptr = reinterpret_cast<' + use_class_name +'*>(this);\n'
#     func_def += (base_indent+indent)*' ' + 'return *ptr;\n'
#     func_def += base_indent*' ' + '}\n'

#     return func_def

# ====== END: constrAbstractifyFunction ========



# ====== constrWrapperFunction ========

def constrWrapperFunction(method_el, indent=cfg.indent, n_indents=0, remove_n_args=0, only_declaration=False, include_full_namespace=False):

    # Check if this is an operator function
    is_operator = False
    if method_el.tag == 'OperatorMethod':
        is_operator = True

    # If operator, check that we have a name for it
    if (is_operator) and (method_el.get('name') not in cfg.operator_names.keys()):
        raise Exception('No known name for the operator: %s  -- Add an entry to the following dictionary: cfg.operator_names' % method_el.get('name'))

    # Function name
    if is_operator:
        func_name = 'operator' + method_el.get('name')
    else:
        func_name = method_el.get('name')


    # Function return type
    return_type, return_kw, return_id = utils.findType( cfg.id_dict[method_el.get('returns')] )
    return_el = cfg.id_dict[return_id]

    return_is_loaded_class = utils.isLoadedClass(return_el)
    pointerness, is_ref = utils.pointerAndRefCheck(return_type, byname=True)

    # Function arguments (get list of dicts with argument info)
    args = funcutils.getArgs(method_el)

    # Remove arguments when creating overloaded versions (for dealing with default argument values)
    if remove_n_args > 0:
        args = args[:-remove_n_args]

    # Check constness
    if ('const' in method_el.keys()) and (method_el.get('const')=='1'):
        is_const = True
    else:
        is_const = False

    # Construct wrapper function name
    w_func_name = funcutils.constrWrapperName(method_el, include_full_namespace=include_full_namespace)
    # if remove_n_args > 0:
    #     w_func_name += '_overload_' + str(remove_n_args)


    # Choose wrapper return type
    if return_is_loaded_class:
        if pointerness == 0:
            w_return_type = toAbstractType(return_type, include_namespace=True, add_pointer=True, remove_reference=True)
        else:
            w_return_type = toAbstractType(return_type, include_namespace=True)

    else:
        w_return_type = return_type


    # return_type_base = return_type.replace('*','').replace('&','')
    # if (return_type == 'void') or (return_type.count('*') > 0):
    #     w_return_type = return_type
    # else:
    #     w_return_type = return_type.replace(return_type_base, return_type_base+'*')

    # if return_is_native:
    #     w_return_type = cfg.abstr_class_prefix + return_type
    # else:
    #     w_return_type = return_type

    
    # Construct list of arguments for wrapper function
    w_args = funcutils.constrWrapperArgs(args, add_ref=True)

    # Construct bracket with input arguments for wrapper function
    if only_declaration:
        w_args_bracket = funcutils.constrArgsBracket(w_args, include_arg_name=False)
    else:
        w_args_bracket = funcutils.constrArgsBracket(w_args)

    # Construct declaration line for wrapper function
    w_func_line = funcutils.constrDeclLine(w_return_type, w_func_name, w_args_bracket, keywords=return_kw, is_const=is_const)

    # Construct function body for wrapper function
    if only_declaration:
        pass
    else:
        w_func_body = funcutils.constrWrapperBody(return_type, func_name, args, return_is_loaded_class, keywords=return_kw)

    # Combine code and add indentation
    wrapper_code  = ''
    if only_declaration:
        wrapper_code += utils.addIndentation(w_func_line, n_indents*indent) + ';\n'
    else:
        wrapper_code += utils.addIndentation(w_func_line, n_indents*indent) + '\n'
        wrapper_code += utils.addIndentation(w_func_body, n_indents*indent) + '\n'

    # Return result
    return wrapper_code

# ====== END: constrWrapperFunction ========



# ====== constrVariableRefFunction ========

def constrVariableRefFunction(var_el, virtual=False, indent=cfg.indent, n_indents=0, only_declaration=False, include_full_namespace=False):

    func_code = ''

    var_name = var_el.get('name')
    ref_method_name = var_name + '_ref' + cfg.code_suffix

    if include_full_namespace:
        namespaces = utils.getNamespaces(var_el)
        if len(namespaces) > 0:
            ref_method_name = '::'.join(namespaces) + '::' + ref_method_name

    var_type, var_kw, var_id = utils.findType( cfg.id_dict[var_el.get('type')] )
    var_type_base = var_type.replace('*','').replace('&','')
    
    var_kw_str = ' '.join(var_kw)        
    var_kw_str += ' '*bool(len(var_kw))

    var_el = cfg.id_dict[var_id]
    var_is_loaded_class = utils.isLoadedClass(var_el)
    if var_is_loaded_class:
        if include_full_namespace:
            return_type = getAbstractClassName(var_type_base, prefix=cfg.abstr_class_prefix, short=False)
        else:
            return_type = getAbstractClassName(var_type_base, prefix=cfg.abstr_class_prefix, short=True)
    else:
        return_type = var_type_base

    func_code += ' '*n_indents*indent
    if virtual:
        # func_code += 'virtual ' + var_kw_str + return_type + '& ' + ref_method_name + '() {std::cout << "Called virtual function" << std::endl;};\n'
        func_code += 'virtual ' + var_kw_str + return_type + '& ' + ref_method_name + '() =0;\n'
    else:
        func_code += var_kw_str + return_type + '& ' + ref_method_name + '()' 
        if only_declaration:
            func_code += ';\n'
        else:
            func_code += ' { return ' + var_name  +'; }\n'

    return func_code

# ====== END: constrVariableRefFunction ========



# ====== constrPtrCopyFunc ========

def constrPtrCopyFunc(class_el, abstr_class_name_short, class_name_short, virtual=False, indent=cfg.indent, n_indents=0, only_declaration=False, include_full_namespace=False):

    func_name = 'pointerCopy' + cfg.code_suffix
    class_name = class_name_short
    abstr_class_name = abstr_class_name_short

    if include_full_namespace:
        namespaces_with_self = utils.getNamespaces(class_el, include_self=True)
        namespaces           = utils.getNamespaces(class_el)
        if len(namespaces_with_self) > 0:
            func_name = '::'.join(namespaces_with_self) + '::' + func_name
        if len(namespaces) > 0:
            abstr_class_name = '::'.join(namespaces) + '::' + abstr_class_name
            class_name = '::'.join(namespaces) + '::' + class_name

    ptr_code = ''
    ptr_code += ' '*cfg.indent*n_indents
   
    if virtual:
        ptr_code += 'virtual '+ abstr_class_name + '*' + ' ' + func_name + '() =0;\n'   
    else:
        ptr_code += abstr_class_name + '*' + ' ' + func_name + '()'
        if only_declaration:
            ptr_code += ';\n'
        else:
            ptr_code += ' ' + '{ return new ' + class_name + '(*this); }\n'

    return ptr_code

# ====== END: constrPtrCopyFunc ========



# ====== constrPtrAssignFunc ========

def constrPtrAssignFunc(class_el, abstr_class_name_short, class_name_short, virtual=False, indent=cfg.indent, n_indents=0, only_declaration=False, include_full_namespace=False):

    func_name  = 'pointerAssign' + cfg.code_suffix
    class_name = class_name_short
    abstr_class_name = abstr_class_name_short

    if include_full_namespace:
        namespaces_with_self = utils.getNamespaces(class_el, include_self=True)
        namespaces           = utils.getNamespaces(class_el)
        if len(namespaces_with_self) > 0:
            func_name = '::'.join(namespaces_with_self) + '::' + func_name
        if len(namespaces) > 0:
            abstr_class_name = '::'.join(namespaces) + '::' + abstr_class_name
            class_name = '::'.join(namespaces) + '::' + class_name

    ptr_code = ''
    ptr_code += ' '*cfg.indent*n_indents
    
    if virtual:
        ptr_code += 'virtual void ' + func_name + '(' + abstr_class_name + '*) =0;\n'
    else:
        ptr_code += 'void ' + func_name + '(' + abstr_class_name + '* in)'
        if only_declaration:
            ptr_code += ';\n'
        else:
            ptr_code += ' ' + '{ *this = *dynamic_cast<' + class_name_short + '*>(in); }\n'        

    return ptr_code

# ====== END: constrPtrAssignFunc ========



# ====== checkAssignmentOperator ========

def checkAssignmentOperator(class_el):

    found_assignment_operator = False
    is_artificial = False

    # Get list of all class members
    class_members = utils.getMemberElements(class_el, include_artificial=True)

    # Look for assignment operator
    for mem_el in class_members:
        if (mem_el.tag == 'OperatorMethod') and (mem_el.get('name') == '='):

            # Check that return type is either void or the class type itself (possibly as a reference)
            return_type, return_type_kw, return_type_id = utils.findType(mem_el)
            if (return_type == 'void') or (return_type_id == class_el.get('id')):

                # Check that the only argument is another class instance
                args = funcutils.getArgs(mem_el)
                if (len(args) == 1) and (args[0]['id'] == class_el.get('id')):

                    found_assignment_operator = True
    
                    if ('artificial' in mem_el.keys()) and (mem_el.get('artificial')=='1'):
                        is_artificial = True

    return found_assignment_operator, is_artificial

# ====== END: checkAssignmentOperator ========



# ====== checkCopyConstructor ========

def checkCopyConstructor(class_el, return_id=False):

    found_copy_constructor = False
    copy_constr_id = ''

    # Get list of all class members
    class_members = utils.getMemberElements(class_el, include_artificial=True)

    # Look for copy constructor
    for mem_el in class_members:
        if (mem_el.tag == 'Constructor'):

            # # Check that return type is the class type itself
            # return_type, return_type_kw, return_type_id = utils.findType(mem_el)
            # if return_type_id == class_el.get('id'):

            # Check that the only argument is another class instance
            args = funcutils.getArgs(mem_el)
            if (len(args) == 1) and (args[0]['id'] == class_el.get('id')):

                found_copy_constructor = True
                copy_constr_id = mem_el.get('id')

    if return_id:
        return found_copy_constructor, copy_constr_id
    else:
        return found_copy_constructor

# ====== END: checkCopyConstructor ========



# ====== toWrapperType ========

def toWrapperType(input_type_name, remove_reference=False, remove_pointers=False, use_base_type=False):

    # input_type_name  = NameSpace::SomeType<int>**&
    # output_type_name = SomeType_GAMBIT<int>**&

    type_name = input_type_name

    # Search for '*' and '&'
    n_pointers = type_name.count('*')
    is_ref     = bool('&' in type_name)

    # Remove '*' and '&'
    type_name = type_name.replace('*','').replace('&','')

    # Remove namespace
    type_name = type_name.split('::')[-1]

    if use_base_type:
        abstr_type_name = toAbstractType(input_type_name, include_namespace=True).replace('*','').replace('&','')
        type_name = 'WrapperBase< ' + abstr_type_name + ' >'

    else:
        # Insert wrapper class suffix
        if '<' in type_name:
            type_name_part_one, type_name_part_two = type_name.split('<',1)
            type_name = type_name_part_one + cfg.code_suffix + '<' + type_name_part_two
        else:
            type_name = type_name + cfg.code_suffix

    # Add '*' and '&'
    if remove_pointers:
        pass
    else:
        type_name = type_name + '*'*n_pointers

    if remove_reference:
        pass
    else:
        type_name = type_name + '&'*is_ref

    # Return result
    return type_name

# ====== END: toWrapperType ========



# ====== toAbstractType ========

def toAbstractType(input_type_name, include_namespace=True, add_pointer=False, remove_reference=False, remove_pointers=False):

    # input_type_name  = NameSpace::SomeType<int>**&
    # output_type_name = NameSpace::Abstract__SomeType_GAMBIT<int>**&

    # FIXME:
    # Should this function also translate template argument types?
    # Example: TypeA<TypeB>  -->  Abstract__TypeA<Abstract__TypeB>

    type_name = input_type_name

    # Remove template bracket
    type_name_notempl = utils.removeTemplateBracket(type_name)

    # Search for '*' and '&'
    n_pointers = type_name_notempl.count('*')
    is_ref     = bool('&' in type_name_notempl)

    # Get namespace
    if '::' in type_name_notempl:
        namespace = type_name_notempl.rsplit('::',1)[0]
    else:
        namespace = ''

    # Add abstract class prefix
    # if namespace == '':
    #     type_name = cfg.abstr_class_prefix + type_name
    # else:
    #     type_name_short = utils.removeNamespace(type_name)
    #     type_name = (namespace+'::')*include_namespace  + cfg.abstr_class_prefix + type_name_short


    type_name_short = utils.removeNamespace(type_name)

    if is_ref and remove_reference:
        type_name_short = type_name_short.replace('&','')

    if (n_pointers > 0) and remove_pointers:
        type_name_short = type_name_short.replace('*','')

    if namespace == '':
        type_name = cfg.abstr_class_prefix + type_name_short
    else:
        type_name = (namespace+'::')*include_namespace  + cfg.abstr_class_prefix + type_name_short


    if add_pointer:
        if is_ref and not remove_reference:
            type_name = type_name.rstrip('&') + '*&'
        else:
            type_name = type_name + '*'

    # if remove_reference and is_ref:
    #     type_name = type_name.replace('&','')

    # Return result
    return type_name

# ====== END: toAbstractType ========



# ====== getClassNameDict ========

def getClassNameDict(class_el, abstract=False):

    class_name = {}

    if 'demangled' in class_el.keys():
        class_name['long_templ']  = class_el.get('demangled')
    elif 'name' in class_el.keys():
        class_name['long_templ']  = class_el.get('name')
    else:
        raise Exception('Cannot determine the name of XML element %s' % class_el.get('id'))
    class_name['long']        = class_name['long_templ'].split('<',1)[0]
    class_name['short_templ'] = class_el.get('name')
    class_name['short']       = class_name['short_templ'].split('<',1)[0]

    if abstract:
        abstr_class_name = {}
        abstr_class_name['long_templ']  = getAbstractClassName(class_name['long_templ'], prefix=cfg.abstr_class_prefix)
        abstr_class_name['long']        = abstr_class_name['long_templ'].split('<',1)[0]
        abstr_class_name['short_templ'] = getAbstractClassName(class_name['long_templ'], prefix=cfg.abstr_class_prefix, short=True)
        abstr_class_name['short']       = abstr_class_name['short_templ'].split('<',1)[0]

        return abstr_class_name

    else:
        return class_name

# ====== END: getClassNameDict ========



# # ====== generateWrapperBaseHeader ========

# def generateWrapperBaseHeader(wrapper_base_name):

#     insert_pos = 0

#     # std::shared_ptr lives in the <memory> header
#     code = "#include <memory>\n"
#     code += "\n"

#     # A deleter functor needed by std::shared_ptr
#     code += ("template <typename T>\n"
#              "class Deleter\n"
#              "{\n"
#              "    protected:\n"
#              "        bool member_variable;\n"
#              "\n"
#              "    public:\n"
#              "        Deleter() : member_variable(false) {}\n"
#              "        Deleter(bool memvar_in) : member_variable(memvar_in) {}\n"
#              "\n"
#              "        void setMemberVariable(bool memvar_in)\n"
#              "        {\n"
#              "            member_variable = memvar_in;\n"
#              "        }\n"
#              "\n"
#              "        void operator()(T* ptr)\n"
#              "        {\n"
#              "            if (member_variable==false) { delete ptr; }\n"
#              "        }\n"
#              "};\n")

#     code += "\n"

#     # Wrapper base class code
#     code += ("template <typename T>\n" 
#              "class " + wrapper_base_name + "\n"
#              "{\n"
#              "    public:\n"
#              "        std::shared_ptr<T> BEptr;\n"
#              "\n"
#              "        // Constructor\n"
#              "        WrapperBase(T* BEptr_in, bool memvar_in)\n"
#              "        {\n"
#              "            BEptr.reset(BEptr_in, Deleter<T>(memvar_in));\n"
#              "        }\n"
#              "\n"
#              "        // Special member function to set member_variable in Deleter: \n"
#              "        void _setMemberVariable(bool memvar_in)\n" 
#              "        {\n" 
#              "            std::get_deleter<Deleter<T> >(BEptr)->setMemberVariable(memvar_in);\n"
#              "        }\n"
#              "};\n")
    

#     # Register code and return dict
#     # return_code_dict[fname]['code_tuples'].append( (0, code) )

#     return (insert_pos, code)

# # ====== END: generateWrapperBaseHeader ========



# # ====== getIncludeStatements ========

# def getIncludeStatements(el, convert_loaded_to='none', add_extra_include_path=False, exclude_types=[], input_element='class'):

#     include_statements = []

#     # Check string arguments
#     convert_loaded_to = convert_loaded_to.lower()
#     input_element     = input_element.lower()
#     if convert_loaded_to not in ['none', 'abstract', 'wrapper']:
#         raise Exception("getIncludeStatements: Keyword argument 'convert_loaded_to=' must be either 'none', 'abstract' or 'wrapper'.")
#     if input_element not in ['class', 'function']:
#         raise Exception("getIncludeStatements: Keyword argument 'input_element=' mmust be either 'class' or 'function'.")

#     # Get list of all types used in this class (each entry is a dict)
#     if input_element == 'class':
#         all_types = utils.getAllTypesInClass(el, include_parents=False)
#     elif input_element == 'function':
#         all_types = utils.getAllTypesInFunction(el)

#     # Get file name and line number of the current class
#     class_line_number = int( el.get('line') )
#     class_file_el     = cfg.id_dict[ el.get('file') ]
#     class_file_path   = class_file_el.get('name')

#     # Read file from beginning to position of class definition
#     class_file         = open(class_file_path, 'r')
#     class_file_content = class_file.readlines()[0:class_line_number]
#     class_file_content = ''.join(class_file_content)
#     class_file.close()

#     # Identify included header files from this file (utils.identifyIncludedHeaders returns a dict of the form {header_file_name: xml_id})
#     included_headers_dict = utils.identifyIncludedHeaders(class_file_content, only_native=True)

#     # Move up the header tree and identify all the relevant (native) included headers
#     header_paths = [ cfg.id_dict[file_id].get('name') for file_id in included_headers_dict.values() ]
#     header_paths_done = []

#     while len(header_paths) > 0:

#         header_path = header_paths.pop()

#         # Read header
#         header         = open(header_path, 'r')
#         header_content = header.read()
#         header.close()

#         # Identify new headers
#         new_included_headers = utils.identifyIncludedHeaders(header_content, only_native=True)

#         # Add any new headers to included_headers_dict
#         for file_name, file_id in new_included_headers.items():
#             if file_name not in included_headers_dict.keys():
#                 included_headers_dict[file_name] = file_id

#         # Add any new headers to the list of header files to check
#         new_header_paths = [ cfg.id_dict[header_id].get('name') for header_id in new_included_headers.values() ]
#         for new_path in new_header_paths:
#             if (new_path not in header_paths) and (new_path not in header_paths_done):
#                 header_paths.append(new_path)

#         # Keep track of headers we've done
#         header_paths_done.append(header_path)


#     # Determine what include statements to generate:

#     for type_dict in all_types:

#         type_el   = type_dict['el']
#         type_name = type_dict['class_name']

#         if type_name in exclude_types:
#             continue

#         if utils.isAcceptedType(type_el):

#             if utils.isFundamental(type_el):
#                 pass

#             elif utils.isLoadedClass(type_el):

#                 # For each loaded class used in this class, check whether the corresponding class definition can be
#                 # found in the current file (above current class) or among the included headers. Only include headers 
#                 # for the loaded classes that pass this check. (If no such class definition is found, it must be a 
#                 # case of simply using forward declaration, in which case we should *not* include the corresponding header.) 

#                 type_file_id = type_el.get('file')
#                 type_line_number = int(type_el.get('line'))

#                 if (type_file_id in included_headers_dict.values()) :
#                     type_definition_found = True
#                 elif (type_file_id == el.get('file')) and (type_line_number < class_line_number):
#                     type_definition_found = True
#                 else:
#                     type_definition_found = False

#                 if type_definition_found:
#                     if convert_loaded_to == 'none':

#                         type_file_el = cfg.id_dict[type_file_id]
#                         type_file_full_path = type_file_el.get('name')

#                         if utils.isHeader(type_file_el):
#                             type_file_basename = os.path.basename(type_file_full_path)
#                             if add_extra_include_path:
#                                 include_statements.append('#include "' + os.path.join(cfg.add_path_to_includes, type_file_basename) + '"')
#                             else:
#                                 include_statements.append('#include "' + type_file_basename + '"')
#                         else:
#                             print 'INFO: ' + 'Found declaration of loaded type "%s" in file "%s", but this file is not recognized as a header file. No header file include statement generated.' % (type_name['long_templ'], type_file_full_path)

#                     else:
#                         if add_extra_include_path:
#                             include_statements.append('#include "' + os.path.join(cfg.add_path_to_includes, cfg.new_header_files[type_name['long']][convert_loaded_to]) + '"')
#                         else:
#                             include_statements.append('#include "' + cfg.new_header_files[type_name['long']][convert_loaded_to] + '"')

#                 else:
#                     # This must be a case of a type that is only forward declared. Don't include any header (as this will typically lead to a 'header loop').
#                     continue

#             elif utils.isStdType(type_el):

#                 if type_name['long'] in cfg.known_class_headers:
#                     header_name = cfg.known_class_headers[type_name['long']]
#                     if (header_name[0] == '<') and (header_name[-1] == '>'):
#                         include_statements.append('#include ' + cfg.known_class_headers[type_name['long']])
#                     else:
#                         include_statements.append('#include "' + cfg.known_class_headers[type_name['long']] + '"')
#                 else:
#                     # warnings.warn("The standard type '%s' has no specified header file. Please update modules/cfg.py. No header file included." % type_name['long_templ'])
#                     print 'INFO: ' + 'The standard type "%s" has no specified header file. Please update modules/cfg.py. No header file include statement generated.' % type_name['long_templ']

#             else:
#                 if type_name['long'] in cfg.known_class_headers:
#                     include_statements.append('#include "' + cfg.known_class_headers[type_name['long']] + '"')
#                 else:
#                     # warnings.warn("The type '%s' has no specified header file. Please update modules/cfg.py. No header file included." % type_name['long'])
#                     print 'INFO: ' + 'The type "%s" has no specified header file. Please update modules/cfg.py. No header file include statement generated.' % type_name['long']
#         else:
#             # warnings.warn("The type '%s' is unknown. No header file included." % type_name['long'])
#             print 'INFO: ' + 'The type "%s" is unknown. No header file include statement generated.' % type_name['long']

#     # Remove duplicates and return list
#     include_statements = list( set(include_statements) )

#     return include_statements

# # ====== END: getIncludeStatements ========
