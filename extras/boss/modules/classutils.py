####################################
#                                  #
#  Utility functions for handling  # 
#  C++ classes with BOSS           #
#                                  #
####################################

import xml.etree.ElementTree as ET
from collections import OrderedDict
import warnings

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

def constrEmptyTemplClassDecl(short_abstract_class_name, namespaces, template_bracket, indent=4):

    n_indents  = len(namespaces)
    class_decl = ''

    # - Construct the beginning of the namespaces
    class_decl += utils.constrNamespace(namespaces, 'open')
    # for ns in namespaces:
    #     class_decl += ' '*n_indents*indent + 'namespace ' + ns + '\n'
    #     class_decl += ' '*n_indents*indent + '{' + '\n'
    #     n_indents += 1

    class_decl += ' '*n_indents*indent + 'template ' + template_bracket + '\n'
    class_decl += ' '*n_indents*indent + 'class ' + short_abstract_class_name + ' {};\n'

    # - Construct the closing of the namespaces
    class_decl += utils.constrNamespace(namespaces, 'close')
    # for ns in namespaces:
    #     n_indents -= 1
    #     class_decl += ' '*n_indents*indent + '}' + '\n'
    
    class_decl += '\n'

    return class_decl

# ====== END: constrEmptyTemplClassDecl ========



# ====== constrTemplForwDecl ========

def constrTemplForwDecl(short_class_name, namespaces, template_bracket, indent=4):

    n_indents = len(namespaces)
    forw_decl = ''

    # - Construct the beginning of the namespaces
    forw_decl += utils.constrNamespace(namespaces, 'open')
    # for ns in namespaces:
    #     forw_decl += ' '*n_indents*indent + 'namespace ' + ns + '\n'
    #     forw_decl += ' '*n_indents*indent + '{' + '\n'
    #     n_indents += 1

    forw_decl += ' '*n_indents*indent + 'template ' + template_bracket + '\n'
    forw_decl += ' '*n_indents*indent + 'class ' + short_class_name + ';\n'

    # - Construct the closing of the namespaces
    forw_decl += utils.constrNamespace(namespaces, 'close')
    # for ns in namespaces:
    #     n_indents -= 1
    #     forw_decl += ' '*n_indents*indent + '}' + '\n'
    
    forw_decl += '\n'

    return forw_decl

# ====== END: constrTemplForwDecl ========



# ====== constrAbstractClassDecl ========

def constrAbstractClassDecl(class_el, short_class_name, short_abstract_class_name, namespaces, indent=4, template_types=[], 
                            has_copy_constructor=True,  has_assignment_operator=True):

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

    #
    # Construct the abstract class declaration
    #
    
    class_decl = '#pragma GCC diagnostic push\n'
    class_decl += '#pragma GCC diagnostic ignored "-Wunused-parameter"\n'
    class_decl += '#pragma GCC diagnostic ignored "-Wreturn-type"\n'


    # - Construct the name of the abstract class, including full namespace
    # short_class_name    = class_el.get('name')
    # full_class_name     = class_el.get('demangled').replace(' ','')
    # abstract_class_name = getAbstractClassName(full_class_name, prefix=cfg.abstr_class_prefix)
    # short_abstract_class_name = abstract_class_name.rsplit('::',1)[-1]
    # namespaces = abstract_class_name.split('::')[:-1]

    # - Construct the beginning of the namespaces
    class_decl += utils.constrNamespace(namespaces, 'open')
    # for ns in namespaces:
    #     class_decl += ' '*n_indents*indent + 'namespace ' + ns + '\n'
    #     class_decl += ' '*n_indents*indent + '{' + '\n'
    #     n_indents += 1

    # # - Forward declare the child class (so that it can be used for input arguments)
    # #   (for template classes this is already taken care of)
    # if not is_template:
    #     class_decl += ' '*n_indents*indent + 'class ' + short_class_name + ';\n' #.split('<',1)[0] + ';\n'

    # - If this class is a template specialization, add 'template <>' at the top
    if is_template == True:
        class_decl += ' '*n_indents*indent + 'template <>\n'

    # - Construct the declaration line, with inheritance of abstract classes
    base_el_list = class_el.findall('Base')
    if len(base_el_list) > 0:
        inheritance_line = ''
        for base_el in base_el_list:
            base_id            = base_el.get('type')
            base_name          = cfg.id_dict[base_id].get('demangled')

            if base_name in cfg.loaded_classes:
                abstract_base_name = getAbstractClassName(base_name, prefix=cfg.abstr_class_prefix)
                base_access        = base_el.get('access')
                inheritance_line  += base_access + ' ' + 'virtual' + ' ' + abstract_base_name + ', '
            else: 
                pass
        # - Remove trailing whitespace and comma, add colon
        if inheritance_line != '':
            inheritance_line = ' : ' + inheritance_line.rstrip(', ')

    else:
        inheritance_line = ''
    class_decl += ' '*n_indents*indent
    
    if is_template:
        class_decl += 'class ' + short_abstract_class_name + '<' + ','.join(template_types) + '>' + inheritance_line + '\n'
    else:
        class_decl += 'class ' + short_abstract_class_name + inheritance_line + '\n'

    # - Construct body of class declaration
    current_access = ''
    class_decl += ' '*n_indents*indent
    class_decl += '{' + '\n'
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

        elif el.tag == 'Method':

            # Check if this member function should be ignored based on unloaded types.
            if funcutils.ignoreFunction(el):
                print 'INFO: Class ' + class_el.get('demangled') + ' : Member function ' + el.get('name') + ' ignored due to non-accepted type(s).'
                continue

            return_type, return_kw, return_id = utils.findType( cfg.id_dict[el.get('returns')] )
            return_type_base = return_type.replace('*','').replace('&','')
            
            return_kw_str = ' '.join(return_kw)        
            return_kw_str += ' '*bool(len(return_kw))

            return_el = cfg.id_dict[return_id]
            return_is_native = utils.isNative(return_el)
            args = funcutils.getArgs(el)

            if (return_type == 'void') or (return_type.count('*') > 0):
                w_return_type = return_type
            else:
                w_return_type = return_type.replace(return_type_base, return_type_base+'*')

            w_args = funcutils.constrWrapperArgs(args)
            w_args_bracket = funcutils.constrArgsBracket(w_args)
            w_args_bracket_notypes = funcutils.constrArgsBracket(w_args, include_arg_type=False)
            w_func_name = el.get('name') + cfg.code_suffix
            

            #
            # If the method makes use of a native type, construct a pair of wrapper methods.
            #
            if funcutils.usesNativeType(el):

                # Construct the virtual member function that is overridden, e.g.:  
                #
                #   virtual X* getX_GAMBIT(arguments) {}
                #
                class_decl += '\n'
                class_decl += ' '*(n_indents+2)*indent
                class_decl += 'virtual ' + return_kw_str + w_return_type + ' ' + w_func_name + w_args_bracket + ' {std::cout << "Called virtual function" << std::endl;};' + '\n'


                # Construct the member function with the original name, wrapping the overridden one, e.g.: 
                #
                #   Abstract__X* getX(arguments)
                #   {
                #     return reinterpret_cast<Abstract__X*>(getX_GAMBIT(arguments));
                #   }
                #
                if return_is_native:
                    # w2_return_type = cfg.abstr_class_prefix + w_return_type
                    w2_return_type = getAbstractClassName(w_return_type)
                else:
                    w2_return_type = w_return_type

                w2_func_name = el.get('name')
                w2_args_bracket = w_args_bracket

                class_decl += ' '*(n_indents+2)*indent 
                class_decl += return_kw_str + w2_return_type + ' ' + w2_func_name + w2_args_bracket +'\n'
                class_decl += ' '*(n_indents+2)*indent + '{\n'
                if (return_type == 'void'):
                    class_decl += ' '*(n_indents+3)*indent + w_func_name + w_args_bracket_notypes + ';\n'
                else:
                    if return_is_native:
                        class_decl += ' '*(n_indents+3)*indent + 'return reinterpret_cast<' + return_kw_str + w2_return_type + '>(' + w_func_name + w_args_bracket_notypes + ');\n'
                    else:
                        class_decl += ' '*(n_indents+3)*indent + 'return ' + w_func_name + w_args_bracket_notypes + ';\n'
                class_decl += ' '*(n_indents+2)*indent + '}\n'

            #
            # If the method makes only makes use of fundamental/standard types, construct a single virtual method
            #
            else:
                class_decl += '\n'
                class_decl += ' '*(n_indents+2)*indent
                class_decl += 'virtual ' + return_kw_str + return_type + ' ' + el.get('name') + w_args_bracket + ' {std::cout << "Called virtual function" << std::endl;};' + '\n'


        #
        # If element is a public member variable, construct virtual method that returns a reference to this variable
        #
        elif (el.tag in ('Field', 'Variable')) and (el.get('access') == 'public'):

            class_decl += '\n' 
            class_decl += constrVariableRefFunction(el, virtual=True, indent=indent, n_indents=n_indents+2)

        else:
            class_decl += ' '*(n_indents+2)*indent
            class_decl += '// UNKNOWN: ' + el.tag + '\n'

    
    # - Construct 'pointerAssign' and 'pointerCopy' functions
    if has_copy_constructor or has_assignment_operator:
        class_decl += '\n'
        class_decl += ' '*(n_indents+1)*indent + 'public:\n'
        if has_assignment_operator:
            class_decl += constrPtrAssignFunc(short_abstract_class_name, short_class_name, virtual=True, indent=indent, n_indents=n_indents+2)
        if has_copy_constructor:
            class_decl += constrPtrCopyFunc(short_abstract_class_name, short_class_name, virtual=True, indent=indent, n_indents=n_indents+2)

    # # - Construct the 'downcast' converter function
    # class_decl += '\n'
    # # class_decl += ' '*(n_indents+1)*indent + 'NEW CODE HERE!\n'
    # if is_template:
    #     template_bracket = '<' + ','.join(template_types) + '>'
    # else:
    #     template_bracket = ''
    # class_decl += ' '*(n_indents+1)*indent + 'public:\n'
    # class_decl += constrDowncastFunction(short_class_name, indent=indent, n_indents=n_indents+2, template_bracket=template_bracket)

    # - Construct an empty virtual destructor
    class_decl += ' '*(n_indents+2)*indent
    class_decl += 'virtual ~' + short_abstract_class_name + '() {};\n'

    # - Close the class body
    class_decl += ' '*n_indents*indent + '};' + '\n'

    # - Construct the closing of the namespaces
    class_decl += utils.constrNamespace(namespaces, 'close')
    # for ns in namespaces:
    #     n_indents -= 1
    #     class_decl += ' '*n_indents*indent + '}' + '\n'


    class_decl += '#pragma GCC diagnostic pop\n'
    return class_decl

# ====== END: constrAbstractClassDecl ========



# ====== constrFactoryFunction ========

def constrFactoryFunction(class_el, full_class_name, indent=4, template_types=[]):

    # Create list of all constructors of the class
    constructor_elements = []

    # Replace '*' and '&' in list of template types
    template_types = [e.replace('*','P').replace('&','R') for e in template_types]

    if 'members' in class_el.keys():
        for mem_id in class_el.get('members').split():
            el = cfg.id_dict[mem_id]
            if (el.tag == 'Constructor') and (el.get('access') == 'public'): #and ('artificial' not in el.keys()):  #(el.get('explicit') == "1"):
                constructor_elements.append(el)

    # Construct factory function definition(s)
    func_def = ''
    for el in constructor_elements:

        # Useful variables
        # full_class_name  = el.get('demangled').rsplit('::',1)[0]
        short_class_name = full_class_name.split('::')[-1].split('<')[0]
        factory_name = 'Factory_' + short_class_name
        if len(template_types) > 0:
            factory_name += '_' + '_'.join(template_types)

        # Identify arguments
        args = funcutils.getArgs(el)

        # Translate argument type of loaded classes
        w_args = funcutils.constrWrapperArgs(args, add_ref=True)

        # Invent argument names if missing
        argc = 1
        for i in range(len(args)):
            if args[i]['name'] == '':
                args[i]['name'] = 'arg_' + str(argc)
                argc += 1

        # Construct bracket with input arguments
        args_bracket = funcutils.constrArgsBracket(w_args, include_namespace=False)
        args_bracket_notypes = funcutils.constrArgsBracket(args, include_arg_type=False, cast_to_original=True)
        
        # Generate declaration line:
        return_type = getAbstractClassName(full_class_name, prefix=cfg.abstr_class_prefix)
        func_def += return_type + '* ' + factory_name + args_bracket + '\n'
        # Generate body
        func_def += '{' + '\n'
        func_def += indent*' ' + 'return new ' + full_class_name + args_bracket_notypes + ';' + '\n'
        func_def += '}' + 2*'\n'

    return func_def

# ====== END: constrFactoryFunction ========



# ====== constrDowncastFunction ========

def constrDowncastFunction(cast_to_class, indent=4, n_indents=0, template_bracket=''):

    # Example of output C++ code:
    #
    #     public: A* downcast() 
    #     {
    #         return reinterpret_cast<A*>(this); 
    #     }
    #

    base_indent = n_indents*indent
    use_class_name = cast_to_class + template_bracket

    func_def  = ''
    func_def += base_indent*' ' + use_class_name + '* downcast()\n'
    func_def += base_indent*' ' + '{\n'
    func_def += (base_indent+indent)*' ' + 'return reinterpret_cast<' + use_class_name +'*>(this);\n'
    func_def += base_indent*' ' + '}\n'

    return func_def

# ====== END: constrDowncastFunction ========



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

def constrWrapperFunction(method_el, indent=cfg.indent, n_indents=0):

    # Function name
    func_name = method_el.get('name')

    # Function return type
    return_type, return_kw, return_id = utils.findType( cfg.id_dict[method_el.get('returns')] )
    return_el = cfg.id_dict[return_id]

    return_is_loaded_class = utils.isLoadedClass(return_el)
    pointerness, is_ref = utils.pointerAndRefCheck(return_type, byname=True)

    # Function arguments (get list of dicts with argument info)
    args = funcutils.getArgs(method_el)

    print '------------'
    print 'FUNC NAME: ', func_name
    print 'TYPE: ', return_type
    print 'POINTER:', pointerness
    print 'KW  : ', return_kw
    print 'ID  : ', return_id
    print 'LOADED:', return_is_loaded_class

    # Construct wrapper function name
    w_func_name = funcutils.constrWrapperName(method_el)

    # Choose wrapper return type
    if return_is_loaded_class:
        if pointerness == 0:
            w_return_type = toAbstractType(return_type, include_namespace=True, add_pointer=True)
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
    w_args_bracket = funcutils.constrArgsBracket(w_args)

    # Construct declaration line for wrapper function
    w_func_line = funcutils.constrDeclLine(w_return_type, w_func_name, w_args_bracket, keywords=return_kw)

    # Construct function body for wrapper function
    w_func_body = funcutils.constrWrapperBody(return_type, func_name, args, return_is_loaded_class, keywords=return_kw)

    # Combine code and add indentation
    wrapper_code  = ''
    wrapper_code += utils.addIndentation(w_func_line, n_indents*indent) + '\n'
    wrapper_code += utils.addIndentation(w_func_body, n_indents*indent) + '\n'

    # Return result
    return wrapper_code

# ====== END: constrWrapperFunction ========



# ====== constrVariableRefFunction ========

def constrVariableRefFunction(var_el, virtual=False, indent=cfg.indent, n_indents=0):

    func_code = ''

    var_name = var_el.get('name')
    ref_method_name = var_name + '_ref' + cfg.code_suffix

    var_type, var_kw, var_id = utils.findType( cfg.id_dict[var_el.get('type')] )
    var_type_base = var_type.replace('*','').replace('&','')
    
    var_kw_str = ' '.join(var_kw)        
    var_kw_str += ' '*bool(len(var_kw))

    var_el = cfg.id_dict[var_id]
    var_is_native = utils.isNative(var_el)
    if var_is_native:
        return_type = getAbstractClassName(var_type_base, prefix=cfg.abstr_class_prefix, short=True)
    else:
        return_type = var_type_base

    func_code += ' '*n_indents*indent
    if virtual:
        func_code += 'virtual ' + var_kw_str + return_type + '& ' + ref_method_name + '() {std::cout << "Called virtual function" << std::endl;};\n'
    else:
        func_code += var_kw_str + return_type + '& ' + ref_method_name + '() { return ' + var_name  +'; }\n'

    return func_code

# ====== END: constrVariableRefFunction ========



# ====== constrPtrCopyFunc ========

def constrPtrCopyFunc(short_abstr_class_name, short_class_name, virtual=False, indent=cfg.indent, n_indents=0):

    ptr_code = ''
    ptr_code += ' '*cfg.indent*n_indents
    
    if virtual:
        ptr_code += 'virtual '+ short_abstr_class_name + '*' + ' pointerCopy' + cfg.code_suffix + '() {std::cout << "Called virtual function" << std::endl;};\n'   
    else:
        ptr_code += short_abstr_class_name + '*' + ' pointerCopy' + cfg.code_suffix + '()'
        ptr_code += ' ' + '{ return new ' + short_class_name + '(*this); }\n'

    return ptr_code

# ====== END: constrPtrCopyFunc ========



# ====== constrPtrAssignFunc ========

def constrPtrAssignFunc(short_abstr_class_name, short_class_name, virtual=False, indent=cfg.indent, n_indents=0):

    ptr_code = ''
    ptr_code += ' '*cfg.indent*n_indents
    
    if virtual:
        ptr_code += 'virtual void pointerAssign' + cfg.code_suffix + '(' + short_abstr_class_name + '* in) {std::cout << "Called virtual function" << std::endl;};\n'
    else:
        ptr_code += 'void pointerAssign' + cfg.code_suffix + '(' + short_abstr_class_name + '* in)'
        ptr_code += ' ' + '{ *this = *dynamic_cast<' + short_class_name + '*>(in); }\n'        

    return ptr_code

# ====== END: constrPtrAssignFunc ========



# ====== checkAssignmentOperator ========

def checkAssignmentOperator(class_el):

    found_assignment_operator = False

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
                    print 'FOUND ASSIGNMENT OPERATOR:', mem_el.get('name'), return_type, args
    
    return found_assignment_operator

# ====== END: checkAssignmentOperator ========



# ====== checkCopyConstructor ========

def checkCopyConstructor(class_el):

    found_copy_constructor = False

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
                print 'FOUND COPY CONSTRUCTOR:', mem_el.get('name'), args
    
    return found_copy_constructor

# ====== END: checkCopyConstructor ========



# ====== toWrapperType ========

def toWrapperType(input_type_name):

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

    # Insert wrapper class suffix
    if '<' in type_name:
        type_name_part_one, type_name_part_two = type_name.split('<',1)[0]
        type_name = type_name_part_one + cfg.code_suffix + '<' + type_name_part_two
    else:
        type_name = type_name + cfg.code_suffix

    # Add '*' and '&'
    type_name = type_name + '*'*n_pointers + '&'*is_ref

    # Return result
    return type_name

# ====== END: toWrapperType ========



# ====== toAbstractType ========

def toAbstractType(input_type_name, include_namespace=True, add_pointer=False):

    # input_type_name  = NameSpace::SomeType<int>**&
    # output_type_name = NameSpace::Abstract__SomeType_GAMBIT<int>**&

    # FIXME:
    # Should this function also translate template argument types?
    # Example: TypeA<TypeB>  -->  Abstract__TypeA<Abstract__TypeB>

    type_name = input_type_name

    # Search for '*' and '&'
    n_pointers = type_name.count('*')
    is_ref     = bool('&' in type_name)

    # Get namespace
    type_name_notemplate = type_name.split('<',1)[0]
    if '::' in type_name_notemplate:
        namespace = type_name.split('<',1)[0].rsplit('::',1)[0]
    else:
        namespace = ''

    # Add abstract class prefix
    if namespace == '':
        type_name = cfg.abstr_class_prefix + type_name
    else:
        type_name_short = type_name.lstrip(namespace)
        type_name = (namespace+'::')*include_namespace  + cfg.abstr_class_prefix + type_name_short

    if add_pointer:
        if is_ref:
            type_name = type_name.rstrip('&') + '*&'
        else:
            type_name = type_name + '*'

    # Return result
    return type_name

# ====== END: toAbstractType ========