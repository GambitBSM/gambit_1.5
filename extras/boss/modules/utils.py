################################
#                              #
#  Utility functions for BOSS  #
#                              #
################################

import xml.etree.ElementTree as ET
from collections import OrderedDict
from operator import itemgetter
import os
import warnings

import modules.cfg as cfg


# ====== isFundamental ========

def isFundamental(el):

    is_fundamental = False

    if el.tag == 'FundamentalType':
        is_fundamental = True

    return is_fundamental

# ====== END: isFundamental ========


# ====== isNative ========

def isNative(el):

    # Makes use of global variables:  accepted_paths

    is_native = False
    can_check_tags = ['Class', 'Constructor', 'Converter', 'Destructor', 'Enumeration', 
                      'Field', 'File', 'Function', 'Method', 'OperatorFunction', 
                      'OperatorMethod', 'Struct', 'Typedef', 'Union', 'Variable']

    if el.tag == 'FundamentalType':
        is_native = False

    elif el.tag in can_check_tags:

        if el.tag == 'File':
            file_el = el
        else:
            file_el = cfg.id_dict[el.get('file')]

        check_path = file_el.get('name')

        is_native = False
        for accepted_path in cfg.accepted_paths:
            if accepted_path in os.path.dirname(check_path):
                is_native = True
                break

        # if os.path.dirname(check_path) in cfg.accepted_paths:
        #     is_native = True
        # else:
        #     is_native = False

    else:
        raise Exception('Cannot check whether XML element with id="%s" and tag "%s" is native.' % (el.get('id'), el.tag))

    return is_native

# ====== END: isNative ========



# ====== isStdType ========

def isStdType(el):

    # Makes use of global variables:  accepted_paths

    # print 'checking:', el.tag
    # if 'name' in el.keys():
    #     print ' -- ', el.get('name')
    # print

    is_std = False
    can_check_tags = ['Class', 'Struct', 'Union', 'Enumeration']

    if el.tag in can_check_tags:

        if 'demangled' in el.keys():
            demangled_name = el.get('demangled')
            if demangled_name[0:5] == 'std::':
                is_std = True

        # file_el = cfg.id_dict[el.get('file')]

        # check_path = file_el.get('name')

        # for std_path in cfg.std_include_paths:
        #     if std_path in os.path.dirname(check_path):
        #         is_std = True
        #         if check_path not in cfg.std_headers_used:
        #             cfg.std_headers_used.append(check_path)
        #         break

    else:
        pass
        # warnings.warn('Cannot check whether XML element with id="%s" and tag "%s" is a standard type.' % (el.get('id'), el.tag))
        # raise Exception('Cannot check whether XML element with id="%s" and tag "%s" is a standard type.' % (el.get('id'), el.tag))

    return is_std

# ====== END: isStdType ========



# ====== getTemplateBracket ========

def getTemplateBracket(el):

    src_file_name = cfg.id_dict[el.get('file')].get('name')
    line_number   = int(el.get('line'))

    print
    print 'ID:       ', el.get('id')
    print 'SRC_FILE: ', src_file_name
    print 'LINE:     ', line_number
    print

    f = open(src_file_name, 'r')
    file_content = f.read()
    f.close()
    file_content_nocomments = removeComments(file_content, insert_blanks=True)

    # Find index of the \n in line number line_number
    count = 0
    prev_pos = 0
    for index,char in enumerate(file_content_nocomments):
        if char=='\n':
            count += 1
        if count == line_number:
            break
        if char=='\n':         # STUPID HACK
            prev_pos = index

    newline_pos = index


    # Find the template parameter bracket, e.g. <typename A, typename B>
    search_content = file_content_nocomments[:newline_pos]

    start_pos = 0
    end_pos = search_content.rfind('>')
    if end_pos != -1:
        balance = -1
        for i in range(end_pos-1, -1, -1):
            char = search_content[i]
            if char == '>':
                balance -= 1
            elif char == '<':
                balance += 1
            if (balance == 0):
                start_pos = i
                break
        template_bracket = search_content[start_pos:end_pos+1]
    else:
        template_bracket = '<>'

    print 'TEMPLATE BRACKET: ', template_bracket
    # # Isolate template variable names (remove keywords 'class' and 'typename')
    # temp_var_list = template_bracket[1:-1].split(',')
    # temp_var_list = [ e.strip() for e in temp_var_list]
    # for i,arg in enumerate(temp_var_list):
    #     arg_split = arg.split()
    #     if 'class' in arg_split:
    #         arg_split.remove('class')
    #     if 'typename' in arg_split:
    #         arg_split.remove('typename')
    #     temp_var_list[i] = ' '.join(arg_split)

    # Isolate only the template variable names (last word in each entry)
    if template_bracket == '<>':
        temp_var_list = []
    else:
        temp_var_list = template_bracket[1:-1].split(',')
        temp_var_list = [ e.strip() for e in temp_var_list]
        temp_var_list = [ e.split()[-1] for e in temp_var_list]

    # Return result
    return template_bracket, temp_var_list

# ====== END: getTemplateBracket ========



# ====== getSpecTemplateTypes ========

def getSpecTemplateTypes(el):

    # Classes and functions must be treated differently
    if el.tag == 'Class':
        spec_types = el.get('name').split('<',1)[1].rsplit('>',1)[0].strip()
    elif el.tag == 'Function':
        spec_types = el.get('demangled').split('<',1)[1].rsplit('>',1)[0].strip()
    else:
        raise Exception("Don't know how to get template types from XML element with tag: %s" % el.tag)

    # Identify the correct commas
    pos = []
    balance = 0
    for i,c in enumerate(spec_types):
        if c == '<':
            balance += 1
        if c == '>':
            balance -= 1

        if (balance==0) and (c == ','):
            pos.append(i)

    # Construct list of arguments
    spec_types_list = []
    prev_p = 0
    for p in pos:
        spec_types_list.append(spec_types[prev_p:p])
        prev_p = p+1
    spec_types_list.append(spec_types[prev_p:])

    # Return result
    return spec_types_list

# ====== END: getSpecTemplateTypes ========



# ====== removeComments ========

def removeComments(content, insert_blanks=False):

    # Prepare list for storing tuples of the form: (start_position, stop_position)
    content_lenght = len(content)
    comment_sections = []

    #
    # Locate comments:
    #

    # -- One-line comments
    temp_startpos = 0
    while True:

        # Find start of comment
        search_pos = content[temp_startpos:].find('//')
        if search_pos == -1:
            break
        else:
            comment_start = temp_startpos + search_pos

            # Find end of comment
            search_pos = content[comment_start:].find('\n')
            if search_pos == -1:
                comment_end = content_lenght - 1
            else:
                comment_end = comment_start + search_pos

            # Store positions
            comment_sections.append( (comment_start, comment_end) )

            # Update loop variable
            temp_startpos = comment_end


    # -- Multi-line comments
    temp_startpos = 0
    while True:

        # Find start of comment
        search_pos = content[temp_startpos:].find('/*')
        if search_pos == -1:
            break
        else:
            comment_start = temp_startpos + search_pos

            # Find end of comment
            search_pos = content[comment_start:].find('*/')
            if search_pos == -1:
                comment_end = content_lenght - 1
            else:
                comment_end = comment_start + search_pos + 1

            # Store positions
            comment_sections.append( (comment_start, comment_end) )

            # Update loop variable
            temp_startpos = comment_end


    # Sort comment_sections from last to first, depending on stop position
    comment_sections = sorted(comment_sections, key=itemgetter(1), reverse=True)


    # Remove comments
    prev_start_pos = 0
    prev_stop_pos  = 0

    for start_pos, stop_pos in comment_sections:
        new_lenght = len(content)

        # Skip if the current comment was contained within the previous removed comment
        if (start_pos > prev_start_pos) and (stop_pos < prev_stop_pos):
            continue
        # If not, go on to remove comment
        else:
            # Insert whitespace?
            if insert_blanks == False:
                content = content.replace( content[start_pos:stop_pos+1], '')
            else:
                # Construct string of spaces and newlines to replace comments
                insert_string = ''
                for char in content[start_pos:stop_pos+1]:
                    insert_string += ' '*(char!='\n') + '\n'*(char=='\n')

                # Perform replacement
                content = content.replace( content[start_pos:stop_pos+1], insert_string )

            # Update loop variables
            prev_start_pos = start_pos
            prev_stop_pos  = stop_pos

    return content

# ====== END: removeComments ========



# ====== findType ========

def findType(el_input):

    check_keywords = ['static', 'const']
    additional_keywords = []
    refs_and_pointers = ''

    el = el_input
    
    if el.tag in ['FundamentalType', 'Class', 'Struct']:
        type_id = el.get('id')
    else:
        type_id = el.get('type')
        while ('type' in el.keys()) or ('returns' in el.keys()):

            # Get xml id to move further through the xml file
            if el.tag in ['FunctionType', 'Function', 'Method', 'OperatorMethod']:
                type_id = el.get('returns')
            else:
                type_id = el.get('type')

            # Check for reference or pointer type
            if el.tag == 'ReferenceType':
                refs_and_pointers = ''.join(['&', refs_and_pointers])
            if el.tag == 'PointerType':
                refs_and_pointers = ''.join(['*', refs_and_pointers])

            # Pick up any extra keywords
            for keyword in check_keywords:
                if keyword in el.keys():
                    additional_keywords.append(keyword)

            # change xlm element 'el'
            el = cfg.id_dict[type_id]

    # Remove duplicates from the additional_keywords list
    additional_keywords = list(set(additional_keywords))

    # When we exit the loop, 'el' is at the final element.
    # The typename is in the 'name'/'demangled' attribute
    if 'demangled' in el.keys():
        typename = el.get('demangled') + refs_and_pointers    
    else:
        typename = el.get('name') + refs_and_pointers

    return typename, additional_keywords, type_id


# ====== END: findType ========



# ====== findNewLinePos ========

def findNewLinePos(content, line_number):
    count = 0
    for index,char in enumerate(content):
        if char=='\n':
            count += 1
        if count == line_number:
            break
    newline_pos = index

    return newline_pos

# ====== END: findNewLinePos ========



# ====== generateAbstractFilePath ========

# def generateAbstractFilePath(el, short_name, prefix='abstract_', suffix='.hpp'):

#     short_abstr_fname  = prefix + short_name.lower() + suffix
#     src_file_el = cfg.id_dict[el.get('file')]
#     dir_name = os.path.split(src_file_el.get('name'))[0]
    
#     return os.path.join(dir_name,short_abstr_fname)


# ====== END: generateAbstractFilePath ========



# ====== getBracketPositions ========

def getBracketPositions(content, delims=['{','}']):

    # Input:  
    # - Content string
    # - List of left and right delimiters
    #
    # Output: 
    # - List of bracket positions: [l_pos, r_pos]
    #

    l_delim, r_delim = delims

    # Cannot handle identical left and right delimiters
    if l_delim == r_delim:
        raise Exception('Left and right delimiters cannot be identical.')

    # Prepare search
    bracket_count  = 0
    start_count = False
    balance     = False
    l_pos = 0
    r_pos = 0

    # Search
    for i,c in enumerate(content):

        if c == l_delim:
            bracket_count += 1
            if not start_count:
                l_pos = i
                start_count = True

        if start_count and c == r_delim:
            bracket_count -= 1
            r_pos = i 

        if start_count and bracket_count == 0:
            balance = True
            break

    # If brackets did not balance, raise exception
    if not balance:
        raise ValueError("No matching right delimiter for the first left delimiter.")
    # Else, return the found bracket positions
    else:
        return [l_pos, r_pos]

# ====== END: getBracketPositions ========



# ====== addIndentation ========

def addIndentation(content, indent):

    lines = content.split('\n')
    new_content = '\n'.join( [' '*indent + line for line in lines] )

    # If the last char in content was a newline, 
    # remove the indentation that was added after that newline
    if lines[-1] == '':
        new_content = new_content[:-indent]

    return new_content

# ====== END: addIndentation ========



# ====== getNamespaces ========

def getNamespaces(xml_el):

    namespaces = []

    if 'demangled' in xml_el.keys():
        full_name = xml_el.get('demangled')
        modified_name = full_name.split('<')[0]

        namespaces = modified_name.split('::')[:-1]

    return namespaces

# ====== END: getNamespaces ========



# ====== removeTemplateBracket ========

def removeTemplateBracket(type_name):
    
    if '<' in type_name:
        type_name_notempl = type_name.split('<',1)[0] + type_name.rsplit('>',1)[-1]
    else:
        type_name_notempl = type_name

    return type_name_notempl

# ====== END: removeTemplateBracket ========



# ====== removeNamespace ========

def removeNamespace(type_name):
    
    type_name_notempl = removeTemplateBracket(type_name)

    if '::' in type_name_notempl:
        namespace = type_name_notempl.rsplit('::',1)[0]
        new_type_name = type_name.replace(namespace+'::','',1)
    else:
        new_type_name = type_name

    return new_type_name

# ====== END: removeNamespace ========



# ====== isAcceptedType ========

def isAcceptedType(input_type, byname=False):
    #
    # By default, input_type is assumed to be an XML element (xml.etree.ElementTree.Element).
    # If byname=True, input_type is a simple string.
    #

    is_accepted_type = False

    # If input_type is a string
    if byname:
        base_input_type = input_type.replace('*', '').replace('&','')
        if base_input_type in cfg.accepted_types:
            is_accepted_type = True

    # If input_type is an XML element
    else:
        type_name, type_kw, type_id = findType(input_type)
        type_el = cfg.id_dict[type_id]

        if type_el.tag in ['Class', 'Struct']:
            demangled_name = type_el.get('demangled')
            if demangled_name in cfg.accepted_types:
                is_accepted_type = True

        elif type_el.tag=='FundamentalType':
            type_name = type_el.get('name')
            if type_name in cfg.accepted_types:
                is_accepted_type = True

        else:
            warnings.warn('Cannot determine if XML element with id="%s" and tag "%s" corresponds to an accepted type. Assuming it does not.' % (input_type.get('id'), input_type.tag))
            is_accepted_type = False

    return is_accepted_type

# ====== END: isAcceptedType ========



# ====== isLoadedClass ========

def isLoadedClass(input_type, byname=False):

    is_loaded_class = False

    if byname:
        type_name = input_type

        # Remove '*' and '&'
        type_name = type_name.replace('*','').replace('&','')

        # Remove template bracket
        type_name = type_name.split('<')[0]

        # Check against cfg.loaded_classes
        if type_name in cfg.loaded_classes:
            is_loaded_class = True

    else:
        type_name, type_kw, type_id = findType(input_type)
        type_el = cfg.id_dict[type_id]

        if type_el.tag in ['Class', 'Struct']:
            demangled_name = type_el.get('demangled')
            if demangled_name in cfg.loaded_classes:
                is_loaded_class = True

    return is_loaded_class

# ====== END: isLoadedClass ========



# # ====== constrForwardDecls ========

# def constrForwardDecls():

#     import modules.classutils as classutils

#     code = ''
#     current_namespaces = []

#     for class_name_long, class_el in cfg.class_dict.items():

#         # print [class_name_full], [class_name_full.split('<',1)[0]], [class_name_full.split('<',1)[0].rsplit('::',1)[-1]]
#         namespaces    = getNamespaces(class_el)
#         has_namespace = bool(len(namespaces))
#         namespace_str = '::'.join(namespaces) + '::'*has_namespace

#         class_name       = classutils.getClassNameDict(class_el)
#         abstr_class_name = classutils.getClassNameDict(class_el, abstract=True)

#         if namespaces != current_namespaces:
#             # close current namespace
#             code += constrNamespace(current_namespaces, 'close', indent=cfg.indent)
#             # open new namespace
#             code += constrNamespace(namespaces, 'open', indent=cfg.indent)
#             # update current namespace
#             current_namespaces = namespaces

#         # Add forward decls
#         n_indents   = len(namespaces)
#         full_indent = ' '*n_indents*cfg.indent
        
#         if '<' in class_name['long_templ']:
#             is_template = True
#         else:
#             is_template = False

#         if is_template:
#             template_bracket = getTemplateBracket(class_el)[0]

#         if is_template:
#             code += full_indent + 'template ' + template_bracket + '\n'
#             code += full_indent + 'class ' + abstr_class_name['short'] + ';\n'
            
#             code += full_indent + 'template ' + template_bracket + '\n'
#             code += full_indent + 'class ' + class_name['short'] + ';\n'
            
#             code += full_indent + 'class ' + abstr_class_name['short_templ'] + ';\n'
#             code += full_indent + 'class ' + class_name['short_templ'] + ';\n'
#         else:
#             code += full_indent + 'class ' + abstr_class_name['short_templ'] + ';\n'
#             code += full_indent + 'class ' + class_name['short_templ'] + ';\n'
    
#     # Close current namespace
#     code += constrNamespace(current_namespaces, 'close', indent=cfg.indent)
#     code += '\n'

#     return code

# # ====== END: constrForwardDecls ========



# ====== constrForwardDeclHeader ========

def constrForwardDeclHeader():

    import modules.classutils as classutils

    code = ''
    insert_pos = 0

    current_namespaces = []

    for class_name_long, class_el in cfg.class_dict.items():

        # print [class_name_full], [class_name_full.split('<',1)[0]], [class_name_full.split('<',1)[0].rsplit('::',1)[-1]]
        namespaces    = getNamespaces(class_el)
        has_namespace = bool(len(namespaces))
        namespace_str = '::'.join(namespaces) + '::'*has_namespace

        # class_name       = classutils.getClassNameDict(class_el)
        abstr_class_name = classutils.getClassNameDict(class_el, abstract=True)

        if namespaces != current_namespaces:
            # close current namespace
            code += constrNamespace(current_namespaces, 'close', indent=cfg.indent)
            # open new namespace
            code += constrNamespace(namespaces, 'open', indent=cfg.indent)
            # update current namespace
            current_namespaces = namespaces

        # Add forward decls
        n_indents   = len(namespaces)
        full_indent = ' '*n_indents*cfg.indent
        
        if '<' in abstr_class_name['long_templ']:
            is_template = True
        else:
            is_template = False

        if is_template:
            template_bracket = getTemplateBracket(class_el)[0]

        if is_template:
            code += full_indent + 'template ' + template_bracket + '\n'
            code += full_indent + 'class ' + abstr_class_name['short'] + ';\n'
            
            # code += full_indent + 'template ' + template_bracket + '\n'
            # code += full_indent + 'class ' + class_name['short'] + ';\n'
            
            code += full_indent + 'class ' + abstr_class_name['short_templ'] + ';\n'
            # code += full_indent + 'class ' + class_name['short_templ'] + ';\n'
        else:
            code += full_indent + 'class ' + abstr_class_name['short_templ'] + ';\n'
            # code += full_indent + 'class ' + class_name['short_templ'] + ';\n'
    
    # Close current namespace
    code += constrNamespace(current_namespaces, 'close', indent=cfg.indent)
    code += '\n'

    code_tuple = (insert_pos, code)

    return code_tuple

# ====== END: constrForwardDeclHeader ========



# ====== constrTypedefHeader ========

def constrTypedefHeader():

    typedef_code = ''
    insert_pos   = 0

    # # - Forward declarations
    # typedef_code += '// Forward declarations:\n'
    # typedef_code += constrForwardDecls()

    # - Typedefs
    typedef_code += '// Typedefs:\n'
    for typedef_name, typedef_el in cfg.typedef_dict.items():
        type_name, type_kw, type_id = findType(typedef_el)
        typedef_code += 'typedef ' + type_name + ' ' + typedef_name + ';\n'

    code_tuple = (insert_pos, typedef_code)

    return code_tuple

# ====== END: constrTypedefHeader ========



# ====== getParentClasses ========

def getParentClasses(class_el, only_native_classes=False, only_loaded_classes=False):

    import modules.classutils as classutils

    parent_classes = []

    sub_el_list = class_el.findall('Base')
    for sub_el in sub_el_list:

        base_id = sub_el.get('type')
        base_el = cfg.id_dict[base_id]

        if (only_loaded_classes) and (not isLoadedClass(base_el)):
            continue
        elif (only_native_classes) and (not isNative(base_el)):
            continue
        else:
            base_access    = sub_el.get('access')
            base_virtual   = bool( int( sub_el.get('virtual') ) )
            
            base_name_dict       = classutils.getClassNameDict(base_el)
            abstr_base_name_dict = classutils.getClassNameDict(base_el, abstract=True)

            is_accepted_type = isAcceptedType(base_el)
            is_native        = isNative(base_el)
            is_fundamental   = isFundamental(base_el)
            is_std           = isStdType(base_el)
            is_loaded_class  = isLoadedClass(base_el)


            temp_dict = {}
            temp_dict['class_name']       = base_name_dict
            temp_dict['abstr_class_name'] = abstr_base_name_dict
            temp_dict['wrapper_name']     = classutils.toWrapperType(base_name_dict['long'])
            temp_dict['access']           = base_access
            temp_dict['virtual']          = base_virtual
            temp_dict['id']               = base_id

            temp_dict['accepted']         = is_accepted_type
            temp_dict['native']           = is_native
            temp_dict['fundamental']      = is_fundamental
            temp_dict['std']              = is_std
            temp_dict['loaded']           = is_loaded_class

            parent_classes.append(temp_dict)

    return parent_classes

# ====== END: getParentClasses ========



# ====== getAllParentClasses ========

def getAllParentClasses(class_el, only_native_classes=True, only_loaded_classes=False):

    parent_classes = []

    temp_class_list = [class_el]
    while len(temp_class_list) > 0:
        current_class = temp_class_list.pop()
        if 'bases' in current_class.keys():
            for parent_class_id in current_class.get('bases').split():

                # Remove accessor info from id, e.g. "private:_123" --> "_123"
                parent_class_id = parent_class_id.split(':')[-1]

                parent_class_el = cfg.id_dict[parent_class_id]

                if only_loaded_classes and not isLoadedClass(parent_class_el):
                    continue
                elif only_native_classes and not isNative(parent_class_el):
                    continue
                else:
                    parent_classes.append(parent_class_el)
                    temp_class_list.append(parent_class_el)

    return parent_classes

# ====== END: getAllParentClasses ========



# ====== getAllTypesInClass ========

def getAllTypesInClass(class_el, include_parents=False):

    import modules.classutils as classutils
    import modules.funcutils as funcutils

    all_types = []

    check_member_elements = getMemberElements(class_el)

    class_id = class_el.get('id')
    for mem_el in check_member_elements:

        if mem_el.tag in ['Constructor', 'Destructor', 'Method', 'OperatorMethod']:
            args_list = funcutils.getArgs(mem_el)
            for arg_dict in args_list:

                arg_type_el   = cfg.id_dict[arg_dict['id']]
                arg_type_name = classutils.getClassNameDict(arg_type_el)                

                arg_type_dict = {}
                arg_type_dict['class_name'] = arg_type_name
                arg_type_dict['el']         = arg_type_el

                all_types.append(arg_type_dict)

        if ('type' in mem_el.keys()) or ('returns' in mem_el.keys()):

            mem_type, mem_type_kw, mem_type_id = findType(mem_el)

            type_el   = cfg.id_dict[mem_type_id]
            type_name = classutils.getClassNameDict(type_el)

            type_dict = {}
            type_dict['class_name'] = type_name
            type_dict['el']         = type_el

            all_types.append(type_dict)


    if include_parents:
        parent_classes = getParentClasses(class_el, only_native_classes=False, only_loaded_classes=False)

        for parent_dict in parent_classes:

            small_parent_dict = {}
            small_parent_dict['class_name'] = parent_dict['class_name']
            small_parent_dict['el']         = cfg.id_dict[parent_dict['id']]

            all_types.append(small_parent_dict)

    return all_types

# ====== END: getAllTypesInClass ========



# ====== getMemberElements ========

def getMemberElements(el, include_artificial=False):

    member_elements = []

    if 'members' in el.keys():
        for mem_id in el.get('members').split():
            mem_el = cfg.id_dict[mem_id]
            if include_artificial:
                member_elements.append(mem_el)
            else:
                if not 'artificial' in mem_el.keys():
                    member_elements.append(mem_el)

    return member_elements

# ====== END: getMemberElements ========



# ====== getMemberFunctions ========

def getMemberFunctions(class_el, include_artificial=False, include_inherited=False, only_accepted=True, limit_pointerness=True, include_operators=False):

    import modules.funcutils as funcutils

    all_classes   = [class_el]
    all_members   = []
    all_functions = []

    # If include_inherited=True, append all (native) parent classes 
    # the list 'all_classes'
    if include_inherited:
        parent_classes = getAllParentClasses(class_el, only_loaded_classes=True)
        all_classes = all_classes + parent_classes

        # temp_class_list = list(all_classes)
        # while len(temp_class_list) > 0:
        #     current_class = temp_class_list.pop()
        #     if 'bases' in current_class.keys():
        #         for parent_class_id in current_class.get('bases').split():
        #             parent_class_el = cfg.id_dict[parent_class_id]
        #             if isLoadedClass(parent_class_el):
        #                 all_classes.append(parent_class_el)
        #                 temp_class_list.append(parent_class_el)

    # Get all member elements
    for el in all_classes:
        class_members = getMemberElements(el, include_artificial=include_artificial)
        all_members = all_members + class_members

    # Extract only regular member functions (no variables, constructors, destructors, ...)
    for mem_el in all_members:
        if (mem_el.tag == 'Method' or (include_operators==True and mem_el.tag=='OperatorMethod')) and (mem_el.get('access') == 'public'):

            if only_accepted and funcutils.ignoreFunction(mem_el, limit_pointerness=limit_pointerness):
                method_name = mem_el.get('name')
                if mem_el.tag=='OperatorMethod':
                    method_name = 'operator' + method_name
                # warnings.warn('The member "%s", belonging to (or inherited by) class "%s", makes use of a non-accepted type and will be ignored.' % (method_name, class_el.get('name')))
                print 'INFO: ' + 'The member "%s", belonging to (or inherited by) class "%s", makes use of a non-accepted type and will be ignored.' % (method_name, class_el.get('name'))

                continue

            else:
                all_functions.append(mem_el)

    return all_functions

# ====== END: getMemberFunctions ========




# ====== constrNamespace ========

def constrNamespace(namespaces, open_or_close, indent=cfg.indent):

    code = ''

    if open_or_close == 'open':
        n_indents = 0
        for ns in namespaces:
            code += ' '*n_indents*indent + 'namespace ' + ns + '\n'
            code += ' '*n_indents*indent + '{' + '\n'
            n_indents += 1

    elif open_or_close == 'close':
        n_indents = len(namespaces)
        for ns in namespaces:
            n_indents -= 1
            code += ' '*n_indents*indent + '}' + '\n'

    else:
        raise ValueError("Second argument must be either 'open' or 'close'.")

    return code

# ====== END: constrNamespace ========



# ====== pointerAndRefCheck ========

def pointerAndRefCheck(input_type, byname=False):

    #
    # Input type should either be an XML element (byname=False)
    # or a string (byname=True)
    #

    if byname:
        type_name = input_type
    else:
        type_name, type_kw, type_id = findType(input_type)

    # Remove template bracket
    if '<' in type_name:
        type_name = type_name.split('<',1)[0] + type_name.rsplit('>',1)[-1]

    # Check pointerness
    pointerness = type_name.count('*')

    # Check reference
    is_reference = bool('&' in type_name)

    return pointerness, is_reference

# ====== END: pointerAndRefCheck ========



# ====== addIncludeGuard ========

def addIncludeGuard(code, file_name):

    guard_var = '__' + file_name.replace('.','_').upper() + '__'

    guard_code_top    = '#ifndef ' + guard_var + '\n' + '#define ' + guard_var + '\n'
    guard_code_bottom = '#endif /* ' + guard_var + ' */\n'
    
    new_code = guard_code_top + '\n' + code + '\n' + guard_code_bottom

    return new_code 

# ====== END: addIncludeGuard ========



# ====== identifyIncludedHeaders ========

def identifyIncludedHeaders(content, only_native=True):
    
    return_dict = OrderedDict()

    # Remove comments
    content      = removeComments(content, insert_blanks=True)
    content_list = content.split('\n')

    # Parse content and identify header file names
    headers_in_file = []
    for line in content_list:

        line = line.strip()

        if line[0:8] == '#include':
            
            header_file_name = line.split()[1]

            # Skip standard headers (of the form: #include <FILENAME>)
            if header_file_name[0] == '<':
                continue
            else:
                header_file_name = header_file_name.strip('"')
                headers_in_file.append(os.path.basename(header_file_name))
        else:
            continue

    # Connect with XML elements
    for check_file_path, file_el in cfg.file_dict.items():

        # - If only_native=True, check for match with cfg.accepted_paths
        if only_native:
            is_native_file = False
            for accepted_path in cfg.accepted_paths:
                if accepted_path in os.path.dirname(check_file_path):
                    is_native_file = True
                    break
            if not is_native_file:
                continue

        # - Cut down to file name only
        check_file_name = os.path.basename(check_file_path)

        # - Keep XML id if the corresponding file name mathces with an identified header
        if check_file_name in headers_in_file:
            return_dict[check_file_name] = file_el.get('id')

    return return_dict

# ====== END: identifyIncludedHeaders ========




