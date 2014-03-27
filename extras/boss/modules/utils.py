################################
#                              #
#  Utility functions for BOSS  #
#                              #
################################

import xml.etree.ElementTree as ET
from collections import OrderedDict
from operator import itemgetter
import os

import modules.cfg as cfg


# ====== isFundamental ========

def isFundamental(el):

    is_fundamental = False

    if el.tag == 'FundamentalType':
        is_fundamental = True

    return is_fundamental

# ====== END: isNative ========


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

    print
    print 'AT END OF LINE: ', file_content_nocomments[prev_pos:newline_pos]
    print
    # # HACK
    # # Then, find position of class name
    # short_class_name = el.get('name').rsplit('::')[-1]
    # print 'short class name: ', short_class_name
    # search_limit = newline_pos
    # while search_limit > -1:
    #     pos = file_content_nocomments[:search_limit].rfind(short_class_name)
    #     pre_char  = file_content_nocomments[pos-1]
    #     post_char = file_content_nocomments[pos+len(short_class_name)]
    #     if (pre_char in [' ','\n','\t']) and (post_char in [' ', ':', '\n', '<', '{']):
    #         break
    #     else:
    #         search_limit = pos

    # class_name_pos = pos


    # Find the template parameter bracket, e.g. <typename A, typename B>
    search_content = file_content_nocomments[:newline_pos]
    # search_content = file_content_nocomments[:class_name_pos]
    # HACK END
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
        spec_types = el.get('name').split('<',1)[1].rsplit('>',1)[0]
    elif el.tag == 'Function':
        spec_types = el.get('demangled').split('<',1)[1].rsplit('>',1)[0]
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

def findType(el):

    check_keywords = ['static', 'const']
    additional_keywords = []
    refs_and_pointers = ''

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
        modified_name = full_name.replace(' ','').split('<')[0]

        namespaces = modified_name.split('::')[:-1]

    return namespaces

# ====== END: getNamespaces ========



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
            demangled_name = type_el.get('demangled').replace(' ','')
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



# ====== constrForwardDecls ========

def constrForwardDecls():

    import modules.classutils as classutils

    code = ''
    current_namespaces = []

    for class_name_full, class_el in cfg.class_dict.items():

        # print [class_name_full], [class_name_full.split('<',1)[0]], [class_name_full.split('<',1)[0].rsplit('::',1)[-1]]
        namespaces    = getNamespaces(class_el)
        has_namespace = bool(len(namespaces))
        namespace_str = '::'.join(namespaces) + '::'*has_namespace

        class_name_short                  = class_name_full.replace(namespace_str,'',1)
        class_name_short_notemplate       = class_name_short.split('<',1)[0]
        abstr_class_name_short            = classutils.getAbstractClassName(class_name_full, short=True)
        abstr_class_name_short_notemplate = abstr_class_name_short.split('<',1)[0]

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
        
        if '<' in class_name_full:
            is_template = True
        else:
            is_template = False

        if is_template:
            template_bracket = utils.getTemplateBracket(class_el)[0]

        if is_template:
            code += full_indent + 'template ' + template_bracket + '\n'
            code += full_indent + 'class ' + abstr_class_name_short_notemplate + ';\n'
            
            code += full_indent + 'template ' + template_bracket + '\n'
            code += full_indent + 'class ' + class_name_short_notemplate + ';\n'
            
            code += full_indent + 'class ' + abstr_class_name_short + ';\n'
            code += full_indent + 'class ' + class_name_short + ';\n'
        else:
            code += full_indent + 'class ' + abstr_class_name_short + ';\n'
            code += full_indent + 'class ' + class_name_short + ';\n'
    
    # Close current namespace
    code += constrNamespace(current_namespaces, 'close', indent=cfg.indent)
    code += '\n'

    return code

# ====== END: constrForwardDecls ========



# ====== constrTypedefHeader ========

def constrTypedefHeader():

    typedef_code = ''
    insert_pos   = 0

    # - Forward declarations
    typedef_code += '// Forward declarations:\n'
    typedef_code += constrForwardDecls()

    # - Typedefs
    typedef_code += '// Typedefs:\n'
    for typedef_name, typedef_el in cfg.typedef_dict.items():
        type_name, type_kw, type_id = findType(typedef_el)
        typedef_code += 'typedef ' + type_name + ' ' + typedef_name + ';\n'

    code_tuple = (insert_pos, typedef_code)

    return code_tuple

# ====== END: constrTypedefHeader ========



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



