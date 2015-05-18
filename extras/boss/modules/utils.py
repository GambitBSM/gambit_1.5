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
import subprocess

import modules.cfg as cfg
import modules.gb as gb
import modules.shelltimeout as shelltimeout
import modules.exceptions as exceptions

import inspect


# ====== isComplete ========

def isComplete(class_el):

    is_complete = True

    # Check that class is complete (not only forward declared):
    if ('incomplete' in class_el.keys()) and (class_el.get('incomplete') == '1'):
        is_complete = False

    return is_complete

# ====== END: isComplete ========



# ====== isLoadable ========

def isLoadable(class_el, print_warning=False):

    import modules.classutils as classutils

    class_name  = classutils.getClassNameDict(class_el) 
    is_loadable = True

    # - Check that class is complete (not only forward declared):
    if not isComplete(class_el):
        is_loadable = False
        if print_warning:
            print 'INFO: ' + 'Class %s is not complete.' % (class_name['long_templ'])

    # - Check for pure virtual members:
    pure_virtual_members = classutils.pureVirtualMembers(class_el)
    if len(pure_virtual_members) > 0:
        is_loadable = False
        if print_warning:
            print 'INFO: ' + 'The following member functions in class %s are pure virtual: %s' % (class_name['long_templ'], ', '.join(pure_virtual_members))
    
    if (not is_loadable) and (print_warning):
        print 'INFO: ' + 'Class %s cannot be loaded and will be ignored.' % (class_name['long_templ'])

    return is_loadable

# ====== END: isLoadable ========



# ====== isFundamental ========

def isFundamental(el):

    is_fundamental = False

    if el.tag == 'FundamentalType':
        is_fundamental = True

    return is_fundamental

# ====== END: isFundamental ========



# ====== isEnumeration ========

def isEnumeration(el):

    is_enumeration = False

    if el.tag == 'Enumeration':
        is_enumeration = True

    return is_enumeration

# ====== END: isEnumeration ========



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
            file_el = gb.id_dict[el.get('file')]

        check_path = file_el.get('name')

        is_native = False
        for accepted_path in cfg.accepted_paths:
            if accepted_path in os.path.dirname(check_path):
                is_native = True
                break

    else:
        raise Exception('Cannot check whether XML element with id="%s" and tag "%s" is native.' % (el.get('id'), el.tag))

    return is_native

# ====== END: isNative ========



# ====== isStdType ========

def isStdType(el):

    # Makes use of global variables:  accepted_paths

    is_std = False
    can_check_tags = ['Class', 'Struct', 'Union', 'Enumeration']

    if el.tag in can_check_tags:

        if 'demangled' in el.keys():
            demangled_name = el.get('demangled')
            if demangled_name[0:5] == 'std::':
                is_std = True

    else:
        is_std = False
        # warnings.warn('Cannot check whether XML element with id="%s" and tag "%s" is a standard type.' % (el.get('id'), el.tag))
        # raise Exception('Cannot check whether XML element with id="%s" and tag "%s" is a standard type.' % (el.get('id'), el.tag))

    return is_std

# ====== END: isStdType ========



# ====== isConstFunction ========

def isConstFunction(func_el):

    is_const_func = False

    if ('const' in func_el.keys()) and (func_el.get('const')=='1'):
        is_const_func = True

    return is_const_func

# ====== END: isConstFunction ========



# ====== getTemplateBracket ========

def getTemplateBracket(el):

    src_file_name = gb.id_dict[el.get('file')].get('name')
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

def getSpecTemplateTypes(input_type, byname=False):

    # If input is a string
    if byname:
        input_name = input_type

    # If input is an xml element
    else:
        el = input_type

        # Classes and functions must be treated differently
        if el.tag in ['Class', 'Struct']:
            input_name = el.get('name')
        elif el.tag in ['Function', 'Method', 'OperatorMethod', 'OperatorFunction']:
            input_name = el.get('demangled')
        else:
            raise Exception("Don't know how to get template types from XML element with tag: %s" % el.tag)

    input_name_no_templ, templ_bracket = removeTemplateBracket(input_name, return_bracket=True)
    spec_types = templ_bracket.strip().lstrip('<').rstrip('>').strip()


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
        spec_types_list.append(spec_types[prev_p:p].strip())
        prev_p = p+1
    spec_types_list.append(spec_types[prev_p:].strip())

    # Return result
    if spec_types_list == ['']:
        return []
    else:
        return spec_types_list

# ====== END: getSpecTemplateTypes ========



# ====== unpackAllSpecTemplateTypes ========

def unpackAllSpecTemplateTypes(input_bracket, result_list):

    spec_types = getSpecTemplateTypes(input_bracket, byname=True)

    for type_name in spec_types:
        
        if '<' in type_name:
            unpackAllSpecTemplateTypes(type_name, result_list=result_list)
            
        else:
            result_list.append(type_name)

# ====== END: unpackAllSpecTemplateTypes ========



# ====== getBasicTypeName ========

def getBasicTypeName(type_name):

    # If type name contains a template brackets
    if '<' in type_name:
        type_name_notempl, templ_bracket = removeTemplateBracket(type_name, return_bracket=True)
        before_bracket, after_bracket = type_name.split(templ_bracket)

        if (len(after_bracket) > 0) and (after_bracket[0] == ' '):
            space_after_bracket = True
        else:
            space_after_bracket = False

        # Remove asterix and/or ampersand
        before_bracket = before_bracket.replace('*', '').replace('&', '')
        after_bracket  = after_bracket.replace('*', '').replace('&', '')

        # Remove 'const' and 'volatile'
        before_bracket_list = before_bracket.split()
        before_bracket_list = [item for item in before_bracket_list if item != 'const']
        before_bracket_list = [item for item in before_bracket_list if item != 'volatile']
        before_bracket = ' '.join(before_bracket_list)

        after_bracket_list = after_bracket.split()
        after_bracket_list = [item for item in after_bracket_list if item != 'const']
        after_bracket_list = [item for item in after_bracket_list if item != 'volatile']
        after_bracket = ' '.join(after_bracket_list)

        basic_type_name = before_bracket + templ_bracket + ' '*space_after_bracket + after_bracket

    # If no template bracket
    else:

        basic_type_name = type_name

        # Remove asterix and/or ampersand
        basic_type_name = basic_type_name.replace('*', '').replace('&', '')
        
        # Remove 'const' and 'volatile'
        basic_type_name_list = basic_type_name.split()
        basic_type_name_list = [item for item in basic_type_name_list if item != 'const']
        basic_type_name_list = [item for item in basic_type_name_list if item != 'volatile']
        basic_type_name = ' '.join(basic_type_name_list)

    # Return result
    return basic_type_name

# ====== END: getBasicTypeName ========


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

    # check_keywords = ['const']
    # additional_keywords = []

    cv_qualifiers = []
    is_reference = False
    pointerness = 0
    found_function_pointer = False
    is_array = False
    array_limits = []

    el = el_input
    
    if el.tag in ['FundamentalType', 'Class', 'Struct', 'Enumeration']:
        type_id = el.get('id')

    elif el.tag in ['Constructor']:
        type_id = el.get('context')
        el = gb.id_dict[type_id]

    else:
        type_id = el.get('type')
        prev_tag = ''
        while ('type' in el.keys()) or ('returns' in el.keys()):

            # Get xml id to move further through the xml file
            if el.tag in ['FunctionType', 'Function', 'Method', 'OperatorMethod']:

                if (el.tag == 'FunctionType') and (prev_tag == 'PointerType'):
                    found_function_pointer = True
                type_id = el.get('returns')
            else:
                type_id = el.get('type')

            # Check for reference or pointer type
            if el.tag == 'ReferenceType':
                is_reference = True
            if el.tag == 'PointerType':
                pointerness += 1

            # Pick up any extra keywords (e.g. 'const') from non-function xml elements
            if el.tag == 'CvQualifiedType':
                if 'const' in el.keys():
                    if 'const' not in cv_qualifiers:
                        cv_qualifiers.append('const')
                if 'volatile' in el.keys():
                    if 'volatile' not in cv_qualifiers:
                        cv_qualifiers.append('volatile')

            # Pick up any array indices
            if el.tag == 'ArrayType':
                is_array = True
                if el.get('max') != '':
                    max_index_str = el.get('max')
                    max_index = int(max_index_str.strip('u'))
                    array_limits.append(max_index + 1)


            # Store tag (to identify function pointers)
            prev_tag = el.tag

            # change xlm element 'el'
            el = gb.id_dict[type_id]

    # When we exit the loop, 'el' is at the final element.
    # Now get the full name, including any namespaces.
    name_and_namespaces = getNamespaces(el, include_self=True)
    typename = '::'.join(name_and_namespaces)

    # if check_function_pointer:
    #     return typename, additional_keywords, type_id, found_function_pointer
    # else:
    #     return typename, additional_keywords, type_id


    type_dict = {}
    type_dict['name']                = typename
    type_dict['cv_qualifiers']       = cv_qualifiers
    type_dict['is_reference']        = is_reference
    type_dict['pointerness']         = pointerness
    type_dict['id']                  = type_id
    type_dict['el']                  = el
    type_dict['is_function_pointer'] = found_function_pointer
    type_dict['is_array']            = is_array
    type_dict['array_limits']        = tuple(array_limits)

    return type_dict

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

    if indent == 0:
        new_content = content

    else:
        lines = content.split('\n')

        new_content = '\n'.join( [' '*indent + line for line in lines] )

        # If the last char in content was a newline, 
        # remove the indentation that was added after that newline
        if lines[-1] == '':
            new_content = new_content[:-indent]

    return new_content

# ====== END: addIndentation ========



# ====== getNamespaces ========

def getNamespaces(xml_el, include_self=False):

    namespaces = []

    if include_self:
        namespaces.append(xml_el.get('name'))

    current_xml_el = xml_el
    while 'context' in current_xml_el.keys():
        context_id = current_xml_el.get('context')
        context_xml_el = gb.id_dict[context_id]
        context_name = context_xml_el.get('name')
        namespaces.append(context_name)
        current_xml_el = context_xml_el

    namespaces.reverse()

    # If present, remove the first, default namespace in the XML file from the list (with name='::')
    if (len(namespaces) > 0) and (namespaces[0] == '::'):
        namespaces = namespaces[1:]

    return namespaces

# ====== END: getNamespaces ========



# ====== removeTemplateBracket ========

def removeTemplateBracket(type_name, return_bracket=False):
    
    if '<' in type_name:
        # type_name_notempl = type_name.split('<',1)[0] + type_name.rsplit('>',1)[-1]
        l_pos = type_name.find('<')
        r_pos = type_name.rfind('>')
        
        type_name_notempl = type_name[:l_pos] + type_name[r_pos+1:]
        template_bracket  = type_name[l_pos:r_pos+1]
    else:
        type_name_notempl = type_name
        template_bracket  = ''

    if return_bracket:
        return type_name_notempl, template_bracket
    else:
        return type_name_notempl

# ====== END: removeTemplateBracket ========



# ====== removeNamespace ========

def removeNamespace(type_name, return_namespace=False):
    
    type_name_notempl = removeTemplateBracket(type_name)

    if '::' in type_name_notempl:
        namespace = type_name_notempl.rsplit('::',1)[0]
        new_type_name = type_name.replace(namespace+'::','',1)
    else:
        new_type_name = type_name
        namespace = ''

    if return_namespace:
        return namespace, new_type_name
    else:
        return new_type_name

# ====== END: removeNamespace ========



# ====== isAcceptedType ========

def isAcceptedType(input, byname=False):
    #
    # By default, input is assumed to be an XML element (xml.etree.ElementTree.Element).
    # If byname=True, input is a simple string.
    #

    is_accepted_type = False


    # If input is a string
    if byname:

        base_input_name = getBasicTypeName(input)
     
        if base_input_name in gb.accepted_types:
            is_accepted_type = True


    # If input is an XML element
    else:

        type_dict = findType(input)
        type_el      = type_dict['el']
        type_name    = type_dict['name']
        pointerness  = type_dict['pointerness']
        is_ref       = type_dict['is_reference']
        is_array     = type_dict['is_array']
        array_limits = type_dict['array_limits']


        # # BOSS cannot yet handle pointers to loaded types
        # if (pointerness > 0) and isLoadedClass(type_el):
        #     print 'INFO: ' + 'Problem with "%s". BOSS cannot yet handle pointer-to-loaded-type.' % (input.get('name'))
        #     is_accepted_type = False
        #     return is_accepted_type

        # BOSS cannot yet handle delarations of arrays with unspecified length
        if (is_array) and (len(array_limits) == 0):
            print 'INFO: ' + 'Problem with "%s". BOSS cannot yet handle arrays declared with unspecified length.' % (input.get('name'))
            is_accepted_type = False
            return is_accepted_type

        # BOSS cannot yet handle arrays of loaded types
        if isLoadedClass(type_el) and is_array:
            print 'INFO: ' + 'Problem with "%s". BOSS cannot yet handle arrays of a loaded type.' % (input.get('name'))
            is_accepted_type = False
            return is_accepted_type

        # BOSS cannot yet handle function pointers
        if type_dict['is_function_pointer']:
            print 'INFO: ' + 'Problem with "%s". BOSS cannot yet handle function pointers.' % (input.get('name'))
            is_accepted_type = False
            return is_accepted_type


        if type_el.tag in ['Class', 'Struct']:
            demangled_name = type_el.get('demangled')
            if demangled_name in gb.accepted_types:
                is_accepted_type = True

        elif type_el.tag in ['FundamentalType', 'Enumeration']:
            type_name = type_el.get('name')
            if type_name in gb.accepted_types:
                is_accepted_type = True

        else:
            warnings.warn('Cannot determine if XML element with id="%s" and tag "%s" corresponds to an accepted type. Assuming it does not.' % (input.get('id'), input.tag))
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
        type_dict = findType(input_type)
        type_el = type_dict['el']

        if type_el.tag in ['Class', 'Struct']:
            demangled_name = type_el.get('demangled')
            if demangled_name in cfg.loaded_classes:
                is_loaded_class = True

    return is_loaded_class

# ====== END: isLoadedClass ========



# ====== constrAbsForwardDeclHeader ========

def constrAbsForwardDeclHeader(file_output_path):

    import modules.classutils as classutils

    # If this is the first time this function is executed, read initial code from header_templates/ folder
    if file_output_path not in gb.new_code.keys():
        f = open('header_templates/standard_header_template.hpp')
        initial_code = f.read()
        f.close()
        initial_code_tuple = (0, initial_code)
        gb.new_code[file_output_path] = {'code_tuples':[initial_code_tuple], 'add_include_guard':True}


    current_code = gb.new_code[file_output_path]['code_tuples'][0][1]

    insert_code = ''
    tag_pos = current_code.find('__INSERT_CODE_HERE__')

    current_namespaces = []
    for class_name_long, class_el in gb.loaded_classes_in_xml.items():

        # print [class_name_full], [class_name_full.split('<',1)[0]], [class_name_full.split('<',1)[0].rsplit('::',1)[-1]]
        namespaces    = getNamespaces(class_el)
        has_namespace = bool(len(namespaces))
        namespace_str = '::'.join(namespaces) + '::'*has_namespace

        # class_name       = classutils.getClassNameDict(class_el)
        abstr_class_name = classutils.getClassNameDict(class_el, abstract=True)

        if namespaces != current_namespaces:
            # close current namespace
            insert_code += constrNamespace(current_namespaces, 'close', indent=cfg.indent)
            # open new namespace
            insert_code += constrNamespace(namespaces, 'open', indent=cfg.indent)
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
            insert_code += full_indent + 'template ' + template_bracket + '\n'
            insert_code += full_indent + 'class ' + abstr_class_name['short'] + ';\n'
            insert_code += full_indent + 'class ' + abstr_class_name['short_templ'] + ';\n'
        else:
            insert_code += full_indent + 'class ' + abstr_class_name['short_templ'] + ';\n'
    
    # Close current namespace
    insert_code += constrNamespace(current_namespaces, 'close', indent=cfg.indent)
    insert_code += '\n'

    new_code = current_code[:tag_pos] + insert_code + current_code[tag_pos:]
    
    
    # Replace other code tags
    new_code = replaceCodeTags(new_code)

    new_code_tuple = (0,new_code)
    
    # Overwrite existing code tuple
    gb.new_code[file_output_path]['code_tuples'] = [(new_code_tuple)]

# ====== END: constrAbsForwardDeclHeader ========



# ====== constrWrpForwardDeclHeader ========

def constrWrpForwardDeclHeader(file_output_path):

    import modules.classutils as classutils

    # If this is the first time this function is executed, read initial code from header_templates/ folder
    if file_output_path not in gb.new_code.keys():
        f = open('header_templates/standard_header_template.hpp')
        initial_code = f.read()
        f.close()
        initial_code_tuple = (0, initial_code)
        gb.new_code[file_output_path] = {'code_tuples':[initial_code_tuple], 'add_include_guard':True}

    current_code = gb.new_code[file_output_path]['code_tuples'][0][1]

    insert_code = ''
    tag_pos = current_code.find('__INSERT_CODE_HERE__')

    for class_name in gb.loaded_classes_in_xml.keys():

        namespace, class_name_short = removeNamespace(class_name, return_namespace=True)

        if namespace == '':
            namespace_list = []
        else:
            namespace_list = namespace.split('::')

        n_indents = len(namespace_list)
        

        # - Open namespace
        insert_code += constrNamespace(namespace_list, 'open')

        # - Forward declaration
        insert_code += ' '*n_indents*cfg.indent + 'class ' + class_name_short + ';\n'

        # - Close namespace
        insert_code += constrNamespace(namespace_list, 'close')


    new_code = current_code[:tag_pos] + insert_code + current_code[tag_pos:]
    

    # Replace other code tags
    new_code = replaceCodeTags(new_code)


    new_code_tuple = (0,new_code)
    
    # Overwrite existing code tuple
    gb.new_code[file_output_path]['code_tuples'] = [(new_code_tuple)]

# ====== END: constrWrpForwardDeclHeader ========



# # ====== constrWrpForwardDeclHeader ========

# def constrWrpForwardDeclHeader():

#     import modules.classutils as classutils

#     code = ''
#     insert_pos = 0

#     for class_name in cfg.loaded_classes:

#         namespace, class_name_short = removeNamespace(class_name, return_namespace=True)

#         if namespace == '':
#             namespace_list = []
#         else:
#             namespace_list = namespace.split('::')

#         n_indents = len(namespace_list)
        

#         # - Open namespace
#         code += constrNamespace(namespace_list, 'open')

#         # - Forward declaration
#         code += ' '*n_indents*cfg.indent + 'class ' + class_name_short + ';\n'

#         # - Close namespace
#         code += constrNamespace(namespace_list, 'close')



#     # Insert tags for the GAMBIT namespace
#     code = '\n__START_GAMBIT_NAMESPACE__\n\n' + code + '\n__END_GAMBIT_NAMESPACE__\n'


#     # Insert include statements needed by GAMBIT
#     code = '\n' + '#include "identification.hpp"\n' + code + '\n#include "backend_undefs.hpp"\n'


#     code_tuple = (insert_pos, code)

#     return code_tuple

# # ====== END: constrWrpForwardDeclHeader ========



# ====== constrTypedefHeader ========

def constrTypedefHeader():

    typedef_code = ''
    insert_pos   = 0

    # - Typedefs
    typedef_code += '// Typedefs:\n'
    for typedef_name, typedef_el in gb.typedef_dict.items():

        type_dict = findType(typedef_el)
        type_name = type_dict['name'] + '*'*type_dict['pointerness'] + '&'*type_dict['is_reference']

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
        base_el = gb.id_dict[base_id]

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

                parent_class_el = gb.id_dict[parent_class_id]

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

                arg_type_el   = gb.id_dict[arg_dict['id']]
                arg_type_name = classutils.getClassNameDict(arg_type_el)                

                arg_type_dict = {}
                arg_type_dict['class_name'] = arg_type_name
                arg_type_dict['el']         = arg_type_el

                all_types.append(arg_type_dict)

        if ('type' in mem_el.keys()) or ('returns' in mem_el.keys()):

            mem_type_dict = findType(mem_el)
            type_el = mem_type_dict['el'] 
            
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
            small_parent_dict['el']         = gb.id_dict[parent_dict['id']]

            all_types.append(small_parent_dict)

    return all_types

# ====== END: getAllTypesInClass ========



# ====== getMemberElements ========

def getMemberElements(el, include_artificial=False):

    member_elements = []

    if 'members' in el.keys():
        for mem_id in el.get('members').split():
            mem_el = gb.id_dict[mem_id]

            # Check if this member element should be ditched
            if mem_el.tag in ['Method', 'OperatorMethod', 'Field', 'Variable']:
                if 'name' in mem_el.keys():
                    namespaces_list = getNamespaces(mem_el, include_self=True)
                    full_name = '::'.join(namespaces_list)
                    if full_name in cfg.ditch:
                        print 'DITCH:', full_name
                        continue


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
        #             parent_class_el = gb.id_dict[parent_class_id]
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



# ====== getAllTypesInFunction ========

def getAllTypesInFunction(func_el):

    import modules.classutils as classutils
    import modules.funcutils as funcutils

    all_types = []

    func_id = func_el.get('id')

    if func_el.tag in ['Function', 'Constructor']:
        args_list = funcutils.getArgs(func_el)
        for arg_dict in args_list:

            arg_type_el   = gb.id_dict[arg_dict['id']]
            arg_type_name = classutils.getClassNameDict(arg_type_el)                

            arg_type_dict = {}
            arg_type_dict['class_name'] = arg_type_name
            arg_type_dict['el']         = arg_type_el

            all_types.append(arg_type_dict)

    if ('type' in func_el.keys()) or ('returns' in func_el.keys()) or (func_el.tag=='Constructor' and 'context' in func_el.keys()):

        mem_type_dict = findType(func_el)
        type_el = mem_type_dict['el'] 

        type_name = classutils.getClassNameDict(type_el)

        type_dict = {}
        type_dict['class_name'] = type_name
        type_dict['el']         = type_el

        all_types.append(type_dict)

    return all_types

# ====== END: getAllTypesInFunction ========



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

        # Remove template bracket
        if '<' in type_name:
            type_name = type_name.split('<',1)[0] + type_name.rsplit('>',1)[-1]

        # Check pointerness
        pointerness = type_name.count('*')

        # Check reference
        is_reference = bool('&' in type_name)

    else:
        type_dict = findType(input_type)
        pointerness  = type_dict['pointerness']
        is_reference = type_dict['is_reference']


    return pointerness, is_reference

# ====== END: pointerAndRefCheck ========



# ====== addIncludeGuard ========

def addIncludeGuard(code, file_name, prefix='', suffix='', uppercase=False):

    # if suffix == '':
    #     guard_var = '__' + file_name.replace('.','_') + '__'
    # else:
    #     file_name_no_ext, file_ext = os.path.splitext(file_name)
    #     guard_var = '__' + file_name_no_ext.replace('.','_') + '_' + suffix + file_ext.replace('.','_')  + '__'

    if suffix == '':
        guard_var = '__' + (prefix + '__')*bool(len(prefix)) + file_name.replace('.','_') + '__'
    else:
        file_name_no_ext, file_ext = os.path.splitext(file_name)
        guard_var = '__' + (prefix + '__')*bool(len(prefix)) + file_name_no_ext.replace('.','_') + '_' + suffix + file_ext.replace('.','_')  + '__'




    if uppercase:
        guard_var = guard_var.upper()


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
    for check_file_path, file_el in gb.file_dict.items():

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



# ====== isHeader ========

def isHeader(file_el):

    is_header = False

    file_name = file_el.get('name')
    if '.' in file_name:
        extension = file_name.rsplit('.',1)[1]
        if extension.lower() in ['hpp', 'h']:
            is_header = True

    return is_header

# ====== END: isHeader ========



# ====== getIncludeStatements ========

def getIncludeStatements(input_el, convert_loaded_to='none', add_extra_include_path=False, exclude_types=[],
                         input_element='class', forward_declared='exclude',
                         use_full_path=False, include_parents=False):

    include_statements = []

    # Check string arguments
    convert_loaded_to = convert_loaded_to.lower()
    input_element     = input_element.lower()
    forward_declared  = forward_declared.lower()
    if convert_loaded_to not in ['none', 'abstract', 'wrapper', 'wrapper_decl', 'wrapper_def']:
        raise Exception("getIncludeStatements: Keyword argument 'convert_loaded_to=' must be either 'none', 'abstract', 'wrapper', 'wrapper_decl' or 'wrapper_def'.")
    if input_element not in ['class', 'function']:
        raise Exception("getIncludeStatements: Keyword argument 'input_element=' must be either 'class' or 'function'.")
    if forward_declared not in ['include', 'exclude', 'only']:
        raise Exception("getIncludeStatements: Keyword argument 'forward_declared=' must be either 'include', 'exclude' or 'only'.")


    # Get list of all types used in this class/function (each entry is a dict)
    if input_element == 'class':
        all_types = getAllTypesInClass(input_el, include_parents=include_parents)
    elif input_element == 'function':
        all_types = getAllTypesInFunction(input_el)

    # Get file name and line number of the current class/function
    start_line_number = int( input_el.get('line') )
    start_file_el     = gb.id_dict[ input_el.get('file') ]
    start_file_path   = start_file_el.get('name')

    # Read file from beginning to position of class/function definition
    start_file         = open(start_file_path, 'r')
    start_file_content = start_file.readlines()[0:start_line_number]
    start_file_content = ''.join(start_file_content)
    start_file.close()

    # Identify included header files from this file (utils.identifyIncludedHeaders returns a dict of the form {header_file_name: xml_id})
    included_headers_dict = identifyIncludedHeaders(start_file_content, only_native=True)

    # Move up the header tree and identify all the relevant (native) included headers
    header_paths = [ gb.id_dict[file_id].get('name') for file_id in included_headers_dict.values() ]
    header_paths_done = []

    while len(header_paths) > 0:

        header_path = header_paths.pop()

        # Read header
        header         = open(header_path, 'r')
        header_content = header.read()
        header.close()

        # Identify new headers
        new_included_headers = identifyIncludedHeaders(header_content, only_native=True)

        # Add any new headers to included_headers_dict
        for file_name, file_id in new_included_headers.items():
            if file_name not in included_headers_dict.keys():
                included_headers_dict[file_name] = file_id

        # Add any new headers to the list of header files to check
        new_header_paths = [ gb.id_dict[header_id].get('name') for header_id in new_included_headers.values() ]
        for new_path in new_header_paths:
            if (new_path not in header_paths) and (new_path not in header_paths_done):
                header_paths.append(new_path)

        # Keep track of headers we've done
        header_paths_done.append(header_path)


    # Determine what include statements to generate:

    for type_dict in all_types:

        type_el   = type_dict['el']
        type_name = type_dict['class_name']

        if type_name in exclude_types:
            continue

        if isAcceptedType(type_el):

            if isFundamental(type_el):
                pass

            elif isEnumeration(type_el):
                pass

            elif isLoadedClass(type_el):

                # For each loaded class used in this class/function, check whether the corresponding class definition can be
                # found in the current file (above current class/function) or among the included headers. If no such class 
                # definition is found, it must be a case of simply using forward declaration.

                type_file_id = type_el.get('file')
                type_line_number = int(type_el.get('line'))

                if (type_file_id in included_headers_dict.values()) :
                    type_definition_found = True
                elif (type_file_id == input_el.get('file')) and (type_line_number < start_line_number):
                    type_definition_found = True
                else:
                    type_definition_found = False

                if (not type_definition_found) and (forward_declared=='exclude'):
                    # This must be a case of a type that is only forward declared. Don't include any header (as this will typically lead to a 'header loop').
                    continue

                elif (type_definition_found) and (forward_declared=='only'):
                    # This must be a case of a type that *is* fully declared, so we ignore it if only_forward:declared=True.
                    continue
                else:
                    if convert_loaded_to == 'none':

                        type_file_el = gb.id_dict[type_file_id]
                        type_file_full_path = type_file_el.get('name')

                        if isHeader(type_file_el):
                            type_file_basename = os.path.basename(type_file_full_path)
                            if add_extra_include_path:
                                include_statements.append('#include "' + os.path.join(cfg.add_path_to_includes, type_file_basename) + '"')
                            else:
                                include_statements.append('#include "' + type_file_basename + '"')
                        else:
                            print 'INFO: ' + 'Found declaration of loaded type "%s" in file "%s", but this file is not recognized as a header file. No header file include statement generated.' % (type_name['long_templ'], type_file_full_path)

                    else:
                        if use_full_path:
                            header_key = convert_loaded_to + '_fullpath'
                        else:
                            header_key = convert_loaded_to

                        if add_extra_include_path:
                            include_statements.append('#include "' + os.path.join(cfg.add_path_to_includes, gb.new_header_files[type_name['long']][header_key]) + '"')
                        else:
                            include_statements.append('#include "' + gb.new_header_files[type_name['long']][header_key] + '"')

            elif isStdType(type_el):

                if type_name['long'] in cfg.known_class_headers:
                    header_name = cfg.known_class_headers[type_name['long']]
                    if (header_name[0] == '<') and (header_name[-1] == '>'):
                        include_statements.append('#include ' + cfg.known_class_headers[type_name['long']])
                    else:
                        include_statements.append('#include "' + cfg.known_class_headers[type_name['long']] + '"')
                else:
                    # warnings.warn("The standard type '%s' has no specified header file. Please update modules/cfg.py. No header file included." % type_name['long_templ'])
                    print 'INFO: ' + 'The standard type "%s" has no specified header file. Please update modules/cfg.py. No header file include statement generated.' % type_name['long_templ']

            else:
                if type_name['long'] in cfg.known_class_headers:
                    include_statements.append('#include "' + cfg.known_class_headers[type_name['long']] + '"')
                else:
                    # warnings.warn("The type '%s' has no specified header file. Please update modules/cfg.py. No header file included." % type_name['long'])
                    print 'INFO: ' + 'The type "%s" has no specified header file. Please update modules/cfg.py. No header file include statement generated.' % type_name['long_templ']
        else:
            # warnings.warn("The type '%s' is unknown. No header file included." % type_name['long'])
            print 'INFO: ' + 'The type "%s" is unknown or not loadable. No header file include statement generated.' % type_name['long_templ']

    # Remove duplicates and return list
    include_statements = list( OrderedDict.fromkeys(include_statements) )

    return include_statements

# ====== END: getIncludeStatements ========



# ====== getOriginalHeaderPath ========

def getOriginalHeaderPath(el, full_path=False):

    file_id = el.get('file')
    file_el = gb.id_dict[file_id]

    file_full_path = file_el.get('name')

    if isHeader(file_el):

        if full_path:
            header_path = file_full_path
        else:
            header_path = os.path.basename(file_full_path)

    else:
        raise exceptions.ReturnError("The file %s is not recognized as a header file." % file_full_path)

    return header_path

# ====== END: getOriginalHeaderPath ========



# ====== constrNamespaceFromTags ========

def constrNamespaceFromTags(content, new_namespace, open_tag, close_tag):

    new_namespace_list = new_namespace.split('::')

    new_content = content
    while True:

        prev_content = new_content

        # Find tag positions
        open_pos  = prev_content.find(open_tag)
        close_pos = prev_content.find(close_tag)

        if (open_pos == -1) and (close_pos == -1):
            # No tags found. Do nothing more.
            # return content
            break
        elif (open_pos == -1) or (close_pos == -1):
            raise Exception('Matching pair of namespace tags %s and %s not found in given content.' % (open_tag, close_tag))
        else:
            pass

        # Split content into three parts
        content_before = prev_content[:open_pos]
        content_within = prev_content[open_pos:close_pos]
        content_after  = prev_content[close_pos:]

        # Remove the namespace tags
        content_within = content_within.replace(open_tag, '')
        content_after  = content_after.replace(close_tag, '', 1)

        # Add indentation to middle part
        content_within = addIndentation(content_within, cfg.indent*len(new_namespace_list))

        # Contruct new namespace and combine code
        open_new_namespace_code  = constrNamespace(new_namespace_list, 'open', indent=cfg.indent)
        close_new_namespace_code = constrNamespace(new_namespace_list, 'close', indent=cfg.indent)

        new_content = content_before + open_new_namespace_code + content_within + close_new_namespace_code + content_after

    return new_content

# ====== END: constrNamespaceFromTags ========



# ====== replaceCodeTags ========

def replaceCodeTags(input, file_input=False, write_file=False):

    # Input is either a file name or a string with content
    if file_input:
        f = open(input, 'r')
        new_content = f.read()
        f.close()
    else:
        new_content = input

    # Replace various tags in template code with code specific for the current backend
    new_content = new_content.replace('__BACKEND_NAME__'         ,  cfg.gambit_backend_name)
    new_content = new_content.replace('__BACKEND_VERSION__'      ,  cfg.gambit_backend_version)
    new_content = new_content.replace('__BACKEND_SAFE_VERSION__' ,  gb.gambit_backend_safeversion)
    new_content = new_content.replace('__CODE_SUFFIX__'          ,  gb.code_suffix)
    
    new_content = new_content.replace('__PATH_TO_FRWD_DECLS_ABS_CLASSES_HEADER__', os.path.join(gb.gambit_backend_types_basedir, gb.gambit_backend_name_full, gb.frwd_decls_abs_fname + cfg.header_extension))
    new_content = new_content.replace('__PATH_TO_IDENTIFICATION_HEADER__'        , os.path.join(gb.gambit_backend_types_basedir, gb.gambit_backend_name_full, 'identification.hpp'))
    new_content = new_content.replace('__PATH_TO_BACKEND_UNDEFS_HEADER__'        , os.path.join(gb.gambit_backend_incl_dir, "backend_undefs.hpp"))


    # Should a new file be written?
    if file_input and write_file:
        f.open(input,'w')
        f.write(new_content)
        f.close()

    # Return code
    return new_content

# ====== END: replaceCodeTags ========



# ====== removeCodeTags ========

def removeCodeTags(content, remove_tags_list):

    new_content = content
    # Remove tags
    for tag in remove_tags_list:
        new_content = new_content.replace(tag, '')

    # Return code
    return new_content

# ====== END: removeCodeTags ========



# ====== constrLoadedTypesHeaderContent ======

def constrLoadedTypesHeaderContent():

    # 
    # Construct the code lines for the loaded classes, containg all the factory symbols and argument brackets for that class
    #
    class_lines = []

    for class_name in gb.classes_done:

        if not class_name['long'] in gb.factory_info.keys():
            print "WARNING: No factory function has been generated for class '%s'." % class_name['long']

        else:

            class_line = '  (( /*class*/'

            namespace, class_name_short = removeNamespace(class_name['long'], return_namespace=True)

            if namespace == '':
                namespace_list = []
            else:
                namespace_list = namespace.split('::')

            for ns_part in namespace_list:
                class_line += '(' + ns_part + ')'

            class_line += '(' + class_name['short'] + '),'


            class_line += '    /*constructors*/'

            # for factory_name, args_bracket, symbol in factory_dict[class_name['long']]:
            #     class_line += '(("' + symbol + '",' + args_bracket + ')) '

            for info_dict in gb.factory_info[ class_name['long'] ]:
                class_line += '(("' + info_dict['symbol'] + '",' + info_dict['args_bracket'] + ')) '

            class_line += ')) \\'
            class_lines.append(class_line)

    class_lines_code  = ''
    class_lines_code += '#define ' + gb.gambit_backend_name_full + '_all_data \\\n'
    class_lines_code += '\n'.join(class_lines) + '\n'


    #
    # Construct include guards with additional  ' 1' appended to the line starting with #define 
    #
    # incl_guard = addIncludeGuard('', 'loaded_types.hpp', extra_string=gb.gambit_backend_name_full)
    incl_guard = addIncludeGuard('', 'loaded_types.hpp', prefix='', suffix=gb.gambit_backend_name_full)
    incl_guard_lines = incl_guard.split('\n')

    incl_guard_start = '\n'.join(incl_guard_lines[:2]) + ' 1\n'
    incl_guard_end   = incl_guard_lines[-2] + '\n'


    #
    # Construct include statements
    #
    incl_statements_code = ''
    for class_name in gb.classes_done:
        if class_name['long'] in gb.factory_info.keys():
            namespace, class_name_short = removeNamespace(class_name['long'], return_namespace=True)
            incl_statements_code += '#include "' + cfg.wrapper_header_prefix + class_name['short'] + cfg.header_extension + '"\n'
    incl_statements_code += '#include "identification.hpp"\n'
    # incl_statements_code += '#include "' + os.path.join(gb.gambit_backend_types_basedir, gb.gambit_backend_name_full, 'identification.hpp') + '"\n'


    #
    # Combine everything to construct header code
    #
    code  = ''
    code += incl_guard_start

    code += '\n'
    code += incl_statements_code

    code += '\n'
    code += '// Indicate which types are provided by this backend, and what the symbols of their factories are.\n'
    code += class_lines_code

    code += '\n'
    code += '// If the default version has been loaded, set it as default.\n'
    code += '#if ALREADY_LOADED(CAT_3(BACKENDNAME,_,CAT(Default_,BACKENDNAME)))\n'
    code += '  SET_DEFAULT_VERSION_FOR_LOADING_TYPES(BACKENDNAME,SAFE_VERSION,CAT(Default_,BACKENDNAME))\n'
    code += '#endif\n'

    code += '\n'
    code += '// Undefine macros to avoid conflict with other backends.\n'
    code += '#include "' + os.path.join(gb.gambit_backend_incl_dir, "backend_undefs.hpp") + '"\n'

    code += '\n'
    code += incl_guard_end

    return code
    
# ====== END: constrLoadedTypesHeaderContent ======



# ====== constrEnumDeclHeader ========

def constrEnumDeclHeader(enum_el_list, file_output_path):

    import modules.classutils as classutils

    # If this is the first time this function is executed, read initial code from header_templates/ folder
    if file_output_path not in gb.new_code.keys():
        f = open('header_templates/standard_header_template.hpp')
        initial_code = f.read()
        f.close()
        initial_code_tuple = (0, initial_code)
        gb.new_code[file_output_path] = {'code_tuples':[initial_code_tuple], 'add_include_guard':True}


    current_code = gb.new_code[file_output_path]['code_tuples'][0][1]

    insert_code = ''
    tag_pos = current_code.find('__INSERT_CODE_HERE__')

    for enum_el in enum_el_list:

        # Skip any enumerated type that is not native to the external code
        if not isNative(enum_el):
            continue

        # Get enum names and values
        enum_members_list = []
        for sub_el in enum_el.getchildren():
            if sub_el.tag == 'EnumValue':
                member_string = sub_el.get('name') + '=' + sub_el.get('init')
                enum_members_list.append(member_string)
        all_members_string = ', '.join(enum_members_list)

        # Get namespace list
        namespace_list = getNamespaces(enum_el)
        
        n_indents = len(namespace_list)

        # - Open namespace
        insert_code += constrNamespace(namespace_list, 'open')

        insert_code += ' '*n_indents*cfg.indent + 'typedef enum {' + all_members_string + '} ' + enum_el.get('name') + ';\n'

        # - Close namespace
        insert_code += constrNamespace(namespace_list, 'close')

    new_code = current_code[:tag_pos] + insert_code + current_code[tag_pos:]
   

    # Replace other code tags
    new_code = replaceCodeTags(new_code)


    new_code_tuple = (0,new_code)

    # Overwrite existing code tuple
    gb.new_code[file_output_path]['code_tuples'] = [(new_code_tuple)]

# ====== END: constrEnumDeclHeader ========



# ====== gccxmlRunner ========

# Calls gccxml from the shell (via modules.shelltimeout).

def gccxmlRunner(source_path, include_paths_list, xml_output_path, timeout_limit=30., poll_interval=0.5):

    # Construct gccxml command to run
    gccxml_cmd = 'gccxml '

    # - Add include paths
    for incl_path in include_paths_list:
        gccxml_cmd += '-I' + incl_path + ' '

    # - Add the factory source file (full path)
    gccxml_cmd += source_path + ' '

    # - Add gccxml option that specifies the xml output file: input_file_short_name.xml
    gccxml_cmd += '-fxml=' + xml_output_path

    # Run gccxml
    print 'Runing command: ' + gccxml_cmd
    proc, output, timed_out = shelltimeout.shrun(gccxml_cmd, timeout_limit, use_exec=True, poll_interval=poll_interval)

    did_fail = False

    # Check for timeout or error
    if timed_out:
        print 'ERROR: gccxml timed out.'
        did_fail = True
    elif proc.returncode != 0:        
        print 'ERROR: gccxml failed.'
        did_fail = True

    # Print error report
    if did_fail:
        print
        print '==== START GCCXML OUTPUT ===='
        print
        print output
        print '==== END GCCXML OUTPUT ===='
        print
        raise Exception('gccxml failed')
    
    else:
        print 'Command finished successfully.'
    print

# ====== END: gccxmlRunner ========


