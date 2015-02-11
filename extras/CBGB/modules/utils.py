#
#   ==================================
#   |                                |
#   |   Utility functions for CBGB   |
#   |                                |
#   ==================================
#

from collections import OrderedDict
from modules import gb
from modules import cfg


# ====== removeComments ========

# Takes a list of code lines and removes comments.
# Lines starting with 'c' are replaced by an empty line, 
# while for lines containing '!' everything after '!'
# is removed.

def removeComments(code_lines):

    code_lines_nocomment = []

    for line in code_lines:

        if len(line) == 0:
            code_lines_nocomment.append('')
            continue

        if (cfg.format == 'fixed') and (line[0] == 'c'):
            new_line = ''

        elif '!' in line:
            pos = line.find('!')
            new_line = line[:pos]

        else:
            new_line = line

        code_lines_nocomment.append(new_line)

    return code_lines_nocomment

# ====== END: removeComments ========



# ====== removeBlankLines ========

# Removes any empty (all whitespace) strings from a list of strings.

def removeBlankLines(code_lines):

    # Walk through the list of code lines backwards and discard 
    # any lines that contain nothing but whitespace.
    for i in range(len(code_lines))[::-1]:
        if code_lines[i].strip() == '':
            code_lines.pop(i)

    return code_lines

# ====== END: removeBlankLines ========



# ====== removeLeadingTrailingBlanks ========

# Removes leading and trailing blanks from the strings
# in a list of strings.

def removeLeadingTrailingBlanks(code_lines):

    for i in range(len(code_lines)):
        code_lines[i] = code_lines[i].lstrip().rstrip()

    return code_lines

# ====== END: removeLeadingTrailingBlanks ========



# ====== removeStatementLabels ========

# Replaces statement labels with empty spaces.
# (A statement label is a number given as the first 
# non-blank part of a statement.)

def removeStatementLabels(code_lines):

    for i in range(len(code_lines)):

        line = code_lines[i]

        if cfg.format == 'fixed':
            label = line[0:5].strip()
            if label.isdigit():
                code_lines[i] = line.replace(label, ' '*len(label), 1)

        elif cfg.format == 'free':

            line_list = line.split()
            if (len(line_list) > 0):
                label = line_list[0]
                if label.isdigit():
                    code_lines[i] = line.replace(label, ' '*len(label), 1)                    

        else:
            raise Exception("cfg.format must be set to either 'fixed' or 'free'.")

    return code_lines

# ====== END: removeStatementLabels ========



# ====== joinContinuedLines ========

def joinContinuedLines(code_lines):

    joined_code_lines = ['']

    if cfg.format == 'fixed':

        for line in code_lines:
            
            # Check for line continuation (any character at column 6).
            # (This assumes that len(line) >= 6 for all lines in code_lines,
            # which should be OK due to prior code formatting.)

            # - If found, append to previous line.
            try:
                if line[5] not in [' ','\t']:
                    joined_code_lines[-1] += line[6:]
            except:
                print [line]
                raise
            
            # - If not found, store current_line and start constructing a new.
            else:
                joined_code_lines.append(line)


    elif cfg.format == 'free':

        continue_line = False
        for line in code_lines:

            if continue_line:
                if line.lstrip()[0] == '&':
                    joined_code_lines[-1] += line.lstrip()[1:].rstrip().rstrip('&')
                else:
                    joined_code_lines[-1] += line.rstrip().rstrip('&')
            else:
                joined_code_lines.append(line.rstrip().rstrip('&'))


            # Check for line continuation. (Line ends with '&'.)
            if line.rstrip()[-1] == '&':
                continue_line = True
            else:
                continue_line = False

    else:
        raise Exception("cfg.format must be set to either 'fixed' or 'free'.")

    if joined_code_lines[0] == '':
        joined_code_lines.pop(0)

    return joined_code_lines


# ====== END: joinContinuedLines ========




# ====== getCodeParts ========

def getCodeParts(code_lines):

    code_parts_dict = OrderedDict()

    unnamed_part_counter = 1
    start_line  = 0
    end_line    = 0
    found_start = False
    for i, line in enumerate(code_lines):

        #
        # Detect start of a code part
        #
        at_start_line = False

        if not found_start:
            if ('subroutine ' in line) or ('function ' in line) or ('program ' in line) or (line[-7:] == 'program'):
                found_start = True
                at_start_line = True

        if at_start_line:
            
            start_line = i
            
            if 'subroutine' in line:
                keyword = 'subroutine'
            elif 'function' in line:
                keyword = 'function'
            elif 'program' in line:
                keyword = 'program'
            else:
                print 'WARNING: Neither of "subroutine", "function" nor "program" was found in line %i: %s' % (i+1,line)


            line_list = line.split()
            keyword_index = line_list.index(keyword)

            if len(line_list) == keyword_index+1:
                name = 'unnamed_' + keyword + '_' + str(unnamed_part_counter)
                unnamed_part_counter += 1

            else:
                name_item = line_list[line_list.index(keyword)+1]
                if '(' in name_item:
                    name = name_item[:name_item.find('(')]
                else:
                    name = name_item
        

        #
        # Detect end of a code part
        #
        elif line.replace(' ','').strip() in ['end','end'+keyword, 'end'+keyword+name]:

            if not found_start:
                print 'WARNING: Unbalanced "end" statement at line %i' % (i+1)
                found_start = False
                continue

            end_line = i

            # Store in dict
            code_parts_dict[name] = { 
                                      'category'   : keyword,
                                      'code_lines' : code_lines[start_line:end_line+1] 
                                    }
            # Reset flag
            found_start = False

    return code_parts_dict

# ====== END: getCodeParts ========




# ====== getImplicitDefs ========

# Return a dict with the following structure: 
#   { 
#     'a': ('double precision',1), 
#     'b': ('real',8),
#     'c': (None,None),
#     ...
#   }
#

def getImplicitDefs(code_lines):

    implicit_defs = gb.default_implicit_types

    for i,line in enumerate(code_lines):

        # Split line into words
        line_list = line.split()

        # Look for 'implicit' statement
        if line_list[0] == 'implicit':

            # If 'implicit none', then no other 'implicit' statements are allowed
            if line_list[1] == 'none':
                return dict.fromkeys(gb.alphabet,(None,None))

            # Remove the 'implicit' keyword
            typedef_line   = ' '.join(line_list[1:])

            # If there are multiple implicit statements on a single line,
            # split them up and treat them separately.
            for temp_line in typedef_line.split(')'):

                # Do a bunch of string manipulations to identify
                # the type name (e.g. 'double precision') and 
                # character specifications (e.g. 'a-z').
                if temp_line == '':
                    continue

                temp_line = temp_line.replace('(','')
                temp_line = temp_line.replace(',',' ')
                temp_line = temp_line.strip()
                while ' -' in temp_line:
                    temp_line = temp_line.replace(' -','-')
                while '- ' in temp_line:
                    temp_line = temp_line.replace('- ','-')
                temp_line = ' '.join(temp_line.split())
                temp_line_list = temp_line.split()

                char_list = []
                type_name_list = []
                for entry in temp_line_list:
                    if ((len(entry)==1) and (entry in gb.alphabet)) or (len(entry)==3 and (entry[1]=='-')):
                        char_list.append(entry)
                    else:
                        type_name_list.append(entry)

                full_type_name = ''.join(type_name_list)
                if '*' in full_type_name:
                    type_name, type_size_str = full_type_name.split('*')
                    type_size = int(type_size_str)
                else:
                    type_name = full_type_name
                    type_size = 1


                # Loop through the character specifiers in char_list
                # and set the correct types in the implicit_defs dict

                for char in char_list:

                    if (len(char)==1) and (char in gb.alphabet):
                        implicit_defs[char] = (type_name,type_size)

                    elif (len(char)==3 ) and (char[1]=='-'):
                        start_char = char[0]
                        end_char = char[2]

                        for key_char in implicit_defs.keys():
                            if (key_char >= start_char) and (key_char <= end_char):
                                implicit_defs[key_char] = (type_name,type_size)

    return implicit_defs

# ====== END: getImplicitDefs ========




# ====== getCommonBlockDicts ========

def getCommonBlockDicts(code_lines):

    cb_dicts = []

    for line in code_lines:

        # Remove whitespaces
        line = line.replace(' ','')

        # Ignore lines that don't start with 'common/'
        if (len(line) < 7) or (line[:7] != 'common/'):
            continue

        # Identify common block name and names of member variables
        line_list   = line.split('/')
        cb_name     = line_list[1]
        var_seq_str = line_list[2] 

        var_dicts = parseVariableSequence(var_seq_str)
        var_names = var_dicts.keys()

        cb_dicts.append( {'name':cb_name, 'member_names':var_names} )

    return cb_dicts

# ====== END: getCommonBlockDicts ========




# ====== isVariableDecl ========

def isVariableDecl(line_in, return_type=False):

    is_variable_decl = False
    type_name = ''
    type_size = 1

    line = line_in
    line = line.replace(',',' ').replace('*',' * ')
    line = line.replace('(', ' (').replace(')',') ')
    line = ' '.join(line.split())

    line_list = line.split()
    for i in [3,2,1]:
        check_type = ''.join(line_list[:i])

        # Check that we can deal with this Fortran type.
        if check_type in gb.type_translation_dict.keys():

            # If type is 'character*', identify the integer that specifies the 
            # string length.
            if check_type=='character':
                if (line_list[1] == '*') and (line_list[2].isdigit()):
                    check_type += '*' + line_list[2]

            if '*' in check_type:
                type_name, type_size_str = check_type.split('*')
                type_size = int(type_size_str)
            else:
                type_name = check_type

            is_variable_decl = True
            break

    if return_type:
        return is_variable_decl, type_name, type_size
    else:
        return is_variable_decl

# ====== END: isVariableDecl ========




# ====== isDimensionStatement ========

def isDimensionStatement(line_in):

    is_dim_stmnt = False

    line = line_in
    line_list = line.split()

    if (len(line_list) > 1) and (line_list[0] == 'dimension'):
        is_dim_stmnt = True

    return is_dim_stmnt

# ====== END: isDimensionStatement ========



# ====== getArrayIndicesTuples ========

# Example:
# Input:  '-2:10,7,1:2'
# Output: [(-2,7), (1,7), (1,2)]

def getArrayIndicesTuples(dimensions_str):

    indicies_tuples = []

    # Check for empty dimensions string
    if dimensions_str == '':
        return indicies_tuples

    # Loop over comma-separated entries in dimensions_str
    for dim_str in dimensions_str.split(','):

        if ':' in dim_str:
            start_index, end_index = [int(s) for s in dim_str.split(':')]

        else:
            start_index = 1
            end_index   = int(dim_str)

        indicies_tuples.append( (start_index,end_index) )

    return indicies_tuples

# ====== END: getArrayIndicesTuples ========



# ====== getVariablesDict ========

def getVariablesDict(code_lines, get_variables):

    return_var_dicts = OrderedDict.fromkeys(get_variables, value=None)

    implicit_defs = getImplicitDefs(code_lines)


    for line in code_lines:

        #
        # First, make use of all variable type declaration lines
        #
        
        is_var_decl, type_name, type_size = isVariableDecl(line, return_type=True)

        if is_var_decl:

            # Remove type name from beginning of line so that
            # only the list of variable names remain.
            full_type_name = type_name + '*' + str(type_size)

            line_list = line.split()
            i = 1
            while i <= len(line_list):
                if ''.join(line_list[:i]) in full_type_name:
                    i += 1
                    continue
                else:
                    break
            var_seq = ''.join(line_list[i-1:])

            # Parse line to extract info on the different variables
            var_dicts = parseVariableSequence(var_seq)

            # Append type_name and type_size to var_dicts
            for var_name in var_dicts.keys():

                # - Add type name
                var_dicts[var_name]['type'] = type_name
                # - Use the maximum of the sizes specified in the type name and in the variable sequence
                #   (Normally one of these should be 1 by default.)
                var_dicts[var_name]['size'] = max(type_size,var_dicts[var_name]['size'])

                # Check for character array type:
                if (var_dicts[var_name]['type'] == 'character'): 
                    dim_str = var_dicts[var_name]['dimension']
                    size    = var_dicts[var_name]['size']
                    if (dim_str == '') and (size > 1):
                        var_dicts[var_name]['dimension'] = '1:%i' % size

            # For requested variables, append the variable dicts to return_var_dicts
            for var_name in var_dicts.keys():
                if var_name in get_variables:
                    return_var_dicts[var_name] = var_dicts[var_name]



        #
        # Then, check all the 'dimension' statements
        #

        is_dim_stmnt = isDimensionStatement(line)

        if is_dim_stmnt:

            # Remove whitespace and 'dimension' keyword
            line = line.replace(' ','')
            line = line.replace('dimension','',1)

            # Parse line to extract info on the different variables
            dim_var_dicts = parseVariableSequence(line)

            # For variables that already exist in return_var_dicts, simply
            # update the 'dimension'. For variables that don't exist in 
            # return_var_dicts, create a new entry based on implicit types.
            for var_name in dim_var_dicts.keys():

                if var_name in get_variables:

                    # If info on this variable has not yet been added to return_var_dicts,
                    # insert a complete dict
                    if return_var_dicts[var_name] == None:
                        # Get type from implicit types
                        first_char = var_name[0]
                        type_name, type_size = implicit_defs[first_char]

                        if type_name == None or type_size == None:
                            raise Exception("No type declaration (neither explicit nor implicit) was found for variable '%s'." % var_name)

                        return_var_dicts[var_name] = { 
                                                       'type'     : type_name,
                                                       'dimension': dim_var_dicts[var_name]['dimension'],
                                                       'size'     : type_size
                                                     }

                    # If info on this variable already exists, simply update the 'dimension' entry in the
                    # correct dict
                    else:
                        return_var_dicts[var_name]['dimension'] = dim_var_dicts[var_name]['dimension']


    #
    # END: Loop over code lines
    #

    #
    # Finally, add any missing variables that have not appeared in explicit type
    # declarations or 'dimension' statements
    #

    for get_var_name in get_variables:

        if return_var_dicts[get_var_name] == None:

            # Get type from implicit types
            first_char = get_var_name[0]
            type_name, type_size = implicit_defs[first_char]

            if type_name == None or type_size == None:
                raise Exception("No type declaration (neither explicit nor implicit) was found for variable '%s'." % get_var_name)

            return_var_dicts[get_var_name] = { 
                                              'type'     : type_name,
                                              'dimension': '',
                                              'size'     : type_size
                                             }


    return return_var_dicts

# ====== END: getVariablesDict ========




# ====== parseVariableSequence ========

# Input : "var1*100, var2(1:20)*20, var3"
#
# Output: { 
#           'var1': { 'size': 100, 'dimension': ''       }, 
#           'var2': { 'size': 20,  'dimension': '(1:20)' },
#           'var3': { 'size': 1,   'dimension': ''       }
#         }

def parseVariableSequence(var_seq_str):

    result_dict = OrderedDict()

    line = var_seq_str

    # Remove all whitespace
    line = line.replace(' ','')

    # Split into separate variables by detecting commas 
    # (excluding commas inside brackets).
    i = 0
    bracket_balance = 0
    while i < len(line):
        
        char = line[i]

        # Keep track of the brackets
        if char == '(':
            bracket_balance += 1
        elif char == ')':
            bracket_balance -= 1

        # If a comma is found, replace it with a whitespace
        if (char == ',') and (bracket_balance == 0):
            line = line[:i] + ' ' + line[i+1:]

        # Increment index
        i += 1

    # Split line at whitespaces
    var_str_list = line.split()

    for var_str in var_str_list:

        # Check for dimension bracket and size integer
        has_dim_bracket = bool('(' in var_str and ')' in var_str)
        has_size_int    = bool('*' in var_str)

        # Insert whitespace to separate variable name, dimension bracket and size integer
        var_str = var_str.replace('(',' ').replace(')',' ').replace('*',' ')

        # Split at whitespace
        var_str_list = var_str.split()

        # Identify name, dimension, size
        if has_dim_bracket and has_size_int:
            var_name     = var_str_list[0]
            var_dim_str  = var_str_list[1]
            var_size     = int(var_str_list[2])
        elif has_dim_bracket and not has_size_int:
            var_name     = var_str_list[0]
            var_dim_str  = var_str_list[1]
            var_size     = 1
        elif has_size_int and not has_dim_bracket:
            var_name     = var_str_list[0]
            var_dim_str  = ''
            var_size     = int(var_str_list[1])
        else:
            var_name     = var_str_list[0]
            var_dim_str  = '' 
            var_size     = 1

        # Append to result_dict
        result_dict[var_name] = {'dimension': var_dim_str, 'size': var_size}

    return result_dict

# ====== END: parseVariableSequence ========



# ====== generateGambitCommonBlockDecl ========

def generateGambitCommonBlockDecl(cb_dict, var_info_dict):

    indent = ' '*4

    code = ''

    cb_name = cb_dict['name']
    cb_type_name = cb_name + '_type'

    code += 'struct %s\n' % cb_type_name
    code += '{\n'

    for var_name, var_dict in var_info_dict.items():

        c_type_name = getCTypeName(var_dict)

        code += indent + c_type_name + ' ' + var_name + ';\n'

        # array_indices_tuples = getArrayIndicesTuples(var_dict['dimension'])

        # is_array = bool(len(array_indices_tuples))
        
        # if is_array:
    
        #     if c_type_name == 'char':
        #         code += indent + '// Fstring, FstringArray: REMAINS TO BE IMPLEMENTED!\n'
    
        #     else:
        #         all_indices_list = [i for tpl in array_indices_tuples for i in tpl]
        #         all_indices_str = ','.join( map(str,all_indices_list) )
        #         template_bracket = '< %s,%s >' % (c_type_name, all_indices_str)
        #         code += indent + 'Farray' + template_bracket + ' ' + var_name + ';\n'
        # else:
        #     code += indent + c_type_name + ' ' + var_name + ';\n'

    code += '};\n'

    return code

# ====== END: generateGambitCommonBlockDecl ========



# ====== generateGambitFrontendCode ========

def generateGambitFrontendCode(cb_dict):

    code = ''

    cb_name            = cb_dict['name']
    cb_type_name       = cb_name + '_type'
    cb_capability_name = cfg.capability_prefix + cb_name + cfg.capability_suffix
    cb_mangled_symbol  = cb_name + '_'

    code += 'BE_VARIABLE(%s, %s, "%s", "%s")\n' % (cb_type_name, cb_name, cb_mangled_symbol, cb_capability_name)

    return code

# ====== END: generateGambitFrontendCode ========



# ====== getCTypeName ========

def getCTypeName(var_dict):

    fortran_type_name = var_dict['type']
    c_type_base_name = gb.type_translation_dict[fortran_type_name]

    array_indices_tuples = getArrayIndicesTuples(var_dict['dimension'])

    # Is this variable an array?
    if (fortran_type_name != 'character') and (len(array_indices_tuples) > 0):
        is_array = True
    elif (fortran_type_name == 'character') and (len(array_indices_tuples) > 1):
        is_array = True
    else:
        is_array = False

    # For arrays, construct a string of comma-separated array indices
    if is_array:
        all_indices_list = [i for tpl in array_indices_tuples for i in tpl]
        all_indices_str = ','.join( map(str,all_indices_list) )

    #
    # Determine the correct C++ type name
    #

    # Special treatment for the character type
    if (fortran_type_name == 'character') and (var_dict['size'] > 1):
        if is_array:
            template_bracket = '< %i,%s >' % (var_dict['size'], all_indices_str)
            c_type_name = 'FstringArray' + template_bracket
        else:
            c_type_name = 'Fstring<%i>' % var_dict['size'] 

    # All other types
    else:
        if is_array:
            template_bracket = '< %s,%s >' % (c_type_base_name, all_indices_str)
            c_type_name = 'Farray' + template_bracket
        else:
            c_type_name = c_type_base_name

    # Return result
    return c_type_name

# ====== END: getCTypeName ========



# ====== addNamespace ========

# Encapsulate code string in a namespace

def addNamespace(code, namespace_name, indent=4):

    # Add indentation
    code_lines = [' '*indent + line for line in code.splitlines()]
    code = '\n'.join(code_lines)

    # Add namespace
    code = 'namespace ' + namespace_name + '\n' + '{\n' + code + '\n}\n'

    return code

# ====== END: addNamespace ========


