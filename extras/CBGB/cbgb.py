#!/usr/bin/python

#
#   =========================================================
#   |                                                       |
#   |   CBGB - Common Block harvester for GAMBIT Backends   |
#   |                                                       |
#   =========================================================
#
#   Author: Anders Kvellestad (anders.kvellestad@fys.uio.no)
#

from modules import utils
from modules import gb
from modules import cfg
import sys
from collections import OrderedDict


print
print
print '  ~~~~  CBGB  -  Common Block harvester for GAMBIT Backends  ~~~~'
print
print

print '  Input source file:'
print '  ------------------'
print
print '    %s' % cfg.src_file_path
print 
print
print '  Requested common blocks:'
print '  ------------------------'
print

cb_listing = '    '
for i, cb_name in enumerate(cfg.load_common_blocks,1):
    cb_listing += cb_name + ', '

    if i%4==0:
        cb_listing += '\n    '
print cb_listing.rstrip().rstrip(',')
print
print


# Read source file:
src_file = open(cfg.src_file_path, 'r')
src_content = src_file.read()
src_file.close()

# Open output files
out_file_be_types = open(gb.output_file_path_be_types, 'w')
out_file_frontent = open(gb.output_file_path_frontent, 'w')


#
# Do some reformatting of the source code:
#
# - Convert all code to lower-case. (Fortran is not case-sensitive.)
# - Replace tabs with spaces.
# - Split source code into a list of lines.
# - Remove all comments and blank lines.
# - Remove all statement labels.
# - Combine continued source lines into single lines.
# - Remove leading and trailing blanks,
#

# Convert all source code to lower-case. (Fortran is not case-sensitive.)
src_content = src_content.lower()

# Convert tabs to spaces.
src_content = src_content.replace('\t', ' '*cfg.tabs_to_n_spaces)

# Split source code into a list of code lines.
src_lines = src_content.splitlines()

# Remove comments.
src_lines = utils.removeComments(src_lines)

# Remove statement labels.
src_lines = utils.removeStatementLabels(src_lines)

# Remove blank lines.
src_lines = utils.removeBlankLines(src_lines)

# Join continued lines.
src_lines = utils.joinContinuedLines(src_lines)

# Remove leading and trailing blanks
src_lines = utils.removeLeadingTrailingBlanks(src_lines)



# Identify the various parts of the code: program, functions and subroutines.
# Return a dict with the following structure: 
#  {
#   'some_subroutine_name' : { 'category'  : 'subroutine', 
#                              'code_lines': [line1, line2, ...] }, 
#   'some_function_name'   : { ... }, 
#    ...
#  }
code_parts_dict = utils.getCodeParts(src_lines)


# Remove potential duplicates in list of common blocks to load.
cfg.load_common_blocks = list( OrderedDict.fromkeys(cfg.load_common_blocks) )

# Create a copy of common block list and convert all common block names to lower-case.
common_blocks_left = [cb_name.lower() for cb_name in cfg.load_common_blocks]


#
# Loop over code parts to extract info on common blocks.
#

print '  Searching for common blocks:'
print '  ----------------------------'
print
for code_part_name, code_dict in code_parts_dict.items():

    code_lines    = code_dict['code_lines']
    code_category = code_dict['category']

    # Get list of dicts with info on all common blocks in this code part.
    cb_dicts = utils.getCommonBlockDicts(code_lines)

    # Loop over the common blocks found.
    for cb_dict in cb_dicts:

        if cb_dict['name'] in common_blocks_left:

            print "    In %s '%s': Found common block: '%s'" % (code_category, code_part_name, cb_dict['name'])

            # Get dict of dicts with info on the member variables for this common block.
            var_info_dict = utils.getVariablesDict(code_lines, cb_dict['member_names'])

            # Generate code for the backend types header.
            out_file_be_types.write('\n')            
            out_file_be_types.write( utils.generateGambitCommonBlockDecl(cb_dict, var_info_dict) )

            # Generate code for the frontend header.
            out_file_frontent.write('\n')            
            out_file_frontent.write( utils.generateGambitFrontendCode(cb_dict) )

            # Remove common block from list of blocks remaining.
            common_blocks_left.remove(cb_dict['name'])

    #
    # END: loop over common blocks in this code part
    #

    # Break out if all reqested blocks are done.
    if len(common_blocks_left) == 0:
        break

#
# END: loop over code parts
#

print
print
print '  Summary:'
print '  --------'

# Check that all requested common blocks where found.
if len(common_blocks_left) > 0:
    print 
    for cb_name in common_blocks_left:
        print "    WARNING: Reqested common block '%s' was not found." % cb_name 


# Close output files
out_file_be_types.close()
out_file_frontent.close()
print 
print '    Generated code for GAMBIT written to files: %s, %s.' % (gb.output_file_path_be_types, gb.output_file_path_frontent)


print
print
print '  ~~~~  Done!  ~~~~'
print
print 

