#include "common.h"

typedef char FileArg[_ARGUMENT_LENGTH_MAX_];

/* after reading a given file, all relevant information stored in this structure, in view of being processed later*/
struct file_content {
  char * filename;
  int size;
  FileArg * name;  /**< list of (size) names */
  FileArg * value; /**< list of (size) values */
  short * read;    /**< set to _TRUE_ if this parameter is effectively read */
};
