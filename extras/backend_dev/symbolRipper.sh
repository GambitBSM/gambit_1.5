#!/bin/sh
# use nm to get the symbol list, use awk to grab -only- the names
nm -g -P --defined-only libfirst.so | awk '
{print $1}
' > libfirst.sym
# do that again, but use c++filt to unmangle the names
#   (maybe unnecessary for linux, but necessary on mac)
nm -g -P --defined-only libfirst.so | awk '
{print $1}
' | c++filt > libfirst.unmangled
# compare those files with diff and use awk to create a file which
#   acts as a dictionary between the unmangled and mangled names
diff -y --suppress-common-lines libfirst.unmangled libfirst.sym | awk '
{print $1, $3}
' > libfirst.dict
# remove temporary files
rm libfirst.unmangled libfirst.sym
