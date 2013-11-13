#!/usr/bin/perl -w
#
# This script checks if the first argument is present in the files given
# as the remaining arguments. If it is not present, the filenames are
# printed.

$text=shift;
print "Now going through the files searching for '$text'\n";

while(defined($file=shift)) {
    open(FILE,$file) || die "Can't open $file\n";
    $found=0;
    while(defined($line=<FILE>)) {
        $found++ if ($line =~ /$text/i);
    }
    close(FILE);
    if ($found==0) {
        print "Text NOT found in $file\n";
    }
}
