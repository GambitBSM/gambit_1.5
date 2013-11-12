#!/usr/bin/perl
#
# Script to search for a variable in Fortran files

$var=shift;

while (defined($file=shift)) {
    $lno=0;
    $printed=0;
    open(IN,"$file") || die "Can't open $file for reading.\n";
    while (defined($line=<IN>)) {
	$lno++;
        if ($line =~ /\W${var}\W/ ) {
	    print "File: $file\n" if $printed==0;
	    print "Line $lno: $line";
            $printed++;
	}
    }  
}
