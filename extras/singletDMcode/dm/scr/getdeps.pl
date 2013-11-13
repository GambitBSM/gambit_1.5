#!/usr/bin/perl -w
#
# Script to go through and find all needed dependencies in
# the given files.
#
# Author: Joakim Edsjo, edsjo@physto.se
# Date: August 29, 2000.

while(defined($file=shift)) {
    open(FILE,"$file");
    while(defined($line=<FILE>)) {
	if ($line =~ /include\s+'(.+)'/i) {
	    $deps{$1}++;
	}
    }
    close(FILE);
}

foreach $dep (sort keys %deps ) {
    print "$dep\n";
}
