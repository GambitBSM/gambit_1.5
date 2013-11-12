#!/usr/bin/perl -w
#
# Simple script to fix lines like x**-y to x**(-y)
# Joakim Edsjo

$file=shift;
$outfile="$file-new";
open(FILE,$file) || die "Can't open $file\n";
open(OUT,">$outfile") || die "Can't open $outfile for writing.\n";
while(defined($line=<FILE>)) {
    $line =~ s/(\*\*)(-\d)(\D)/$1\($2\)$3/g;
    print OUT $line;
}

