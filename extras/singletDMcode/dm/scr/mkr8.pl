#!/usr/bin/perl -w
#
# This script goes through the files given as arguments and changes explicit
# declarations of real to real*8 and complex to complex*16.

# First define list of names to change
# If the new name is '-', the new name is set to 'ds' plus the old name
%names=qw(real      real*8
          complex   complex*16
);  

while(defined($file=shift)) {
    print "Now fixing file $file\n";
    $outfile="tmp.f";

    open(FILE,$file) || die "Can't open $file\n";
    open(OUT,">$outfile") || die "Can't open $outfile for writing.\n";
    $ch=0;
    $i=0;
    while(defined($line=<FILE>)) {
        $i++;
        chomp($line);
        $ch += ($line =~ s/^      REAL /      REAL*8 /i);
        $ch += ($line =~ s/^      COMPLEX /      COMPLEX*16 /i);
        if (length($line)>72 and substr($line,0,1)=~ /\s/
           and not(substr($line,0,72) =~ /!/)) {
	    print "***Line $i is too long\n";
	}
   	print OUT "$line\n";
    }

    close(FILE);
    close(OUT);
    if ($ch>0) {
	print "  $file changed ($ch place", 
        ($ch==1) ? '' : "s",
        ") - replacing with new version\n";
        rename($outfile,$file);
        $ch=0;
    } else {
	unlink($outfile);
    }
}

