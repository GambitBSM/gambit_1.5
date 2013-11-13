#!/usr/bin/perl -w
#
# This script goes through the files given as arguments and changes names
# according to the table in the hash %names below. The names on the left
# are changed to the names on the right. If the newname is set to - this
# has a special meaning; the new name is 'ds' plus the old name appended to
# it.

# First define list of names to change
# If the new name is '-', the new name is set to 'ds' plus the old name
%names=qw(edsjo@cpfa.berkeley.edu    edsjo@physto.se
);  
# Now go through name changes and add ds where appropriate
foreach $name (keys %names) {
    $newname=$names{$name};
    if ($newname eq "-") {
        $newname = "ds$name";
        $names{$name}=$newname;
    }
}


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
        $ch += ($line =~ s/edsjo\@cfpa.berkeley.edu/edsjo\@physto.se/);
        $ch += ($line =~ s/edsjo\@teorfys.uu.se/edsjo\@physto.se/);
#	 $ch += ($line =~ s/'halo.h'/'hmcom.h'/);
#	 $ch += ($line =~ s/'ddcom.h'/'ddset.h'/);
#        $ch += ($line =~ s/'susy.h'/'dssusy.h'/);
#        $ch += ($line =~ s/double precision/real*8/i);
#        $ch += ($line =~ s/(\W)gef_int(\W|$)/$1dsf_int$2/ig);
#        $ch += ($line =~ s/(\W)gef_int2(\W|$)/$1dsf_int2$2/ig);
#        $ch += ($line =~ s/(\W)absq/$1dsabsq/gi);
#        $ch += ($line =~ s/\s+$//);
#Now change names to ds...
        foreach $name (keys %names) {
	    $newname=$names{$name};
            $ii=1;
            $jj=0;
            while ($ii > 0) {
		$ii = ($line =~ s/(\W|^)$name(\W|$)/$1$newname$2/gi);
		$jj += $ii;
	    }
	    $ch += $jj;
        }
#        $line=lc($line);   # uncomment if you want lowercase
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
        $newfile=$file;
        $ch=0;
        foreach $name (keys %names) {
	    $newname=$names{$name};
            $ch += ($newfile =~ s/(\W|^)$name(\W|$)/$1$newname$2/gi);
        }
        if ($ch>0) {
            print "  filename changed to $newfile\n";
            rename($file,$newfile);
	}

    } else {
	unlink($outfile);
    }
}

