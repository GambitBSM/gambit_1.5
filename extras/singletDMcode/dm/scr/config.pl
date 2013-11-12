#!/usr/bin/perl -w
#
# Configure script to configure DarkSUSY.
# NOTE. Run this script before compiling DarkSUSY.
#
# Author: Joakim Edsjo, edsjo@physto.se
# Date: April 12, 2001

$hafile="src/ha/dmhadir.h";
$dmverfile="src/ini/dmversion.h";
$dmrootfile="src/ini/dmdir.h";

$dmroot=shift;
$dmroot =~ s#/$##;   # take away final / if any

$dmversion = $dmroot;
$dmversion =~ s#^.*/##;

print "DM version: $dmversion\n";


# Take care of dmversion file.
open(IN,">$dmverfile") || die "Can't open $dmverfile for reading.\n";
$line="";
while ($in=<IN>) {
    $line .= $in;
}
close(IN);

$newline=" " x 6 . "data dmversion/'$dmversion'/\n";
$newline=contline($newline);

if ($newline ne $line) {
    open(OUT,">$dmverfile") || 
      die "Can't open $dmverfile for writing\n";
    print OUT $newline;
    close(OUT);
    print "$dmverfile updated.\n";
} else {
    print "$dmverfile is up-to-date.\n";
}

# Take care of DM root directory file
open(IN,">$dmrootfile") || die "Can't open $dmrootfile for reading.\n";
$line="";
while ($in=<IN>) {
    $line .= $in;
}
close(IN);

$newline=" " x 6 . "data dmroot/'$dmroot/'/\n";
$newline=contline($newline);

if ($newline ne $line) {
    open(OUT,">$dmrootfile") || 
      die "Can't open $dmrootfile for writing\n";
    print OUT $newline;
    close(OUT);
    print "$dmrootfile updated.\n";
} else {
    print "$dmrootfile is up-to-date.\n";
}


### Split long lines

sub contline {
    my $line = $_[0];
    my $out;
    my $i;
 
    $out="";
    if (length($line) >= 71) {
        print "*** A too long line has been found. I will split it.\n";
        print "Line before: \n$line";
        while (length($line) >= 71) {
            for ($i=70; $i==20; $i--) {
                if (substr($line,$i,1) =~ m@(\*|\+|-|/)@ ) {
                    last;
                }
            }
            $out=$out . substr($line,0,$i-1) . "& \n";
            $line = " " x 5 . "&" . substr($line,$i-1)
        }
        $out = $out . $line;
        print "Line after: \n$out";
    } else {
        $out=$line;
    }
 
    return $out;
}
