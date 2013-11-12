#!/usr/bin/perl -w
# Script to install/convert yield tables

$file=shift;
$dir=shift;

if ($file =~ /.gz$/) {    # Convert ascii to binary
    $newfile = $file;
    $newfile =~ s/.dat.gz$/.bin/;
#    if ($newfile =~ /int/) {
#	$type="i";
#    } else {
#	$type="d";
#    }
    $tmpfile="tmp.dat";
    system("gunzip -v -c $file > $tmpfile");
# Now we are ready to run ascii2bin
    open(PROG,"|./ascii2bin") ||
	die "Can't start ascii2bin.\n";
#    print PROG "$type\n";
    print PROG "$tmpfile\n";
    print PROG "$newfile\n";
    close(PROG);
    unlink($tmpfile);
#    system("mv $newfile $dir/") &&
#	die "Can't move file $newfile to directory $dir\n";
    system("mv $newfile $dir/");
    system("chmod u+rw $dir/$newfile");
    system("chmod og+r $dir/$newfile");
} else {
    system("cp $file $dir") &&    
	die "Can't copy file $file to directory $dir\n";;
    system("chmod u+rw $dir/$file");
    system("chmod og+r $dir/$file");
}
