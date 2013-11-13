#!/usr/bin/perl -w 

# Script by Joakim Edsjo to tar the DarkSUSY distribution The output is
# written to ~/ut/ Call this script from the root of the DarkSUSY
# distribution as as scr/tarit.pl. The script will determine the
# DarkSUSY version from the name of the DarkSUSY root directory, which
# will give the name to the tar file.

# Determine name of current DarkSUSY directory

$dir=`pwd`;
chomp($dir);
$dir =~ s#^.*/##;
$tarfile="$dir.tar.gz";

# Check if ~/ut exists
mkdir("~/ut",0777) unless -d "~/ut";

# Tar it up
chdir "..";
print "Tarring it up into ~/ut/$tarfile\n";
system("tar zcvf ~/ut/$tarfile $dir/README* $dir/HISTORY"
       ." $dir/inc $dir/scr $dir/lib $dir/dat"
       ." $dir/make* $dir/src/make* $dir/src/*/make* "
       ." $dir/src/*/*.f $dir/src/*/*.h $dir/test");

exit


