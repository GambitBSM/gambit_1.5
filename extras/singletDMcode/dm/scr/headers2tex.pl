#!/usr/bin/perl -w
#
# Script to go through the files in DarkSUSY and produce a tex file
# with the headers of all the routines
# Author: Joakim Edsjo, edsjo@physto.se
# Date: October 23, 2000.
#

@dirlist=qw(ac an an1l anstu dd ep ge ha hm hr ini mu nt pb rd 
   rn su xcern xcmlib);
$texfile="docs/headers.tex";

# Determine name of current DarkSUSY directory
$dsver=`pwd`;
chomp($dsver);
$dsver =~ s#^.*/##;
$dsver =~ s#^ds-##;
$date=localtime;

open(TEX,">$texfile") || die "Can't open $texfile for writing.\n";

print_texheader($dsver,$date);

# Start going through the directories and files

chdir("src") || die "Can't cd to src\n";
foreach $dir (@dirlist) {
    chdir($dir) || die "Can't cd to $dir\n";
    print "Directory: $dir\n";
    print_texdir($dir);
    @files=<*.f>;
    foreach $file (@files) {
        print "   File: $file\n";
	print_texfile($file);
    }
    chdir("..") || die "Can't cd to ..\n";
}

print_texend();
print "Done.\n";
exit;

##########
sub print_texheader {
    my $ver=$_[0];   # DarkSUSY version
    my $date=$_[1];

print TEX <<'END';
\documentclass{article}
% This LaTeX file is automatically created by the script
% headers2tex.pl that can be found in the scr directory in the
% DarkSUSY distribution.
% The script is written by Joakim Edsjo, edsjo@physto.se

\addtolength{\textwidth}{3cm}
\addtolength{\oddsidemargin}{-1.5cm}
\addtolength{\textheight}{1cm}
%\addtolength{\topmargin}{-0.5cm}

%\newcommand{\addtoindex}[1]{
%\addtocontents{toc}{#1\dotfill\arabic{page}\hskip2in\break}}

%\newenvironment{routine}[1]{\newpage\section*{#1}
%\hrule\vspace*{0.3in}\bgroup\addtoindex{#1}}
%{\egroup\vspace*{0.5in}}

\newenvironment{routine}[1]{\subsection*{#1}
\hrule\vspace*{1ex}\addcontentsline{toc}{subsection}{#1}}

\pagestyle{headings}

\begin{document}

\null
\bigskip
\bigskip

END

print TEX <<END;
\\centerline{\\LARGE \\bf DarkSUSY $ver}
END

print TEX <<'END';
\bigskip

\centerline{\LARGE \bf Routine headers}

\bigskip
\bigskip

\centerline{\Large Created automatically by headers2tex.pl}
\smallskip
END

print TEX <<END;
\\centerline{$date}
END

print TEX <<'END';
\bigskip

\newpage

\tableofcontents

END
}


##########
sub print_texdir{
    my $dir=$_[0];

print TEX <<END;

\\newpage
\\section{Directory src/$dir}

END
}


##########
sub print_texfile{
    my $file=$_[0];
    my $tfile=$file;
    my $line;
    my $tline;
    my $i;
    $tfile =~ s#\_#\\\_#g;

print TEX <<END;

%%%%% routine $file %%%%%
\\begin{routine}{$tfile}
\\begin{verbatim}
END

# Now go through the file

    open(IN,"$file") || die "Can't open $file for reading.\n";
    $i=0;
    while(defined($line=<IN>)) {
        if (substr($line,0,1) ne " ") {
	    print TEX $line;
            $i++;
            next;
	}
        next if $line =~ /subroutine/i;
        next if $line =~ /function/i;
        $tline=chomp($line);
	$tline =~ s/\s*//g;
        next if length($tline)==0;
        next if substr($line,5,1) ne " ";
        last;  # OK, now we have come to the code
    }
    close(IN);
    print TEX "No header found.\n" if $i==0;

print TEX<<'END';
\end{verbatim}
\end{routine}
END

}


##########
sub print_texend{

print TEX <<'END';
\end{document}
END

}


