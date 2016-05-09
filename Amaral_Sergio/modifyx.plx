#!/usr/bin/perl
# modifyx.plx
use warnings;
use strict;

#Written by Bogdan Constantin
#This code will "xplexify" existing fortran files by:
#	1: adding "USE XPLEXIFY" in the USE statement block
#       2: adding "USE MYTYPE" in the USE statement block
#NOTE!!! This might 99% cause errors in compilation! ( don't worry! very fixable very fast!) 


foreach (@ARGV){
# Get the name of the file currently read from the command prompt
my $file = $_;
print "$file\n";

# Create a temporary file in which the modified code will be written
my $temp = "$file.tmp";
print "$temp\n";

# Open the file to read 
open FILE, "<$file" or die $!;
# Open the file to write 
open TEMP, ">$temp" or die $!;

# initialize counter
my $lineno=1;

# loop over the read file lines and write it to the temporary file 
# with the modified lines :)
while (<FILE>) {
	my $string = $_;
		if ($string =~ m/IMPLICIT NONE/i ) {
			$string =~ s/IMPLICIT NONE/USE XPLEXIFY\n       USE MYTYPE/i;
			print TEMP "$string";
		}
		if ($string =~ m/REAL\*8/i){
			$string =~ s/REAL\*8/TYPE (XPLEX)/i;
			$_ = $string;
		}
		if ($string =~ m/REAL\*4/i){
                        $string =~ s/REAL\*4/TYPE (XPLEX)/i;
                        $_ = $string;
                }
                if ($string =~ m/double precision/i){
                        $string =~ s/double precision/TYPE (XPLEX)/i;
                        $_ = $string;
                }
                if ($string =~ m/real\*8/i){
                        $string =~ s/real\*8/TYPE (XPLEX)/i;
                        $_ = $string;
                }
                if ($string =~ m/real\*4/i){
                        $string =~ s/real\*4/TYPE (XPLEX)/i;
                        $_ = $string;
                }
                if ($string =~ m/\Qreal(kind=dp)\E/i){
                        $string =~ s/\Qreal(kind=dp)\E/\QTYPE (XPLEX)\E/i;
                        $_ = $string;
                }
                if ($string =~ m/\QREAL(kind=dp)\E/i){
                        $string =~ s/\QREAL(kind=dp)\E/\QTYPE (XPLEX)\E/i;
                        $_ = $string;
                }
                if ($string =~ m/\QREAL(KIND=dp)\E/i){
                        $string =~ s/\QREAL(KIND=dp)\E/\QTYPE (XPLEX)\E/i;
                        $_ = $string;
                }
	       	if ($string =~ m/\Qreal(kind=8)\E/i){
                        $string =~ s/\Qreal(kind=8)\E/\QTYPE (XPLEX)\E/i;
                        $_ = $string;
                }
                if ($string =~ m/\QREAL(kind=8)\E/i){
                        $string =~ s/\QREAL(kind=8)\E/\QTYPE (XPLEX)\E/i;
                        $_ = $string;
                }
                if ($string =~ m/\QREAL(KIND=8)\E/i){
                        $string =~ s/\QREAL(KIND=8)\E/\QTYPE (XPLEX)\E/i;
                        $_ = $string;
                }
	print TEMP "$_";
}
# rename the temporary file as the original one! 
rename ("$temp","$file");

#close the files
close FILE;
close TEMP;

}
