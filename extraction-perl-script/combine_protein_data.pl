#!/usr/bin/perl-w
use strict;
use warnings;
use Cwd;
my $cwd = cwd();
print "The current working directory is:\n $cwd\n\n\n";
chdir "C:/../../../output";

{
	my @protein_files = glob('*.csv*');

print "\n\n=======Final output file format .csv======\n\n\n";
my $finalproteinfile = "Protein_file.csv" ;
open(OUTPUT,">".$finalproteinfile ) or die "If cant open file$!";

foreach my $csvfile(@protein_files)

{open (INPUT, $csvfile ) or die "If cant open file$!";
print OUTPUT <INPUT>;
close (INPUT);
}
close (OUTPUT);

print "=====!!The file is printed and located in the current working directory:====\n\n"; 


print "\n $cwd\n\n";


print "======== Please press Enter======= \n\n\n\n";