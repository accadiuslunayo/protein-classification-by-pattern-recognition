#!/usr/bin/perl-w
use strict;
use warnings;
use Cwd;
use File::Path;

my $cwd = cwd();
print "The current working directory is:\n $cwd\n\n\n";

my $output_dir = "output";
unless (-e $output_dir and -d $output_dir) {
    mkdir $output_dir or die "Unable to create output directory: $!";
    print "Output directory created in the current working directory\n";
}
chdir $output_dir;

#defining the format of the files needed
print "\n====**File type =.txt** ========= sequence format = .fasta \n\n";
#Get the name of the file with the protein sequence data.
print "Please enter filename of the protein sequence data:\n\n";
$filename = <>;
print "Filetype == protein from the $filename\n\n";
my $proteinfile = <>;
#Remove the newline from the protein filename 
chomp $proteinfile;
#Open the file, or exit. Filehandle PROTEINFILE chosen for readability.
open (my $PROTEINFILE,"$proteinfile") or die ("Error in opening results file");
#Initialize the amino_acids using the amino acid symbols.
    my $Ala = 0;my $Arg = 0;my $Asn = 0;my $Asp = 0;my $Cys = 0;
    my $Gln = 0;my $Glu = 0;my $Gly = 0;my $His = 0;my $Ile = 0;
    my $Leu = 0;my $Lys = 0;my $Met = 0;my $Phe = 0;my $Pro = 0;
    my $Ser = 0;my $Thr = 0;my $Trp = 0;my $Tyr = 0;my $Val = 0;
#skipping the first sline of each sequence beginning with > sign
$/ = ">";
#store the information 
my $junk = (<$PROTEINFILE>);
#store the location of each file
#initialize the location === a change is required every time you enter data
#proteinfile source is  data
print "Please enter the location from where the data is obtained"
my $name = '<>';
#create the outfile name and format
open(POS,">Result_.csv") or die "If cant open file$!";
while ( my $record = (<$PROTEINFILE>) ) {chomp $record; 
	my ($defLine, @seqLines) = split / /, $record;
	my $sequence = join('',@seqLines);            
	if ($Ala=$sequence =~ s/A//g){}if ($Arg=$sequence =~ s/R//g){}if ($Asn=$sequence=~ s/N//g) {}
    if ($Asp=$sequence =~ s/D//g){}if ($Cys=$sequence =~ s/C//g){}if ($Gln=$sequence =~ s/Q//g){}
	if ($Glu=$sequence =~ s/E//g){}if ($Gly=$sequence =~ s/G//g){}if ($His=$sequence =~ s/H//g){}
    if ($Ile=$sequence=~ s/I//g) {}if ($Leu=$sequence =~ s/L//g){}if ($Lys=$sequence =~ s/K//g){}
	if ($Met=$sequence=~ s/M//g) {}if ($Phe=$sequence=~ s/F//g) {}if ($Pro=$sequence =~ s/P//g){}
    if ($Ser=$sequence=~ s/S//g) {}if ($Thr=$sequence =~ s/T//g){}if ($Trp=$sequence =~ s/W//g){}
    if ($Tyr=$sequence =~ s/Y//g){}if ($Val=$sequence=~ s/V//g)  {}
	                       
#calculate the sequence length
my $Seq_Length =$Ala + $Arg + $Asn + $Asp + $Cys + $Gln +
                $Glu + $Gly + $His + $Ile + $Leu + $Lys + 
                $Met + $Phe + $Pro + $Ser + $Thr + $Trp + 
                $Tyr + $Val;
# find amino acid percentage frequency
my $AlaP =($Ala/$Seq_Length)*100;my $ArgP =($Arg/$Seq_Length)*100;my $AsnP =($Asn/$Seq_Length)*100;
my $AspP =($Asp/$Seq_Length)*100;my $CysP =($Cys/$Seq_Length)*100;my $GlnP =($Gln/$Seq_Length)*100;
my $GluP =($Glu/$Seq_Length)*100;my $GlyP =($Gly/$Seq_Length)*100;my $HisP =($His/$Seq_Length)*100;
my $IleP =($Ile/$Seq_Length)*100;my $LeuP =($Leu/$Seq_Length)*100;my $LysP =($Lys/$Seq_Length)*100;
my $MetP =($Met/$Seq_Length)*100;my $PheP =($Phe/$Seq_Length)*100;my $ProP =($Pro/$Seq_Length)*100;
my $SerP =($Ser/$Seq_Length)*100;my $ThrP =($Thr/$Seq_Length)*100;my $TrpP =($Trp/$Seq_Length)*100;
my $TyrP =($Tyr/$Seq_Length)*100;my $ValP =($Val/$Seq_Length)*100;
#catecorize the percentage amino acids based on polarity
my $group1_AlNP = 0;$group1_AlNP = $AlaP+$IleP+$LeuP+$ValP;
my $group2_ArNP = 0;$group2_ArNP = $PheP+$TrpP+$TyrP;
my $group3_NeP  = 0;$group3_NeP  = $AsnP+$CysP+$GlnP+$MetP+$SerP+$ThrP;
my $group4_AP   = 0;$group4_AP   = $GluP+$AsnP;
my $group5_BP   = 0;$group5_BP   = $ArgP+$HisP+$LysP;
my $group6_Uniq = 0;$group6_Uniq = $GlyP+$ProP;
#To find out amino acid molecular weight
my $mole_weight= ($Ala*89)+($Arg*174)+($Asn*132)+($Asp*133)+($Cys*121)+($Gln*146)
				+($Glu*147)+($Gly*75)+($His*155)+($Ile*131)+($Leu*131)+($Lys*146)
				+($Met*149)+($Phe*165)+($Pro*115)+($Ser*105)+($Thr*119)+($Trp*204)
				+($Tyr*181)+($Val*117);
#find unique protein compositon
my $unique = 0; $unique = ($Gly * -0.4)+($Pro * -1.6);
#To find the hydrophobicity index for the sequences
my $hydrophilic=0;
$hydrophilic = ($Met * 1.9)+($Ser * -0.8)+($Thr * -1.3)+($Cys * 2.5)+($Asn * -3.5)	
			   +($Gln * -3.5)+($Lys * -3.9)+($His * -3.2)+($Arg * -4.5)+($Asp * -3.5)+($Glu * -3.5);
#to find the hydrophicility index for the sequences
my $hydrophob = 0;
$hydrophob =($Ala * 1.8)+($Leu * 3.8)+($Ile * 4.5)+($Val * 4.2)+($Phe * 2.8)+($Tyr * -1.3)+($Trp * -0.9);
#To find the isoelectricity index for the sequences
 my $isoel=0;
$isoel = ($Gly * (5.97))+($Ala* 6.01)+($Pro * 6.48)+($Val * 5.97)+($Leu * 5.98)+($Ile * (6.02)) +($Met * 5.74)		
		+($Phe * 5.48)+($Tyr * 5.66)+($Trp * 5.89)+($Ser * 5.68)+($Thr * 5.87)+($Cys * 5.07)+($Asn * 5.41)			
		+($Gln * 5.65)+($Lys * 9.74)+($His * 7.59)+($Arg * 10.76)+($Asp * 2.77)+($Glu * 3.22);
print  POS $defLine.",".","."$Seq_Length".","."$mole_weight".","."$hydrophob".","."$hydrophilic".","."$unique".","."$isoel".",";
print  POS $group1_AlNP.",".$group2_ArNP.",".$group3_NeP.",".$group4_AP.","$group5_BP.","$group6_Uniq.",".$name;
							print POS "\n";
}
print "=====!!The file is printed and located in the current working directory:====\n\n"; 
print "\n $cwd\n\n";
print "======== Please press Enter======= \n\n\n\n";
