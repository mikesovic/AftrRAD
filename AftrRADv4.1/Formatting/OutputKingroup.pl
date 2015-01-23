#! /usr/bin/perl
use warnings;
use strict;


my $FileName;

opendir GENOS, "../Output/Genotypes/";
my @AllFileNames = grep { $_ ne '.' && $_ ne '..' && $_ ne '.DS_Store' } readdir(GENOS);
close GENOS;

my @HaplotypeFileNames = ();

for my $name (@AllFileNames)  {
	if ($name =~ /Haplo/) {
		push (@HaplotypeFileNames, $name);
	}
}

my $NumHaplotypeFiles = @HaplotypeFileNames;

if ($NumHaplotypeFiles == 1)  {
	$FileName = $HaplotypeFileNames[0];
}

else {
	print "The following Haplotypes files are available in the Output/Genotypes directory...\n";
	for my $file (@HaplotypeFileNames) {
		print "$file\n";
	}
	print "\nEnter the name of the Haplotypes file you want to use to create the Kingroup infile\n";

	$FileName = <STDIN>;
	chomp($FileName);
}	



# print "Enter the name of the haplotypes file to use\n";
# my $FileName = <STDIN>;
# chomp($FileName);


open HAPLOTYPES, "../Output/Genotypes/$FileName" or die$!;
open INFILE, ">Kingroup_Infile.txt" or die$!;

my $LineCounter = 0;
my @Array1;
my @Array2;

while(<HAPLOTYPES>)  {
	
	chomp($_);
	
	if ($_ =~ /Locus/)  {
		print INFILE "$_\n";
		next;
	}
	
	if ($LineCounter == 0) {
		
		@Array1 = split (/\t/, $_);
		$LineCounter=1;
		next;
	}	
	
	if ($LineCounter == 1)  {
		
		
		@Array2 = split (/\t/, $_);
		$LineCounter=0;
		
		print INFILE "$Array1[0]\t";
		
		my $ArrayLength = @Array1;
		
		
		foreach my $value (1..$ArrayLength-1)  {
			
			
			if ($Array1[$value] =~ /NA/) {
				print INFILE "0\t";
			}	
			
			else {
				print INFILE "$Array1[$value]";
				print INFILE "/";
				print INFILE "$Array2[$value]\t";
			
			}
		}
		
		print INFILE "\n";
	}
}	
			
		
		
	
	

	
