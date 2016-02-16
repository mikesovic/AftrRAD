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
	print "\nEnter the name of the Haplotypes file you want to use to create the Structure infile\n";

	$FileName = <STDIN>;
	chomp($FileName);
}	


# print "Enter the name of the haplotype matrix file (in Results directory) you want to use. This file will be used to determine what loci are included in the Structure input. (include '.txt' extension)\n";
# my $FileName = <STDIN>;
# chomp($FileName);

#Replace any NA's in Haplotypes file with "N"

my $Counter1 = 0;

open HAPLOTYPEFILEWITHNA, "../Output/Genotypes/$FileName" or die$!;
open STRUCTUREFILENONA, ">OutputStructure_Infile_Temp.txt" or die$!;

	while (<HAPLOTYPEFILEWITHNA>)  {
		chomp($_);	

		if ($Counter1 == 0)  {
			$Counter1++;
			print STRUCTUREFILENONA "$_";
		}
		
		else {
			my @TempArray = split(/\t/, $_);
			my $Length = @TempArray;
			print STRUCTUREFILENONA "$TempArray[0]\t";
			
			foreach my $allele (@TempArray[1..$Length-1])  {
				$allele =~ s/NA/N/g;
				print STRUCTUREFILENONA "$allele\t";
			}
		}
		
		print STRUCTUREFILENONA "\n";
	}

close HAPLOTYPEFILEWITHNA;
close STRUCTUREFILENONA;


#Remove extra tabs at the end of each line. Final file is named "OutputStructure_Infile.txt".

open STRUCTURETEMP, "OutputStructure_Infile_Temp.txt" or die$!;
open STRUCTUREOUT, ">OutputStructure_Infile.txt" or die$!;

while (<STRUCTURETEMP>)  {
	$_ =~ s/\t$//;
	print STRUCTUREOUT "$_";
}

close STRUCTURETEMP;
close STRUCTUREOUT;




system "R --vanilla --slave < ../RScripts/OutputStructure.R";


open INFILE, "StructureInput_Temp.txt" or die$!;
open OUTFILE, ">StructureInput_Temp2.txt" or die$!;

my $LineCounter = 0;

while(<INFILE>) {
	if ($LineCounter == 0)  {
		$_ =~ s/"//g;
		$_ =~ s/^\s//;
		print OUTFILE "$_";
		$LineCounter++;
	}
	
	else {
		$_ =~ s/"//g;
		print OUTFILE "$_";
	}
}

close INFILE;
close OUTFILE;
	

system "rm OutputStructure_Infile_Temp.txt";
system "rm OutputStructure_Infile.txt";
system "rm StructureInput_Temp.txt";



	
	
	
open INFILE, "StructureInput_Temp2.txt" or die$!;
open OUTFILE, ">StructureInput.txt" or die$!;
while (<INFILE>)  {
	chomp($_);
	$_ =~ s/Individual//;
	print OUTFILE "$_\n";
}	
close INFILE;
close OUTFILE;

system "rm StructureInput_Temp2.txt";	
	
