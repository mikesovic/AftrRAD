#! /usr/bin/perl
use warnings;
use strict;


#This script creates a Genepop input file from a Haplotypes file.  Note that the 'POP' flags required in a Genepop input file must be added by hand after running this script.

#First, check to see how many Haplotypes files exist.  If only one, use this.  If more than one, prompt the user to choose one.

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
	print "\nEnter the name of the Haplotypes file you want to use to create the Genepop infile\n";

	$FileName = <STDIN>;
	chomp($FileName);
}




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
open OUTFILE, ">GenepopInputTemp.txt" or die$!;

print OUTFILE "Genepop Input for file $FileName\n";

my $LineCounter = 0;
my $NumLoci;
my @TempArray1;
my @TempArray2;
my $FirstLine = 1;

while(<INFILE>) {
	
	chomp($_);
	
	if ($FirstLine == 1)  {		#On the line with locus names.
		$_ =~ s/"//g;
		$_ =~ s/^\s//;
		my @TempArrayNames = split(/\t/,$_);
		
		foreach my $locusname (@TempArrayNames)  {
			print OUTFILE "$locusname, ";
		}
		
		print OUTFILE "\n";
		
		$FirstLine = 0;
		
		next;
	}	
		
	
	
	if ($LineCounter == 0)  {	#Each (diploid) sample is represented with two lines.  Currently on the first of the two.
		$_ =~ s/"//g;
		$_ =~ s/^\s//;
		@TempArray1 = split(/\t/,$_);
		$NumLoci = @TempArray1;
		$LineCounter++;
		next;
	}

	
	if ($LineCounter == 1)  {	#Now on the second line/allele of the current sample.
		$_ =~ s/"//g;
		$_ =~ s/^\s//;
		@TempArray2 = split(/\t/,$_);
		
		print OUTFILE "$TempArray1[0], ";
		
		foreach my $locus (1..$NumLoci-1)  {
			
			my $TempValue1 = $TempArray1[$locus];
			my $TempValue1ToPrint = $TempValue1+1;
			
			my $TempValue2 = $TempArray2[$locus];
			my $TempValue2ToPrint = $TempValue2+1;
		
		
			if ($TempValue1ToPrint == -8)  {	#adjusted this from -9 because 1 was added above.
				print OUTFILE "00";
			}	
				
			elsif ($TempValue1ToPrint < 10)  {	#every allele in the genepop file is represented by two digits.
				print OUTFILE "0$TempValue1ToPrint";
			}
		
			else {
				print OUTFILE "$TempValue1ToPrint";
			}
		
		
			if ($TempValue2ToPrint == -8)  {
				print OUTFILE "00\t";
			}
			
			elsif ($TempValue2ToPrint < 10)  {
				print OUTFILE "0$TempValue2ToPrint\t";
			}
		
			else {
				print OUTFILE "$TempValue2ToPrint\t";
			}
		
		
		
		}

		$LineCounter = 0;
		print OUTFILE "\n";
		
	}
}

close INFILE;
close OUTFILE;


open INFILE, "GenepopInputTemp.txt" or die$!;
open OUTFILE, ">GenepopInput_Temp2.txt" or die$!;

my $LineNumber = 0;

while(<INFILE>) {
	
	if ($LineNumber == 1) {
		chomp($_);
		$_ =~ s/\s$//;
		$_ =~ s/,$//;
		print OUTFILE "$_\n"
	}

	else {
		print OUTFILE "$_";
	}
	
	$LineNumber++;
}

close INFILE;
close OUTFILE;

system "rm GenepopInputTemp.txt";
system "rm StructureInput_Temp.txt";
system "rm OutputStructure_Infile_Temp.txt";
system "rm OutputStructure_Infile.txt";


my $SortCounter = 0;


	
	open INFILE, "GenepopInput_Temp2.txt" or die$!;
	open OUTFILE, ">GenepopInput.txt" or die$!;
	
	my @UnsortedArray = ();
	
	while (<INFILE>) {
		if ($SortCounter <2) {
			print OUTFILE "$_";
			$SortCounter++;
		}
		
		else {
			chomp($_);
			my($first, $rest) = split(/\t/, $_, 2);
			push (@UnsortedArray, $first);
			push (@UnsortedArray, $rest);
		}
		
	}
	
	

	my %HashToSort = @UnsortedArray;

	my @SortedNames = sort { lc($a) cmp lc($b) } keys %HashToSort;

	close INFILE;
	
	
	foreach my $name (@SortedNames)  {
		
		open INFILE, "GenepopInput_Temp2.txt" or die$!;
		
		while (<INFILE>)  {
			
			chomp($_);
			my @TempArray = split(/\t/, $_);
			
			if ($TempArray[0] eq $name)  {
				$_ =~ s/Individual//;
				print OUTFILE "$_\n";
			}
		}
		
		close INFILE;
	}	
	
	close OUTFILE;

system "rm GenepopInput_Temp2.txt";	
	



