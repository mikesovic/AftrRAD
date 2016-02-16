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
	print "\nEnter the name of the Haplotypes file you want to use\n";

	$FileName = <STDIN>;
	chomp($FileName);
}



#Put pairs of replicate samples in hash. 

open REPLICATES, "Replicates.txt" or die$!;

my %ReplicateNames = ();

while (<REPLICATES>)  {
	
	if ($_ =~ /[a-zA-Z]/)  {
		
		chomp($_);
		my @TempArray = split(/\t/,$_);
		$ReplicateNames{$TempArray[0]} = $TempArray[1];
	}
}









################################################################################################################################
#Create file (named SlashFile) with each individual in one line with haplotypes separated by a /.

open HAPLOTYPES, "../Output/Genotypes/$FileName" or die$!;
system "mkdir TempFiles_Dup";

open INFILE, ">TempFiles_Dup/SlashFile.txt" or die$!;

my $LineCounter = 0;
my @Array1;
my @Array2;

while(<HAPLOTYPES>)  {
	
	chomp($_);
	
	if ($_ =~ /Locus/)  {
		print INFILE "$_\n";
		next;
	}
	
	if ($LineCounter == 0) {	#on the first line for a sample
		
		@Array1 = split (/\t/, $_);
		$LineCounter=1;
		next;
	}	
	
	if ($LineCounter == 1)  {	#on the second line for a sample
		
		
		@Array2 = split (/\t/, $_);
		$LineCounter=0;
		
		print INFILE "$Array1[0]\t";	#print the name of the current sample to the slash file
		
		my $ArrayLength = @Array1;
		
		
		foreach my $value (1..$ArrayLength-1)  {	#go through each locus
			
			
			if ($Array1[$value] =~ /NA/) {		#check to see if the current locus is missing
				print INFILE "NA\t";
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



close INFILE;
close HAPLOTYPES;



################################################################################################################################







system "cp TempFiles_Dup/SlashFile.txt TempFiles_Dup/SlashFile2.txt";


open HAPLOTYPES, "TempFiles_Dup/SlashFile.txt" or die$!;
open INFILE, ">Duplicate_Report.txt" or die$!;
print INFILE "Comparison\tSites_Compared\tNumber_Matches\tProportion_Matches\tMismatches_Homozygous_In_Sample1\tMismatches_Homozygous_In_Sample2\tLoci_With_Both_Reps_Homozygous\tLoci_With_At_Least_One_Het\n";
			

my $TotalNonNAComparisons = 0;
my $TotalMatches = 0;

while(<HAPLOTYPES>)  {
	
	chomp($_);
	
	if ($_ =~ /Locus/)  {
		#print INFILE "$_\n";
		next;
	}


	else {
		
		my @TempArray = split(/\t/,$_);
		my $CurrentName = $TempArray[0];
		$CurrentName =~ s/Individual//;
		
		if (exists($ReplicateNames{$CurrentName})) {		#the current sample in the slashfile is part of one of the replicate comparisons.
			
			print INFILE "$CurrentName/$ReplicateNames{$CurrentName}";
			my $CurrentIndividualNonNAComparisons = 0;
			my $CurrentIndividualMatches = 0;
			my $FirstSampleHomozygousLoci = 0;
			my $SecondSampleHomozygousLoci = 0;
			my $DoublyHomozygousComparisons = 0;
			
			my @FirstDup = split(/\t/,$_);		#store the genotypes from the current sample
			my @SecondDup = ();
			
			my $DupName = $ReplicateNames{$CurrentName};	#Get the name of the replicate sample.
			
			open HAPLOTYPESB, "TempFiles_Dup/SlashFile2.txt" or die$!;
			
			while (<HAPLOTYPESB>)  {
				
				if ($_ =~ /Locus/)  {
					#print INFILE "$_\n";
					next;
				}
				
				else {
					my @TempArray = split(/\t/,$_);
					my $TempName = $TempArray[0];
					$TempName =~ s/Individual//;
					
					if ($TempName eq $DupName)  {	#Find the replicate sample in the copied Slash file.  Assign this line to @SecondDup for direct comparison to @FirstDup.
						
						@SecondDup = @TempArray;
						last;
					}
				}
			}
			
			
			
			my $NumLoci = @FirstDup-1;
			#print "Number of loci is $NumLoci\n";
			
			
			foreach my $element (1..$NumLoci)  {
				
				
				my @CurrentHaplotype1 = split(/\//, $FirstDup[$element]);
				
				my $CurrentHaplotype1a = $CurrentHaplotype1[0];
				my $CurrentHaplotype1b = $CurrentHaplotype1[1];
				
			
				my @CurrentHaplotype2 = split(/\//, $SecondDup[$element]);
				
				my $CurrentHaplotype2a = $CurrentHaplotype2[0];
				my $CurrentHaplotype2b = $CurrentHaplotype2[1];
				
				if (($CurrentHaplotype1a =~ /NA/) || ($CurrentHaplotype2a =~ /NA/))  {	#missing data in at least one of the samples - don't compare
					next;
				}	
				
				$CurrentIndividualNonNAComparisons++;
				$TotalNonNAComparisons++;
				
				my $FirstHaplotype = $CurrentHaplotype1a.$CurrentHaplotype1b;  #use this as reference haplotype
				
				my $SecondHaplotypeA = $CurrentHaplotype2a.$CurrentHaplotype2b;		#Check the second haplotype against the first - have to check in both directions (i.e. A/G and G/A are the same genotype)
				my $SecondHaplotypeB = $CurrentHaplotype2b.$CurrentHaplotype2a;
				
				if (($FirstHaplotype eq $SecondHaplotypeA) || ($FirstHaplotype eq $SecondHaplotypeB))  {
					$CurrentIndividualMatches++;
					$TotalMatches++;
					
					if ($CurrentHaplotype1a eq $CurrentHaplotype1b)  {	#The two replicates match, and the locus is homozygous
						$DoublyHomozygousComparisons++;
					}	
				}
				
				elsif ($CurrentHaplotype1a eq $CurrentHaplotype1b)  {		#There is a mismatch, and the first replicate is homozygous
					$FirstSampleHomozygousLoci++;
				}
				
				elsif ($CurrentHaplotype2a eq $CurrentHaplotype2b)  {		#There is a mismatch, and the second replicate is homozygous
					$SecondSampleHomozygousLoci++;
				}	
					
					
				
			}
			
			
			#print INFILE "Comparison\tSites_Compared\tNumber_Matches\tProportion_Matches\tMismatches_Homozygous_In_Sample1\tMismatches_Homozygous_In_Sample2\tLoci_With_Both_Reps_Homozygous\tLoci_With_At_Least_One_Het\n";
			#Sites compared ($CurrentIndividualNonNAComparisons) is the number of sites that were genotyped in both samples (neither had missing data).
			#Number_Matches ($CurrentIndividualMatches) is the total number of sites (of those compared - no missing data) that were scored the same.
			#Proportion_Matches ($CurrentIndividualProportionCorrect)  Number matches divided by sites compared
			#Mismatches_Homozygous_In_Sample1 ($FirstSampleHomozygousLoci)  Of the mismatched loci, the number homozygous in sample 1
			#Mismatches_Homozygous_In_Sample2 ($SecondSampleHomozygousLoci)  Of the mismatched loci, the number homozygous in sample 2
			#Loci_With_Both_Reps_Homozygous ($DoublyHomozygousComparisons)  Of the loci compared, the number that were homozygous in both samples.
			#Loci_With_At_Least_One_Het ($HeterozygousLociForDropout)  Of the loci compared, the number that were heterozygous in at least one of the two samples.
			
			my $CurrentIndividualProportionCorrect = $CurrentIndividualMatches/$CurrentIndividualNonNAComparisons;
			
			
			print INFILE "$CurrentIndividualNonNAComparisons\t$CurrentIndividualMatches\t$CurrentIndividualProportionCorrect\t$FirstSampleHomozygousLoci\t$SecondSampleHomozygousLoci\t";
				
			my $HeterozygousLociForDropout = $CurrentIndividualNonNAComparisons-$DoublyHomozygousComparisons;
			
			print INFILE "$DoublyHomozygousComparisons\t$HeterozygousLociForDropout\n";	
				
				
		}
	}

}	


my $TotalProportionCorrect = $TotalMatches/$TotalNonNAComparisons;


print INFILE "\n\nAcross all samples, there were $TotalNonNAComparisons loci compared.  Of these, $TotalMatches ($TotalProportionCorrect) were scored the same in the replicate samples.\n";

#print INFILE "$TotalNonNAComparisons\t$TotalMatches\t$TotalProportionCorrect\n";




















	
			
		
		
	
	

	
