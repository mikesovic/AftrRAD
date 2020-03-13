#! /usr/bin/perl
use warnings;
use strict;

#Before running this script, edit the file "Genotypes/Output/SNPMatrix_X.Y.txt" by reordering individuals into populations and keep track of the number of individuals per population.

#Built this script from OutputFastSimCoal.pl.

print "\nInformation for OutputTreeMix.pl\n";
print "\nBefore running this script, make sure...\n";
print "1) samples in the file \"Genotypes/Output/SNPMatrix_X.Y.txt\" are ordered by population\n";
print "2) you know the number of individuals in each population, in the order they occur in the SNPMatrix file\n";


#Set parameters for the run

# -linked Flag indicating whether to include all SNPs, or only unlinked SNPs in the TreeMix file.  Default is 0 (includes only unlinked SNPs).

my %RunArguments = ();

#Defaults
$RunArguments{linked} = 0;

#read in command line arguments

for my $argument (@ARGV)  {
	my @TempArgs = split(/-/,$argument);
	
	#print help information
	if (($TempArgs[1] =~ /^h$/)|| ($TempArgs[1] =~ /^help$/)) {
		
		print "\nInformation for OutputTreeMix.pl...\n\n";
		print "This script creates a file that can be used as input for analyses in TreeMix\n";  
		print "It requires a SNPMatrix file in the Output/Genotypes folder. The samples in this file must be grouped by population.\n\n";
		
		print "\n\nCommand line arguments available for this script...\n\n";
		print "linked\nFlag indicating whether to include all SNPs, or only unlinked SNPs.\n";
		print "Default is 0 (only unlinked SNPs are included in the infile)\n\n";
		
		exit;
		
	}
	
	#update default values with entered parameters as appropriate and continue with run
	elsif (exists($RunArguments{$TempArgs[0]}))  {
		$RunArguments{$TempArgs[0]} = $TempArgs[1];
	}
}


#Set command line arguments for the run and ask which to use if more than 1.

my $Linked = $RunArguments{linked};

my $FileName;

#read the SNPMatrix files available for the run

opendir GENOS, "../Output/Genotypes/";
my @AllFileNames = grep { $_ ne '.' && $_ ne '..' && $_ ne '.DS_Store' } readdir(GENOS);
close GENOS;

my @SNPFileNames = ();

for my $name (@AllFileNames)  {
	if ($name =~ /SNPMatrix/) {
		push (@SNPFileNames, $name);
	}
}

my $NumSNPFiles = @SNPFileNames;

if ($NumSNPFiles == 1)  {
	$FileName = $SNPFileNames[0];
}

else {
	
		print "The following SNPMatrix files are available in the Output/Genotypes directory...\n";
		for my $file (@SNPFileNames) {
			print "$file\n";
		}
		print "\nEnter the name of the SNPMatrix file you want to use to create the TreeMix infile.\n";
	
		$FileName = <STDIN>;
		chomp($FileName);
	
}




print "\nEnter the number of individuals in each population, in the order the populations occur in the SNPMatrix_X.Y.txt file.  Separate each with a tab.\n";
my $PopSizes = <STDIN>;
chomp($PopSizes);



mkdir "TempFiles" unless (-d "TempFiles");

#Create array that contains the number of individuals (assumes diploid) in each population.
my @HapSizes = split(/\t/,$PopSizes);

#Create array that contains the number of chromosomes sampled at each locus in each population.
my @DiploidSizes = ();

foreach my $size (@HapSizes)  {
	my $TempSize = $size*2;
	push (@DiploidSizes,$TempSize);
}

my $NumPops = @DiploidSizes;

if ($NumPops == 1)  {
	
	print "Warning: Only one population recognized\n";
	print "Need multiple populations for TreeMix file.";
	exit;
}

else {
	print "Recognized $NumPops populations.\n";  
}	
	

print "\n\n";


			
			
print "\n\nWorking on creating TreeMix file\n\n";
				
##############################################################################################################################################
				
				
#Clean up the names in SNPMatrix file.
				
				
open FILE, "../Output/Genotypes/$FileName" or die$!;
open OUTFILE, ">TempFiles/SNPMatrix_Edit.txt" or die$!;
				
my @LocusNames = ();
my $LineNumber = 0;
				
while(<FILE>)  {
					
	if ($LineNumber == 0)  {
		@LocusNames = split(/\t/,$_);
		$LineNumber++;
	}
					
	$_ =~ s/Individual//;
	print OUTFILE "$_";
}
				
shift(@LocusNames);
				
my %AllPolymorphicLoci = ();
				
foreach my $name (@LocusNames)  {
	$AllPolymorphicLoci{$name}=1;
}
				
#Get total number of polymorphic loci, which is also the number of unlinked SNPs.
my $NumPolymorphicLoci = keys %AllPolymorphicLoci;
				
#Get total number of SNPs (includes linked SNPs).
my $TotalNumSNPs = @LocusNames;
				
close FILE;
close OUTFILE;
				
				
				
				
				
				
				
#Remove quotes in file.
				
#my $LineCounter = 0;
				
open SNPFILEWITHNA, "TempFiles/SNPMatrix_Edit.txt" or die$!;
open SNPFILENONA, ">TempFiles/SingleSNPsAllRaw.txt" or die$!;
				
while (<SNPFILEWITHNA>)  {
	#chomp($_);	
	$_ =~ s/"//g;
	print SNPFILENONA "$_";
	#if ($LineCounter == 0)  {
	#	$_ =~ s/"//g;
	#	$LineCounter++;
	#	print SNPFILENONA "$_";
	#}
						
	#else {
	#	$_ =~ s/"//g;
	#	my @TempArray = split(/\t/, $_);
	#	my $Length = @TempArray;
	#	print SNPFILENONA "$TempArray[0]\t";
							
	#	foreach my $allele (@TempArray[1..$Length-1])  {
	#		$allele =~ s/NA/N/g;
	#		print SNPFILENONA "$allele\t";
	#	}
	#}
						
	#print SNPFILENONA "\n";
}
				
close SNPFILEWITHNA;
close SNPFILENONA;
				
				
				
#Each line in SingleSNPsAllRaw.txt ends with \t\n and has quotes.  Remove these.
				
open SNPFILE, "TempFiles/SingleSNPsAllRaw.txt" or die$!;
open SNPFILEUPDATE, ">TempFiles/SingleSNPsAll.txt" or die$!;
				 
while (<SNPFILE>)  {
	$_ =~ s/\t\n$/\n/;
	$_ =~ s/"//g;
	print SNPFILEUPDATE "$_";
}
				 
close SNPFILE;
close SNPFILEUPDATE;
				

#Update the appropriate R script with the population sizes.

if ($Linked == 1)  {

	open IN, "../RScripts/OutputTreeMix_Linked.R" or die$!;
	open OUT, ">../RScripts/OutputTreeMix_Linked_Edit.R" or die$!;
	
	while(<IN>)  {
		if ($_ =~ /^PopSizesVector</) {
			print OUT "PopSizesVector<-c(";
			for my $popnumber (0..$NumPops-2) {
				print OUT "$DiploidSizes[$popnumber],";
			}
			print OUT "$DiploidSizes[$NumPops-1])\n";
		}

		else {
			print OUT "$_";
		}
	}

	close IN;
	close OUT;
}

else {
	open IN, "../RScripts/OutputTreeMix_Unlinked.R" or die$!;
	open OUT, ">../RScripts/OutputTreeMix_Unlinked_Edit.R" or die$!;
	
	while(<IN>)  {
		if ($_ =~ /^PopSizesVector</) {
			print OUT "PopSizesVector<-c(";
			for my $popnumber (0..$NumPops-2) {
				print OUT "$DiploidSizes[$popnumber],";
			}
			print OUT "$DiploidSizes[$NumPops-1])\n";
		}

		else {
			print OUT "$_";
		}
	}

	close IN;
	close OUT;
	
	
	
}	

				
#Run appropriate R script "OutputBiallelicSingleSNPs".  This outputs two files: AllBiallelicSNPsRaw.txt and UnlinkedBiallelicSNPs_Raw.txt. A maximum of one SNP is output for each locus in UnlinkedBiallelicSNPs_Raw.txt.

if ($Linked == 1)  {
	system "R --vanilla --slave < ../RScripts/OutputTreeMix_Linked_Edit.R";
}

else {
	system "R --vanilla --slave < ../RScripts/OutputTreeMix_Unlinked_Edit.R";
}

