#! /usr/bin/perl
use warnings;
use strict;


#This script uses a SNPMatrix file located in the Output/Genotypes directory to create a file that can be input into Beauti, which will in turn create a SNAPP input file.
#If you want to include a subset of samples from the SNPMatix file, set the subset argument to '1', and add a text file to the Formatting directory that contains the names of the samples to include.
#The default name of this text file is IncludedSamples.txt, but can be changed with the 'file' argument.



#Set parameters for the run

# -subset  Flag indicating whether to include a subset of the samples in the SNPMatrix file. Default is 0 (all samples are included in the Beauti input file)
# -file  Name of the file containing the samples to include in the dataset.  Default is 'IncludedSamples.txt'.  Only applies if subset (above) is set to '1'.
# -SNPsOnly Flag indicating whether to print only variable sites. Default is '0', which prints all sites (including monomorphic sites) to the Beauti infile.

my %RunArguments = ();

#Defaults
$RunArguments{subset} = 0;
$RunArguments{file} = 'IncludedSamples.txt';
$RunArguments{SNPsOnly} = 0;
$RunArguments{unlinked} = 1;


for my $argument (@ARGV)  {
	my @TempArgs = split(/-/,$argument);
	
	#print help information
	if (($TempArgs[1] =~ /^h$/)|| ($TempArgs[1] =~ /^help$/)) {
		
		print "Information for OutputSNAPP.pl...\n\n";
		print "This script uses a SNPMatrix file located in the Output/Genotypes directory to create a file that can be input into Beauti, which will in turn create a SNAPP input file.\n";  
		
		print "\n\nCommand line arguments available for this script...\n\n";
		print "subset\nFlag indicating whether to include a subset of the samples in the SNPMatrix file.\n";
		print "Set this to 1 to include a subset of the samples.\n";
		print "Default is 0 (all samples are included in the Beauti input file)\n";
		print "file\nName of the file containing the samples to include in the dataset.\n";
		print "Only applies if subset (above) is set to '1'.\n";
		print "Default is 'IncludedSamples.txt'.\n";
		print "\nSNPsOnly\nSet this to '1' to print only variable sites to the SNAPP input file.\n";
		print "Default is '0', which prints all sites\n";
		print "unlinked\nFlag indicating whether to include only unlinked SNPs in the SNAPP input file.\n";
		print "Default is 1 (only unlinked SNPs are included)\n\n";
		
		exit;
		
	}
	
	#update default values with entered parameters as appropriate and continue with run
	elsif (exists($RunArguments{$TempArgs[0]}))  {
		$RunArguments{$TempArgs[0]} = $TempArgs[1];
	}
}





my $AllOrSubset = $RunArguments{subset};
my $SampleFile = $RunArguments{file};
my $Unlinked = $RunArguments{unlinked};


my $FileName;


opendir GENOS, "../Output/Genotypes/";
my @AllFileNames = grep { $_ ne '.' && $_ ne '..' && $_ ne '.DS_Store' } readdir(GENOS);
close GENOS;

my @SNPFileNames = ();

for my $name (@AllFileNames)  {
	if ($name =~ /SNP/) {
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
	print "\nEnter the name of the SNPMatrix file you want to use to create the SNAPP infile\n";

	$FileName = <STDIN>;
	chomp($FileName);
}	





#Clean up the names in SNPMatrix file.


open FILE, "../Output/Genotypes/$FileName" or die$!;
mkdir "TempFiles" unless (-d "TempFiles");
open OUTFILE, ">TempFiles/SNPMatrix_Edit.txt" or die$!;

while(<FILE>)  {
	$_ =~ s/Individual//;
	print OUTFILE "$_";
}

close FILE;
close OUTFILE;



##############################################################################################################################################
#Subset the SNPMatrix file as necessary to include the desired samples.


if ($AllOrSubset == 1) {
	
	#Get the names of the individuals that will be included and store them as keys in hash SamplesWanted.

	my @SamplesWanted = ();

	open SAMPLESWANTED, "$SampleFile" or die$!;

	while (<SAMPLESWANTED>)  {
		chomp($_);
		push (@SamplesWanted, $_);
	}

	print "Samples included in SNAPP infile...@SamplesWanted\n\n";

	close SAMPLESWANTED;

	my %SamplesWanted;

	foreach my $sample (@SamplesWanted)  {
		$SamplesWanted{$sample} = 1;
	}

	
	#Go through SNPMatrix file and print the samples in the hash to file "General_Biallelics_Infile.txt";
	open SNAPPFILE, ">TempFiles/General_Biallelics_Infile.txt" or die$!;
	open ALLSAMPLESFILE, "TempFiles/SNPMatrix_Edit.txt" or die$!;

	my $Counter = 0;

	while (<ALLSAMPLESFILE>)  {
	
		if ($Counter == 0)  {
			print SNAPPFILE "$_";
			$Counter++;
		
		}
	
		else {
			my @TempArray = split(/\t/, $_);
			if (exists($SamplesWanted{$TempArray[0]}))  {
				print SNAPPFILE "$_";
			}
		}		
	}

	close SNAPPFILE;
	close ALLSAMPLESFILE;

}
	

else {
	system "cp TempFiles/SNPMatrix_Edit.txt TempFiles/General_Biallelics_Infile.txt";

	
}






#Replace any NA's in file SNAPP_Infile_NA.txt with "N".

my $SnappCounter = 0;

my $NumPolymorphicLoci;
my $TotalNumSNPs;

open SNAPPFILEWITHNA, "TempFiles/General_Biallelics_Infile.txt" or die$!;
open SNAPPFILENONA, ">TempFiles/SingleSNPsAllRaw.txt" or die$!;

	while (<SNAPPFILEWITHNA>)  {
		chomp($_);	

		if ($SnappCounter == 0)  {		#On first line with locus names - use hash to get number of unique loci
			$_ =~ s/"//g;
			$SnappCounter++;
			print SNAPPFILENONA "$_";
			
			my @LocusNamesArray = split(/\t/, $_);
			$TotalNumSNPs = @LocusNamesArray;
			my %PolymorphicLociHash = ();
			foreach my $locus (@LocusNamesArray)  {
				if ($locus =~ /[A-Za-z0-9]/) {
					$PolymorphicLociHash{$locus} = 1;
				}
			}
			
			$NumPolymorphicLoci = keys(%PolymorphicLociHash);
			
			
		}
		
		else {
			$_ =~ s/"//g;
			my @TempArray = split(/\t/, $_);
			my $Length = @TempArray;
			print SNAPPFILENONA "$TempArray[0]\t";
			
			foreach my $allele (@TempArray[1..$Length-1])  {
				$allele =~ s/NA/N/g;
				print SNAPPFILENONA "$allele\t";
			}
		}
		
		print SNAPPFILENONA "\n";
	}

close SNAPPFILEWITHNA;
close SNAPPFILENONA;

print "Number of polymorphic loci is $NumPolymorphicLoci\n\n";




#Each line in SingleSNPsRaw ends with \t\n and has quotes.  Remove these.

open SNPFILE, "TempFiles/SingleSNPsAllRaw.txt" or die$!;
open SNPFILEUPDATE, ">TempFiles/SingleSNPsAll.txt" or die$!;
 
while (<SNPFILE>)  {
	$_ =~ s/\t\n$/\n/;
	$_ =~ s/"//g;
	print SNPFILEUPDATE "$_";
}
 
close SNPFILE;
close SNPFILEUPDATE;




#Run R script "OutputBiallelicSingleSNPs".  This outputs file "UnlinkedBiallelicSNPs_Raw.txt". A maximum of one SNP is output for each locus. SNAPP assumes SNPs are unlinked.
	
system "R --vanilla --slave < ../RScripts/OutputBiallelicSingleSNPs.R";








#Each line in UnlinkedBiallelicSNPs_Raw.txt ends with \t\n and has quotes.  Remove these.

open SNPFILE, "TempFiles/UnlinkedBiallelicSNPsRaw.txt" or die$!;
open SNPFILEUPDATE, ">TempFiles/FinalBiallelicSNPs.txt" or die$!;


if ($Unlinked == 1)  {
	open SNPFILE, "TempFiles/UnlinkedBiallelicSNPsRaw.txt" or die$!;
}
				
else {
	open SNPFILE, "TempFiles/AllBiallelicSNPsRaw.txt" or die$!;
}


while (<SNPFILE>)  {
	$_ =~ s/\t\n$/\n/;
	$_ =~ s/"//g;
	print SNPFILEUPDATE "$_";
}
 
close SNPFILE;
close SNPFILEUPDATE;


	




#Run R Script OutputSNAPPMatrix.R.  This will convert the individuals to a single row with scores 0,1, or 2.

system "R --vanilla --slave < ../RScripts/OutputSNAPPMatrix.R";






open RAWMATRIX, "TempFiles/SNAPPMatrixRaw.txt" or die$!;
open SNAPPMATRIX, ">TempFiles/SNAPPMatrix.txt" or die$!;
 
while (<RAWMATRIX>)  {
	 
	$_ =~ s/"//g;
	print SNAPPMATRIX "$_";
}


#system "rm TempFiles/SNAPPMatrixRaw.txt";

close RAWMATRIX;
close SNAPPMATRIX;

my $TotalMonomorphicSites;
my $RetainedReadLength;

if ($RunArguments{SNPsOnly} == 0) {	#Add monomorphic sites to the file
	
	#Get number of monomorphic loci
	#First, get the percent of loci scored
	
	my @FirstSplit = split(/_/, $FileName);
	my $SecondElement = $FirstSplit[1];
	my @SecondSplit = split(/\./, $SecondElement);
	my $PctLociScored = $SecondSplit[0];
	$RetainedReadLength = $SecondSplit[1];
	
	if ($RetainedReadLength =~ /All/) {
		open MONOMORPHICINALL, "../Output/Genotypes/Monomorphics_$PctLociScored.txt" or die$!;
		while(<MONOMORPHICINALL>)  {
			if ($_ =~ /Sample/) {
				next;
			}
			
			else {
				if ($_ =~ /[ATGC]/) {
					chomp($_);
					$RetainedReadLength = length($_);
					last;
				}
			}
		}
		close MONOMORPHICINALL;
	}
	
	#print "PctLociScored is $PctLociScored\n\n";
	
	open MONOMORPHICINALL, "../Output/Genotypes/Monomorphics_$PctLociScored.txt" or die$!;	
	
	my @Monomorphics = ();

	while(<MONOMORPHICINALL>)  {
	
		if ($_ =~ /Sequence/) {
			next;
		}	
		
		if ($_ =~ /[ATGC]/)  {
			my @TempArray = split (/\s/, $_);
			my $Seq = $TempArray[0];
			push(@Monomorphics, $Seq);
			
		}
	}	

	my $NumMonomorphicLoci = @Monomorphics;
	my $MonomorphicSites = 	$NumMonomorphicLoci*$RetainedReadLength;
	
	my $TotalLengthOfPolymorphicLoci = $NumPolymorphicLoci*$RetainedReadLength;
	my $NumMonomorphicSitesInPolymorphicLoci = $TotalLengthOfPolymorphicLoci-$TotalNumSNPs;
	
	$TotalMonomorphicSites = $MonomorphicSites+$NumMonomorphicSitesInPolymorphicLoci;
	
	
}

 
 
open MATRIX, "TempFiles/SNAPPMatrix.txt" or die$!;
open NEXUSRAW, ">TempFiles/NexusRaw.txt" or die$!;
 
my $LineCounter2 = 0;
my @Allele1Array;
my @Allele2Array;
my $ArrayLength;
my $NumTaxa = 0;
 
while(<MATRIX>)  {
 
	#Skip the line with the locus names.
	
	if ($LineCounter2 == 0)  {
		$LineCounter2++;
		next;
	}

	#Store the first allele at each locus for the current individual in array "Allele1Array"
	elsif ($LineCounter2 == 1)  {
		 
		@Allele1Array = split (/\t/, $_);
		$LineCounter2 = 2;
		$ArrayLength = @Allele1Array;
		next;
	}
	 
	else {
		 
		#Store the second allele at each locus for the current individual in array "Allele2Array".
		@Allele2Array = split (/\t/, $_);
		$NumTaxa++;
		 
		
		my $LocusElement = 1;
		 
		print NEXUSRAW "$Allele1Array[0]\t";
		 
		foreach my $locus (@Allele1Array[1..$ArrayLength-1])  {
			 
			if ($locus =~ /N/)  {
				print NEXUSRAW "?";
			}	
			
			elsif (($locus == 0) && ($Allele2Array[$LocusElement] == 0))  {
				print NEXUSRAW "0";
			}
			 
			elsif (($locus == 1) && ($Allele2Array[$LocusElement] == 0))  {
				print NEXUSRAW "1";
			}
			 
			elsif (($locus == 0) && ($Allele2Array[$LocusElement] == 1))  {
				print NEXUSRAW "1";
			}
			 
			else {
				print NEXUSRAW "2";
			}
			 
			$LocusElement++;
		}
		 
		print NEXUSRAW "\n";
		 
		@Allele1Array = ();
		@Allele2Array = ();
			 
		$LineCounter2 = 1;
	}
}
 
close MATRIX;
close NEXUSRAW;


#Edit raw nexus file and output a file ready to open in Beauti.

open NEXUSRAW, "TempFiles/NexusRaw.txt" or die$!;
open NEXUS, ">Beauti_Infile.nex" or die$!;
 
my $NumChar = $ArrayLength-1;

if ($RunArguments{SNPsOnly} == 0)  {
	$NumChar = $NumChar+$TotalMonomorphicSites;
}	
 
print NEXUS "#NEXUS\n\n";
print NEXUS "begin data;\n";
 
print NEXUS "dimensions\tntax=$NumTaxa\tnchar=$NumChar;\n";
print NEXUS "format\tdatatype=integerdata\tmissing=\"?\"\tsymbols=\"012\";\n";
print NEXUS "matrix\n";
 
while (<NEXUSRAW>)  {
	if ($RunArguments{SNPsOnly} == 0)  {
		chomp($_);
		print NEXUS "$_";
		my $TempString = 0 x $TotalMonomorphicSites;
		print NEXUS "$TempString\n";
	}	
	
	else {
		print NEXUS "$_";
	}	
	
}
 
print NEXUS ";\n";
print NEXUS "End;";
 
close NEXUS;
close NEXUSRAW;


system "rm TempFiles/*";
system "rmdir TempFiles";
