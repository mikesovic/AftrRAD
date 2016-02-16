#! /usr/bin/perl
use warnings;
use strict;

#Set parameters for the run

# -resamp Flag indicating whether to create resampled datasets.
# -num  Number of resampled datasets to create (default 1).
# -prop Percent of loci sampled to create the resampled datasets (1-99).
# -rep  Flag indicating whether to sample with replacement when generating resampled datasets.  Default is 0 (sample loci without replacement).
# -unlinked Flag indicating whether to include only unlinked SNPs in the site frequency spectrum.  Default is 0 (all SNPs are included).
# -folded  Flag indicating the site frequency spectrum will be folded.  Default is 0 (unfolded spectrum, which requires an outgroup sample).


my %RunArguments = ();

#Defaults
$RunArguments{resamp} = 0;
$RunArguments{num} = 1;
$RunArguments{pct} = 50;
$RunArguments{rep} = 0;
$RunArguments{unlinked} = 0;
$RunArguments{folded} = 0;
$RunArguments{MonoScaled} = 0;


for my $argument (@ARGV)  {
	my @TempArgs = split(/-/,$argument);
	
	#print help information
	if (($TempArgs[1] =~ /^h$/)|| ($TempArgs[1] =~ /^help$/)) {
		
		print "Information for OutputFastSimCoal_SingleSFS.pl...\n\n";
		print "This script creates an unfolded site frequency spectrum from a single population that can be used as input for analyses in FastSimCoal\n";  
		print "It requires a SNPMatrix file with no missing data (SNPMatrix_100.X.txt) in the Output/Genotypes folder, and also requires that an outgroup sample is included as the first sample in this SNPMatrix file.\n";
		
		print "\n\nCommand line arguments available for this script...\n\n";
		print "unlinked\nFlag indicating whether to include only unlinked SNPs in the site frequency spectrum.\n";
		print "Default is 0 (all SNPs are included)\n\n";
		print "resamp\nFlag indicating whether to create resampled datasets.\n";
		print "Set this to 1 to perform resampling.\n\n";
		print "num\nNumber of resampled datasets to create (Default is 1).\n\n";
		
		print "pct\nPercent of loci sampled to create each resampled dataset (1-99, default is 50).\n\n";
		
		print "rep\nFlag indicating whether to sample with replacement when generating resampled datasets.  Default is 0 (sample loci without replacement).\n\n";
		print "folded\nFlag indicating the site frequency spectrum will be folded.\n";
		print "Default is 0 (unfolded spectrum, which requires that an outgroup sample is the first sample in the SNPMatrix file).\n\n";
		print "MonoScaled\nFlag indicating whether to scale the number of monomorphic sites based on the proportion of SNPs retained after removing linked SNPs.\n";
		print "Default is 0, meaning all monomorphic sites are counted.  If set to '1', the proportion of total SNPs that are unlinked is calculated, and the total\n";
		print "number of monomorphic sites is scaled by this proportion.\n\n";
		
		
		exit;
		
	}
	
	#update default values with entered parameters as appropriate and continue with run
	elsif (exists($RunArguments{$TempArgs[0]}))  {
		$RunArguments{$TempArgs[0]} = $TempArgs[1];
	}
}





my $resamp = $RunArguments{resamp};
my $NumResampReps = $RunArguments{num};
my $PctLoci = $RunArguments{pct};
my $Replacement = $RunArguments{rep};
my $Unlinked = $RunArguments{unlinked};
my $Folded = $RunArguments{folded};
my $ScaleMonos = $RunArguments{MonoScaled};

#Get the name of the SNPMatrix file to use to create the SFS.  If there is only one available, that one is used.  
#If there is more than one in the Output/Genotypes directory, the user is prompted to enter the name of the one to use.

my $FileName;


opendir GENOS, "../Output/Genotypes/";
my @AllFileNames = grep { $_ ne '.' && $_ ne '..' && $_ ne '.DS_Store' } readdir(GENOS);
close GENOS;

my @SNPFileNames = ();

for my $name (@AllFileNames)  {
	if ($name =~ /SNPMatrix_100/) {
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
	print "\nEnter the name of the SNPMatrix file you want to use to create the FastSimCoal infile\n";

	$FileName = <STDIN>;
	chomp($FileName);
}




if ($Folded == 1)  {
	
	
	#Get the number of ingroup samples in the dataset
	
	open FILE, "../Output/Genotypes/$FileName" or die$!;
	
	my $NumLinesInFile = 0;
	
	while(<FILE>)  {
		if (($_ =~ /[A-Za-z]/) && ($_ !~ /Locus/)) {
			$NumLinesInFile++;
		}
	}
	
	close FILE;
	
	my $PopSizes = ($NumLinesInFile)/2;
	
	
	mkdir "TempFiles" unless (-d "TempFiles");
	
	my @HapSizes = split(/\t/,$PopSizes);
	
	my @DiploidSizes = ();
	
	foreach my $size (@HapSizes)  {
		my $TempSize = $size*2;
		push (@DiploidSizes,$TempSize);
	}
	
	my $NumPops = @DiploidSizes;	
		
	
	
	print "\n\n";
	
	
	
	
	##############################################################################################################################################
	##############################################################################################################################################
	#If the resamp argument is set to 1, need to generate the resampled datasets.
	
	if ($resamp == 1)  {
	
		print "Arguments entered are...\n";
		
		for (keys %RunArguments) {
			print "$_\t$RunArguments{$_}\n";
		}
		
		
		mkdir "ResampledDatasets" unless (-d "ResampledDatasets");
		mkdir "ResampledSFS" unless (-d "ResampledSFS");
		
		#Get the number of loci in the SNPMatrix file.
		my $TotalNumSNPs;
		
		open SNPMATRIX, "../Output/Genotypes/$FileName" or die$!;
		
		while(<SNPMATRIX>)  {
			my @LocusNames = split(/\t/,$_);
			$TotalNumSNPs = @LocusNames;
			$TotalNumSNPs = $TotalNumSNPs-1;
			last;
		}
	
		close SNPMATRIX;
		
		my $NumSNPsToSample = int(($PctLoci/100)*$TotalNumSNPs);
		
		
		#Edit Resample.R script with updated SNPMatrix name and number of loci to sample.
		open RFILE, "../RScripts/Resample.R" or die$!;
		open ROUT, ">../RScripts/Resample_Edit.R" or die$!;
		
		while(<RFILE>) {
			if (($_ =~ /read/) && ($_ =~ /SNPMatrix/))  {
				print ROUT "DataTableFromFile<-read.table(\"../Output/Genotypes/$FileName\", header=TRUE, sep=\"\t\")\n";
			}
			
			elsif ($_ =~ /NumLociToSample</)  {
				print ROUT "NumLociToSample<-$NumSNPsToSample\n";
					
			}
			
			elsif ($_ =~ /columns<-/)  {
				if ($Replacement == 1)  {
					print ROUT "columns<-c(sample(2:NumColsInMatrix, NumLociToSample, replace=T))\n";
				}
				
				else {
					print ROUT "$_";
				}
			}
			
			
			else {
				print ROUT "$_";
			}
		}
	
	
	
		close RFILE;
		close ROUT;
		
		
		
		for my $RepDataset (1..$NumResampReps)  {
		
			#Create each resampled dataset
	
			print "Creating resampled dataset $RepDataset of $NumResampReps\n";
			system "R --vanilla --slave < ../RScripts/Resample_Edit.R";
	
			my $TempFileName = $FileName;
			$TempFileName =~ s/.txt//;
			system "mv ResampledDatasets/SNPMatrix_resamp.txt ResampledDatasets/$TempFileName.$RepDataset.txt";
			
			
			
	
			open FILE, "ResampledDatasets/$TempFileName.$RepDataset.txt" or die$!;
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
	
			my $NumPolymorphicLoci = keys %AllPolymorphicLoci;
			my $TotalNumSNPs = @LocusNames;
	
			close FILE;
			close OUTFILE;
	
	
	
	
	
	
	
			#Replace any NA's in file SNPMatrix_Edit.txt with "N".
	
			my $LineCounter = 0;
	
			open SNPFILEWITHNA, "TempFiles/SNPMatrix_Edit.txt" or die$!;
			open SNPFILENONA, ">TempFiles/SingleSNPsAllRaw.txt" or die$!;
	
			while (<SNPFILEWITHNA>)  {
				chomp($_);	
	
				if ($LineCounter == 0)  {
					$_ =~ s/"//g;
					$LineCounter++;
					print SNPFILENONA "$_";
				}
			
				else {
					$_ =~ s/"//g;
					my @TempArray = split(/\t/, $_);
					my $Length = @TempArray;
					print SNPFILENONA "$TempArray[0]\t";
				
					foreach my $allele (@TempArray[1..$Length-1])  {
						$allele =~ s/NA/N/g;
						print SNPFILENONA "$allele\t";
					}
				}
			
				print SNPFILENONA "\n";
			}
	
			close SNPFILEWITHNA;
			close SNPFILENONA;
	
	
	
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
	
	
			#Run R script "OutputBiallelicSingleSNPs".  This outputs file "UnlinkedBiallelicSNPs_Raw.txt". A maximum of one SNP is output for each locus.
		
			system "R --vanilla --slave < ../RScripts/OutputBiallelicSingleSNPs.R";
	
	
			my $NumUnlinkedBiallelicLoci;
	
	
	
			#Each line in UnlinkedBiallelicSNPsRaw.txt ends with \t\n and has quotes.  Remove these.
	
			if ($Unlinked == 1)  {
				open SNPFILE, "TempFiles/UnlinkedBiallelicSNPsRaw.txt" or die$!;
			}
		
			else {
				open SNPFILE, "TempFiles/AllBiallelicSNPsRaw.txt" or die$!;
			}	
			
			open SNPFILEUPDATE, ">TempFiles/BiallelicSNPs_SFS.txt" or die$!;
	 
			while (<SNPFILE>)  {
		
				if ($_ =~ /[A-Za-z1-9]/)  {
					$_ =~ s/\t\n$/\n/;
					$_ =~ s/"//g;
					my @TempArray = split(/\t/,$_);
					$NumUnlinkedBiallelicLoci = @TempArray;	#has an empty tab at beginning
					$NumUnlinkedBiallelicLoci = $NumUnlinkedBiallelicLoci-1;
					print SNPFILEUPDATE "$_";
				}	
			}
	 
			close SNPFILE;
			close SNPFILEUPDATE;
	
	
	##############################################################################################################################################
	
	
			print "\n\nNumber of biallelic loci for site frequency spectrum is $NumUnlinkedBiallelicLoci\n\n";
	
			my $TotalMonomorphicSites;
			my $RetainedReadLength;
			
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
			my $ScaledMonomorphicLoci = int(($PctLoci/100)*$NumMonomorphicLoci);   #If resampled dataset includes less than 100% of the loci, need to scale the monomorphic sites.
					
			my $MonomorphicSites = 	$ScaledMonomorphicLoci*$RetainedReadLength;
					
			my $TotalLengthOfPolymorphicLoci = $NumPolymorphicLoci*$RetainedReadLength;
			my $NumMonomorphicSitesInPolymorphicLoci = $TotalLengthOfPolymorphicLoci-$TotalNumSNPs;
				
			$TotalMonomorphicSites = $MonomorphicSites+$NumMonomorphicSitesInPolymorphicLoci;
	
			if ($ScaleMonos == 1)  {
				my $ProportionUnlinkedSNPs = $NumUnlinkedBiallelicLoci/$TotalNumSNPs;
				$TotalMonomorphicSites = int($TotalMonomorphicSites*$ProportionUnlinkedSNPs);
					
			}
			
			##############################################################################################################################################
	
			#Create hash with possible values of the derived allele counts as keys (keys 0 to Number of alleles sampled).
	
			my %DerivedCountsHash = ();
	
			foreach my $value (0..$DiploidSizes[0])  {
		
				$DerivedCountsHash{$value} = 0;
			}
	
	
			
			##############################################################################################################################################
			#Go through file UnlinkedBiallelicSNPs.txt and count ancestral and derived alleles at each locus.
			
			open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
	
			my @UnlinkedBiallelicsNames = ();
	
			while(<TOTALBIALLELICS>) {
				@UnlinkedBiallelicsNames = split(/\t/, $_);
				shift(@UnlinkedBiallelicsNames);
				last;
			}
	
			close TOTALBIALLELICS;
	
			
			
			#First, define ancestral and derived alleles at each locus, as the spectrum is folded.
	
			foreach my $locusname (@UnlinkedBiallelicsNames) {
	
				my $CurrentSiteInTotalBiallelics = 0;
				my @CurrentSiteIndividualIDs = ();
				my @CurrentSiteAlleles = ();
				my $CurrentFirstAllele;
				my $CurrentSecondAllele;
				my $NumFirstAlleles = 1;
				my $NumSecondAlleles = 0;

		
				open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
			
				my $TempCounter = 0;
			
					
				#First, identify the minor allele at the locus
				
				while (<TOTALBIALLELICS>)  {	#get the column number for that locus in GenotypesUpdate.txt
				
				
					if ($TempCounter == 0)  {		#on the first line that has locus names.
						my @TempArray = split(/\t/,$_);
						my $TempArrayElementCounter = 0;
				
						foreach my $name (@TempArray)  {
							if ($locusname eq $name) {	#find the column number of the locus we are currently on
								$CurrentSiteInTotalBiallelics = $TempArrayElementCounter;
								last;
							}	
						
							else {
								$TempArrayElementCounter++;
							}
						
						}
						$TempCounter++;
					
					}
				
					elsif ($TempCounter == 1)  {	#On the first sample - will use this allele as CurrentFirstAllele.
							
						my @TempArray = split(/\t/,$_);
						$CurrentFirstAllele = $TempArray[$CurrentSiteInTotalBiallelics];
					
						$TempCounter++;
					
					}
				
				
				
				
					else {
				
						my @TempArray = split(/\t/,$_);
						my $CurrentRead = $TempArray[$CurrentSiteInTotalBiallelics];
						if ($CurrentRead eq $CurrentFirstAllele)  {
							$NumFirstAlleles++;
						}
							
						else {
							$CurrentSecondAllele = $CurrentRead;
							$NumSecondAlleles++;
						}	
							
							
						
					}
				}
			
				close TOTALBIALLELICS;
				
				
				my @AncestralAllelePopCounts = ();
				my @DerivedAllelePopCounts = ();
				my $EqualFreqFlag = 0;
	
				#Define the derived allele and add the count to the DerivedCountsHash hash
				
				my $CurrentAncestralAllele;
				my $CurrentDerivedAllele;
				
				
				if ($NumFirstAlleles > $NumSecondAlleles)  {
					$CurrentDerivedAllele = $CurrentSecondAllele;
					$DerivedCountsHash{$NumSecondAlleles}++;
					
				}
				
				elsif ($NumSecondAlleles > $NumFirstAlleles)  {
					$CurrentDerivedAllele = $CurrentFirstAllele;
					$DerivedCountsHash{$NumFirstAlleles}++;
				}
				
				elsif ($NumSecondAlleles == $NumFirstAlleles) {
					$DerivedCountsHash{$NumFirstAlleles}++;
				}	
			}
		
		
		
			#Create the SFS	
				
			open OUTFILE, ">ResampledSFS/$RepDataset.resamp_MAFpop0.obs" or die$!;		#Create file for current population to store DFS.
			
			print OUTFILE "1 observations\n";     #\t$NumMonomorphicLoci monomorphic loci scored in all\t$NumPolymorphicLoci polymorphic loci scored in all\t$ReadLength bases in each read\n";
			
			foreach my $Pop1Count (0..$DiploidSizes[0]-1)  {
				print OUTFILE "d0\_$Pop1Count\t";
				
			}
			
			print OUTFILE "d0\_$DiploidSizes[0]\n";
			
			print OUTFILE "$TotalMonomorphicSites\t";
			
			foreach my $Pop1Count (1..$DiploidSizes[0]-1)  {
				print OUTFILE "$DerivedCountsHash{$Pop1Count}\t";
			}
			
			print OUTFILE "$DerivedCountsHash{$DiploidSizes[0]}";
			
			
			system "rm TempFiles/*";
				
			
		
		}
		system "rmdir TempFiles";
		print "Resampled folded site frequency spectra have been generated and printed to folder ResampledSFS in Formatting directory.\n";
	
	}
	
	
	
	
	
	
	
	else {		#folded spectrum without resampling
		
		
		
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
	
			my $NumPolymorphicLoci = keys %AllPolymorphicLoci;
			my $TotalNumSNPs = @LocusNames;
	
			close FILE;
			close OUTFILE;
	
	
	
	
	
	
	
			#Replace any NA's in file SNPMatrix_Edit.txt with "N".
	
			my $LineCounter = 0;
	
			open SNPFILEWITHNA, "TempFiles/SNPMatrix_Edit.txt" or die$!;
			open SNPFILENONA, ">TempFiles/SingleSNPsAllRaw.txt" or die$!;
	
			while (<SNPFILEWITHNA>)  {
				chomp($_);	
	
				if ($LineCounter == 0)  {
					$_ =~ s/"//g;
					$LineCounter++;
					print SNPFILENONA "$_";
				}
			
				else {
					$_ =~ s/"//g;
					my @TempArray = split(/\t/, $_);
					my $Length = @TempArray;
					print SNPFILENONA "$TempArray[0]\t";
				
					foreach my $allele (@TempArray[1..$Length-1])  {
						$allele =~ s/NA/N/g;
						print SNPFILENONA "$allele\t";
					}
				}
			
				print SNPFILENONA "\n";
			}
	
			close SNPFILEWITHNA;
			close SNPFILENONA;
	
	
	
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
	
	
			#Run R script "OutputBiallelicSingleSNPs".  This outputs file "UnlinkedBiallelicSNPs_Raw.txt". A maximum of one SNP is output for each locus.
		
			system "R --vanilla --slave < ../RScripts/OutputBiallelicSingleSNPs.R";
	
	
			my $NumUnlinkedBiallelicLoci;
	
	
	
			#Each line in UnlinkedBiallelicSNPsRaw.txt ends with \t\n and has quotes.  Remove these.
	
			if ($Unlinked == 1)  {
				open SNPFILE, "TempFiles/UnlinkedBiallelicSNPsRaw.txt" or die$!;
			}
		
			else {
				open SNPFILE, "TempFiles/AllBiallelicSNPsRaw.txt" or die$!;
			}	
			
			open SNPFILEUPDATE, ">TempFiles/BiallelicSNPs_SFS.txt" or die$!;
	 
			while (<SNPFILE>)  {
		
				if ($_ =~ /[A-Za-z1-9]/)  {
					$_ =~ s/\t\n$/\n/;
					$_ =~ s/"//g;
					my @TempArray = split(/\t/,$_);
					$NumUnlinkedBiallelicLoci = @TempArray;	#has an empty tab at beginning
					$NumUnlinkedBiallelicLoci = $NumUnlinkedBiallelicLoci-1;
					print SNPFILEUPDATE "$_";
				}	
			}
	 
			close SNPFILE;
			close SNPFILEUPDATE;
	
	
	##############################################################################################################################################
	
	
			print "\n\nNumber of biallelic loci for site frequency spectrum is $NumUnlinkedBiallelicLoci\n\n";
	
			my $TotalMonomorphicSites;
			my $RetainedReadLength;
			
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
				
			if ($ScaleMonos == 1)  {
				my $ProportionUnlinkedSNPs = $NumUnlinkedBiallelicLoci/$TotalNumSNPs;
				$TotalMonomorphicSites = int($TotalMonomorphicSites*$ProportionUnlinkedSNPs);
					
			}
	
			##############################################################################################################################################
	
			#Create hash with possible values of the derived allele counts as keys (keys 0 to Number of alleles sampled).
	
			my %DerivedCountsHash = ();
	
			foreach my $value (0..$DiploidSizes[0])  {
		
				$DerivedCountsHash{$value} = 0;
			}
	
	
			
			
			##############################################################################################################################################
			#Go through file UnlinkedBiallelicSNPs.txt and count ancestral and derived alleles at each locus.
			
			open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
	
			my @UnlinkedBiallelicsNames = ();
	
			while(<TOTALBIALLELICS>) {
				@UnlinkedBiallelicsNames = split(/\t/, $_);
				shift(@UnlinkedBiallelicsNames);
				last;
			}
	
			close TOTALBIALLELICS;
	
			
			
			#First, define ancestral and derived alleles at each locus, as the spectrum is folded.
	
			foreach my $locusname (@UnlinkedBiallelicsNames) {
	
				my $CurrentSiteInTotalBiallelics = 0;
				my @CurrentSiteIndividualIDs = ();
				my @CurrentSiteAlleles = ();
				my $CurrentFirstAllele;
				my $CurrentSecondAllele;
				my $NumFirstAlleles = 1;
				my $NumSecondAlleles = 0;

		
				open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
			
				my $TempCounter = 0;
			
					
				#First, identify the minor allele at the locus
				
				while (<TOTALBIALLELICS>)  {	#get the column number for that locus in GenotypesUpdate.txt
				
				
					if ($TempCounter == 0)  {		#on the first line that has locus names.
						my @TempArray = split(/\t/,$_);
						my $TempArrayElementCounter = 0;
				
						foreach my $name (@TempArray)  {
							if ($locusname eq $name) {	#find the column number of the locus we are currently on
								$CurrentSiteInTotalBiallelics = $TempArrayElementCounter;
								last;
							}	
						
							else {
								$TempArrayElementCounter++;
							}
						
						}
						$TempCounter++;
					
					}
				
					elsif ($TempCounter == 1)  {	#On the first sample - will use this allele as CurrentFirstAllele.
							
						my @TempArray = split(/\t/,$_);
						$CurrentFirstAllele = $TempArray[$CurrentSiteInTotalBiallelics];
					
						$TempCounter++;
					
					}
				
				
						
				
				
					else {
				
						my @TempArray = split(/\t/,$_);
						my $CurrentRead = $TempArray[$CurrentSiteInTotalBiallelics];
						if ($CurrentRead eq $CurrentFirstAllele)  {
							$NumFirstAlleles++;
						}
							
						else {
							$CurrentSecondAllele = $CurrentRead;
							$NumSecondAlleles++;
						}	
							
							
						
					}
				}
			
				close TOTALBIALLELICS;
				
				
				my @AncestralAllelePopCounts = ();
				my @DerivedAllelePopCounts = ();
				my $EqualFreqFlag = 0;
	
				#Define the derived allele and add the count to the DerivedCountsHash hash
				
				my $CurrentAncestralAllele;
				my $CurrentDerivedAllele;
				
				
				if ($NumFirstAlleles > $NumSecondAlleles)  {
					$CurrentDerivedAllele = $CurrentSecondAllele;
					$DerivedCountsHash{$NumSecondAlleles}++;
					
				}
				
				elsif ($NumSecondAlleles > $NumFirstAlleles)  {
					$CurrentDerivedAllele = $CurrentFirstAllele;
					$DerivedCountsHash{$NumFirstAlleles}++;
				}
				
				elsif ($NumSecondAlleles == $NumFirstAlleles) {
					$DerivedCountsHash{$NumFirstAlleles}++;
				}	
			}	
				
		
		#Create the SFS	
			
		open OUTFILE, ">_MAFpop0.obs" or die$!;		#Create file for current population to store DFS.
		
		print OUTFILE "1 observations\n";     #\t$NumMonomorphicLoci monomorphic loci scored in all\t$NumPolymorphicLoci polymorphic loci scored in all\t$ReadLength bases in each read\n";
		
		foreach my $Pop1Count (0..$DiploidSizes[0]-1)  {
			print OUTFILE "d0\_$Pop1Count\t";
			
		}
		
		print OUTFILE "d0\_$DiploidSizes[0]\n";
		
		print OUTFILE "$TotalMonomorphicSites\t";
		
		foreach my $Pop1Count (1..$DiploidSizes[0]-1)  {
			print OUTFILE "$DerivedCountsHash{$Pop1Count}\t";
		}
		
		print OUTFILE "$DerivedCountsHash{$DiploidSizes[0]}";
		
		
		system "rm TempFiles/*";
		system "rmdir TempFiles";			
					
		print "Folded site frequency spectrum has been generated and printed to Formatting directory with name _MAFpop0.obs.\n";			
			
		}
}
			
			
			
	











else {		#SFS is unfolded (outgroup is included)
	
	
	#Get the number of ingroup samples in the dataset
	
	open FILE, "../Output/Genotypes/$FileName" or die$!;
	
	my $NumLinesInFile = 0;
	
	while(<FILE>)  {
		if (($_ =~ /[A-Za-z]/) && ($_ !~ /Locus/)) {
			$NumLinesInFile++;
		}
	}
	
	close FILE;
	
	my $PopSizes = ($NumLinesInFile-2)/2;
	
	
	mkdir "TempFiles" unless (-d "TempFiles");
	
	my @HapSizes = split(/\t/,$PopSizes);
	
	my @DiploidSizes = ();
	
	foreach my $size (@HapSizes)  {
		my $TempSize = $size*2;
		push (@DiploidSizes,$TempSize);
	}
	
	my $NumPops = @DiploidSizes;
	
	if ($NumPops == 1)  {
		
		print "One population recognized\n";
	}
	
	else {
		print "Recognized $NumPops populations.  Because there is more than one population, run OutputFastSimCoal_JointSFS.pl";
		exit;
	}	
		
	
	
	print "\n\n";
		
	
	##############################################################################################################################################
	##############################################################################################################################################
	#If the resamp argument is set to 1, need to generate the resampled datasets.
	
	if ($resamp == 1)  {
	
		print "Arguments entered are...\n";
		
		for (keys %RunArguments) {
			print "$_\t$RunArguments{$_}\n";
		}
		
		
		mkdir "ResampledDatasets" unless (-d "ResampledDatasets");
		mkdir "ResampledSFS" unless (-d "ResampledSFS");
		
		#Get the number of loci in the SNPMatrix file.
		my $TotalNumSNPs;
		
		open SNPMATRIX, "../Output/Genotypes/$FileName" or die$!;
		
		while(<SNPMATRIX>)  {
			my @LocusNames = split(/\t/,$_);
			$TotalNumSNPs = @LocusNames;
			$TotalNumSNPs = $TotalNumSNPs-1;
			last;
		}
	
		close SNPMATRIX;
		
		my $NumSNPsToSample = int(($PctLoci/100)*$TotalNumSNPs);
		
		
		#Edit Resample.R script with updated SNPMatrix name and number of loci to sample.
		open RFILE, "../RScripts/Resample.R" or die$!;
		open ROUT, ">../RScripts/Resample_Edit.R" or die$!;
		
		while(<RFILE>) {
			if (($_ =~ /read/) && ($_ =~ /SNPMatrix/))  {
				print ROUT "DataTableFromFile<-read.table(\"../Output/Genotypes/$FileName\", header=TRUE, sep=\"\t\")\n";
			}
			
			elsif ($_ =~ /NumLociToSample</)  {
				print ROUT "NumLociToSample<-$NumSNPsToSample\n";
					
			}
			
			elsif ($_ =~ /columns<-/)  {
				if ($Replacement == 1)  {
					print ROUT "columns<-c(sample(2:NumColsInMatrix, NumLociToSample, replace=T))\n";
				}
				
				else {
					print ROUT "$_";
				}
			}
			
			
			else {
				print ROUT "$_";
			}
		}
	
	
	
		close RFILE;
		close ROUT;
		
		
		
		for my $RepDataset (1..$NumResampReps)  {
		
			#Create each resampled dataset
	
			print "Creating resampled dataset $RepDataset of $NumResampReps\n";
			system "R --vanilla --slave < ../RScripts/Resample_Edit.R";
	
			my $TempFileName = $FileName;
			$TempFileName =~ s/.txt//;
			system "mv ResampledDatasets/SNPMatrix_resamp.txt ResampledDatasets/$TempFileName.$RepDataset.txt";
			
			
			
	
			open FILE, "ResampledDatasets/$TempFileName.$RepDataset.txt" or die$!;
			open OUTFILE, ">TempFiles/SNPMatrix_Edit.txt" or die$!;
	
			my @LocusNames = ();
			my $LineNumber = 0;
	
			while(<FILE>)  {
		
				if ($LineNumber == 0)  {
					my @TempLocusNames = split(/\t/,$_);
					foreach my $locname(@TempLocusNames)  {
						$locname =~ s/\.[0-9]+//;
						push (@LocusNames, $locname);
					}	
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
	
			my $NumPolymorphicLoci = keys %AllPolymorphicLoci;
			my $TotalNumSNPs = @LocusNames;
	
			close FILE;
			close OUTFILE;
	
	
	
	
	
	
	
			#Replace any NA's in file SNAPP_Infile_NA.txt with "N".
	
			my $LineCounter = 0;
	
			open SNPFILEWITHNA, "TempFiles/SNPMatrix_Edit.txt" or die$!;
			open SNPFILENONA, ">TempFiles/SingleSNPsAllRaw.txt" or die$!;
	
			while (<SNPFILEWITHNA>)  {
				chomp($_);	
	
				if ($LineCounter == 0)  {
					$_ =~ s/"//g;
					$LineCounter++;
					print SNPFILENONA "$_";
				}
			
				else {
					$_ =~ s/"//g;
					my @TempArray = split(/\t/, $_);
					my $Length = @TempArray;
					print SNPFILENONA "$TempArray[0]\t";
				
					foreach my $allele (@TempArray[1..$Length-1])  {
						$allele =~ s/NA/N/g;
						print SNPFILENONA "$allele\t";
					}
				}
			
				print SNPFILENONA "\n";
			}
	
			close SNPFILEWITHNA;
			close SNPFILENONA;
	
	
	
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
	
	
	
	
			#Run R script "OutputBiallelicSingleSNPs".  This outputs file "UnlinkedBiallelicSNPs_Raw.txt". A maximum of one SNP is output for each locus.
		
			system "R --vanilla --slave < ../RScripts/OutputBiallelicSingleSNPs.R";
	
	
			my $NumUnlinkedBiallelicLoci;
	
	
	
			#Each line in UnlinkedBiallelicSNPsRaw.txt ends with \t\n and has quotes.  Remove these.
	
			if ($Unlinked == 1)  {
				open SNPFILE, "TempFiles/UnlinkedBiallelicSNPsRaw.txt" or die$!;
			}
		
			else {
				open SNPFILE, "TempFiles/AllBiallelicSNPsRaw.txt" or die$!;
			}	
			
			open SNPFILEUPDATE, ">TempFiles/BiallelicSNPs_SFS.txt" or die$!;
	 
			while (<SNPFILE>)  {
		
				if ($_ =~ /[A-Za-z1-9]/)  {
					$_ =~ s/\t\n$/\n/;
					$_ =~ s/"//g;
					my @TempArray = split(/\t/,$_);
					$NumUnlinkedBiallelicLoci = @TempArray;	#has an empty tab at beginning
					$NumUnlinkedBiallelicLoci = $NumUnlinkedBiallelicLoci-1;
					print SNPFILEUPDATE "$_";
				}	
			}
	 
			close SNPFILE;
			close SNPFILEUPDATE;
	
	
	
	##############################################################################################################################################
	
	
			print "\n\nNumber of biallelic loci for site frequency spectrum is $NumUnlinkedBiallelicLoci\n\n";
	
			my $TotalMonomorphicSites;
			my $RetainedReadLength;
			
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
			my $ScaledMonomorphicLoci = int(($PctLoci/100)*$NumMonomorphicLoci);   #If resampled dataset includes less than 100% of the loci, need to scale the monomorphic sites.
			my $MonomorphicSites = 	$ScaledMonomorphicLoci*$RetainedReadLength;
			
			
			my $TotalLengthOfPolymorphicLoci = $NumPolymorphicLoci*$RetainedReadLength;
			my $NumMonomorphicSitesInPolymorphicLoci = $TotalLengthOfPolymorphicLoci-$TotalNumSNPs;
				
			$TotalMonomorphicSites = $MonomorphicSites+$NumMonomorphicSitesInPolymorphicLoci;
				
			if ($ScaleMonos == 1)  {
				my $ProportionUnlinkedSNPs = $NumUnlinkedBiallelicLoci/$TotalNumSNPs;
				$TotalMonomorphicSites = int($TotalMonomorphicSites*$ProportionUnlinkedSNPs);
					
			}
	
	
	##############################################################################################################################################
	
			#Create hash with possible values of the derived allele counts as keys (keys 0 to Number of alleles sampled).
	
			
			my %DerivedCountsHash = ();
	
			foreach my $value (0..$DiploidSizes[0])  {
		
				$DerivedCountsHash{$value} = 0;
			}
	
	##############################################################################################################################################
			#Go through file UnlinkedBiallelicSNPs.txt and count ancestral and derived alleles at each locus.
	
	
			my @AncestralAllelePopCounts = ();
			my @DerivedAllelePopCounts = ();
	
	
			my $CurrentAncestralAllele;
			my $CurrentDerivedAllele;
				
			open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
	
			my @UnlinkedBiallelicsNames = ();
	
			while(<TOTALBIALLELICS>) {
				@UnlinkedBiallelicsNames = split(/\t/, $_);
				shift(@UnlinkedBiallelicsNames);
				last;
			}
	
			close TOTALBIALLELICS;
	
	
	
			foreach my $locusname (@UnlinkedBiallelicsNames) {
	
				my $CurrentSiteInTotalBiallelics = 0;
				my @CurrentSiteIndividualIDs = ();
				my @CurrentSiteAlleles = ();
		
				open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
			
				my $TempCounter = 0;
			
					while (<TOTALBIALLELICS>)  {	#get the column number for that locus in GenotypesUpdate.txt
				
				
						if ($TempCounter == 0)  {		#on the first line that has locus names.
							my @TempArray = split(/\t/,$_);
							my $TempArrayElementCounter = 0;
				
							foreach my $name (@TempArray)  {
								if ($locusname eq $name) {	#find the column number of the locus we are currently on
									$CurrentSiteInTotalBiallelics = $TempArrayElementCounter;
									last;
								}	
						
								else {
									$TempArrayElementCounter++;
								}
						
							}
							$TempCounter++;
					
						}
				
						elsif ($TempCounter == 1)  {	#Should be on the outgroup individual - this individual has to be placed in this first position in SNPMatrix_X.Y.txt
							
							my @TempArray = split(/\t/,$_);
							$CurrentAncestralAllele = $TempArray[$CurrentSiteInTotalBiallelics];
					
							$TempCounter++;
					
						}
				
				
						elsif ($TempCounter == 2)  {	#Second line of outgroup - ignore
							$TempCounter++;
					
						}
				
				
				
						else {
					
							my @TempArray = split(/\t/,$_);
							my $CurrentRead = $TempArray[$CurrentSiteInTotalBiallelics];
							push(@CurrentSiteIndividualIDs, $TempCounter);			#keep track of what sample we're working on
							push(@CurrentSiteAlleles, $TempArray[$CurrentSiteInTotalBiallelics]);
						
							if ($TempArray[$CurrentSiteInTotalBiallelics] !~ /$CurrentAncestralAllele/)	{
								$CurrentDerivedAllele = $TempArray[$CurrentSiteInTotalBiallelics];
							}
						
							$TempCounter++;
						
						}
					}
			
					close TOTALBIALLELICS;
				
				
					if ($CurrentDerivedAllele)  {
					
					
						my $PopulationStarter = 0;
					
						foreach my $population (1..$NumPops) {
						
							my $AncestralAlleleCount = 0;
							my $DerivedAlleleCount = 0;
							my $SampleSizeInCurrentPop = $DiploidSizes[$population-1];
						
							foreach my $allele (@CurrentSiteAlleles[$PopulationStarter..($PopulationStarter+$SampleSizeInCurrentPop-1)])  {
								if ($allele eq $CurrentAncestralAllele)  {
									$AncestralAlleleCount++;
								}
						
								elsif ($allele eq $CurrentDerivedAllele)  {
									$DerivedAlleleCount++;
								}
							}
				
							push (@AncestralAllelePopCounts, $AncestralAlleleCount);
			
							push (@DerivedAllelePopCounts, $DerivedAlleleCount);
							
							$PopulationStarter = $PopulationStarter+$SampleSizeInCurrentPop;
						}	
					
				
				
				
				
					
					}
		
				}			
					
		
		
		
		
		
		
				my $LengthDerivedAllelePopCounts = @DerivedAllelePopCounts;
		
				foreach my $count (@DerivedAllelePopCounts)  {
					$DerivedCountsHash{$count}++;
				}
		
		
		
				open OUTFILE, ">ResampledSFS/$RepDataset.resamp_DAFpop0.obs" or die$!;		#Create file for current population to store DFS.
		
				print OUTFILE "1 observations\n";     #\t$NumMonomorphicLoci monomorphic loci scored in all\t$NumPolymorphicLoci polymorphic loci scored in all\t$ReadLength bases in each read\n";
		
				foreach my $Pop1Count (0..$DiploidSizes[0]-1)  {
					print OUTFILE "d0\_$Pop1Count\t";
			
				}
				
				print OUTFILE "d0\_$DiploidSizes[0]\n";
		
				print OUTFILE "$TotalMonomorphicSites\t";
		
				foreach my $Pop1Count (1..$DiploidSizes[0]-1)  {
					print OUTFILE "$DerivedCountsHash{$Pop1Count}\t";
				}
		
				print OUTFILE "$DerivedCountsHash{$DiploidSizes[0]}";
		
		
				system "rm TempFiles/*";
	
		}
	
	
		system "rmdir TempFiles";
		
		print "Resampled unfolded site frequency spectra have been generated and printed to folder ResampledSFS in Formatting directory.\n";
		
	
	}	
	
	
	
	
	
	
		
	
	
	
	else {		#no resampling
	
	
	
		print "Working on creating site frequency spectrum\n\n";
	
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
	
		my $NumPolymorphicLoci = keys %AllPolymorphicLoci;
		my $TotalNumSNPs = @LocusNames;
	
		close FILE;
		close OUTFILE;
	
	
	
	
	
	
	
		#Replace any NA's in file SNAPP_Infile_NA.txt with "N".
	
		my $LineCounter = 0;
	
		open SNPFILEWITHNA, "TempFiles/SNPMatrix_Edit.txt" or die$!;
		open SNPFILENONA, ">TempFiles/SingleSNPsAllRaw.txt" or die$!;
	
			while (<SNPFILEWITHNA>)  {
				chomp($_);	
	
				if ($LineCounter == 0)  {
					$_ =~ s/"//g;
					$LineCounter++;
					print SNPFILENONA "$_";
				}
			
				else {
					$_ =~ s/"//g;
					my @TempArray = split(/\t/, $_);
					my $Length = @TempArray;
					print SNPFILENONA "$TempArray[0]\t";
				
					foreach my $allele (@TempArray[1..$Length-1])  {
						$allele =~ s/NA/N/g;
						print SNPFILENONA "$allele\t";
					}
				}
			
				print SNPFILENONA "\n";
			}
	
			close SNPFILEWITHNA;
			close SNPFILENONA;
	
	
	
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
		
	
	
	
		#Run R script "OutputBiallelicSingleSNPs".  This outputs file "UnlinkedBiallelicSNPs_Raw.txt". A maximum of one SNP is output for each locus.
			
		system "R --vanilla --slave < ../RScripts/OutputBiallelicSingleSNPs.R";
		
		
		my $NumUnlinkedBiallelicLoci;
		
		
		
		#Each line in UnlinkedBiallelicSNPsRaw.txt ends with \t\n and has quotes.  Remove these.
		
		if ($Unlinked == 1)  {
			open SNPFILE, "TempFiles/UnlinkedBiallelicSNPsRaw.txt" or die$!;
		}
		
		else {
			open SNPFILE, "TempFiles/AllBiallelicSNPsRaw.txt" or die$!;
		}
		
		
		open SNPFILEUPDATE, ">TempFiles/BiallelicSNPs_SFS.txt" or die$!;
		
		
		while (<SNPFILE>)  {
			
			if ($_ =~ /[A-Za-z1-9]/)  {
				$_ =~ s/\t\n$/\n/;
				$_ =~ s/"//g;
				my @TempArray = split(/\t/,$_);
				$NumUnlinkedBiallelicLoci = @TempArray;	#has an empty tab at beginning
				$NumUnlinkedBiallelicLoci = $NumUnlinkedBiallelicLoci-1;
				print SNPFILEUPDATE "$_";
			}	
		}
		 
		close SNPFILE;
		close SNPFILEUPDATE;
		
		
		
		##############################################################################################################################################
		
		
		print "\n\nNumber of biallelic loci for site frequency spectrum is $NumUnlinkedBiallelicLoci\n\n";
		
		my $TotalMonomorphicSites;
		my $RetainedReadLength;
		
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
				
			
		if ($ScaleMonos == 1)  {
			my $ProportionUnlinkedSNPs = $NumUnlinkedBiallelicLoci/$TotalNumSNPs;
			$TotalMonomorphicSites = int($TotalMonomorphicSites*$ProportionUnlinkedSNPs);
					
		}
		
		
		##############################################################################################################################################
		
		#Create hash with possible values of the derived allele counts as keys (keys 0 to Number of alleles sampled).
		
		
		my %DerivedCountsHash = ();
		
		foreach my $value (0..$DiploidSizes[0])  {
			
			$DerivedCountsHash{$value} = 0;
		}
		
		##############################################################################################################################################
		#Go through file UnlinkedBiallelicSNPs.txt and count ancestral and derived alleles at each locus.
		
		
		my @AncestralAllelePopCounts = ();
		my @DerivedAllelePopCounts = ();
		
		
		my $CurrentAncestralAllele;
		my $CurrentDerivedAllele;
					
		open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
		
		my @UnlinkedBiallelicsNames = ();
		
		while(<TOTALBIALLELICS>) {
			@UnlinkedBiallelicsNames = split(/\t/, $_);
			shift(@UnlinkedBiallelicsNames);
			last;
		}
		
		close TOTALBIALLELICS;
		
		
		
		foreach my $locusname (@UnlinkedBiallelicsNames) {
		
			my $CurrentSiteInTotalBiallelics = 0;
			my @CurrentSiteIndividualIDs = ();
			my @CurrentSiteAlleles = ();
			
			open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
				
			my $TempCounter = 0;
				
				while (<TOTALBIALLELICS>)  {	#get the column number for that locus in GenotypesUpdate.txt
					
					
				if ($TempCounter == 0)  {		#on the first line that has locus names.
					my @TempArray = split(/\t/,$_);
					my $TempArrayElementCounter = 0;
					
					foreach my $name (@TempArray)  {
						if ($locusname eq $name) {	#find the column number of the locus we are currently on
							$CurrentSiteInTotalBiallelics = $TempArrayElementCounter;
							last;
						}	
							
						else {
							$TempArrayElementCounter++;
						}
						
					}
						$TempCounter++;
						
				}
					
					elsif ($TempCounter == 1)  {	#Should be on the outgroup individual - this individual has to be placed in this first position in SNPMatrix_X.Y.txt
						
						my @TempArray = split(/\t/,$_);
						$CurrentAncestralAllele = $TempArray[$CurrentSiteInTotalBiallelics];
						
						$TempCounter++;
						
					}
					
					
					elsif ($TempCounter == 2)  {	#Second line of outgroup - ignore
						$TempCounter++;
						
					}
					
					
					
					else {
						
						my @TempArray = split(/\t/,$_);
						my $CurrentRead = $TempArray[$CurrentSiteInTotalBiallelics];
						push(@CurrentSiteIndividualIDs, $TempCounter);			#keep track of what sample we're working on
						push(@CurrentSiteAlleles, $TempArray[$CurrentSiteInTotalBiallelics]);
						
						if ($TempArray[$CurrentSiteInTotalBiallelics] !~ /$CurrentAncestralAllele/)	{
							$CurrentDerivedAllele = $TempArray[$CurrentSiteInTotalBiallelics];
						}
						
						$TempCounter++;
						
					}
				}
				
				close TOTALBIALLELICS;
				
				
				if ($CurrentDerivedAllele)  {
					
					
					my $PopulationStarter = 0;
					
					foreach my $population (1..$NumPops) {
					
						my $AncestralAlleleCount = 0;
						my $DerivedAlleleCount = 0;
						my $SampleSizeInCurrentPop = $DiploidSizes[$population-1];
						
						foreach my $allele (@CurrentSiteAlleles[$PopulationStarter..($PopulationStarter+$SampleSizeInCurrentPop-1)])  {
							if ($allele eq $CurrentAncestralAllele)  {
								$AncestralAlleleCount++;
							}
						
							elsif ($allele eq $CurrentDerivedAllele)  {
								$DerivedAlleleCount++;
							}
						}
				
						push (@AncestralAllelePopCounts, $AncestralAlleleCount);
			
						push (@DerivedAllelePopCounts, $DerivedAlleleCount);
						
						$PopulationStarter = $PopulationStarter+$SampleSizeInCurrentPop;
					}	
					
				
				
				
				
					
				}
		
		}			
					
		
		
		
		
		
		
		my $LengthDerivedAllelePopCounts = @DerivedAllelePopCounts;
		
		foreach my $count (@DerivedAllelePopCounts)  {
			$DerivedCountsHash{$count}++;
		}
		
		
		
		open OUTFILE, ">_DAFpop0.obs" or die$!;		#Create file for current population to store DFS.
		
		print OUTFILE "1 observations\n";     #\t$NumMonomorphicLoci monomorphic loci scored in all\t$NumPolymorphicLoci polymorphic loci scored in all\t$ReadLength bases in each read\n";
		
		foreach my $Pop1Count (0..$DiploidSizes[0]-1)  {
			print OUTFILE "d0\_$Pop1Count\t";
			
		}
		
		print OUTFILE "d0\_$DiploidSizes[0]\n";
		
		print OUTFILE "$TotalMonomorphicSites\t";
		
		foreach my $Pop1Count (1..$DiploidSizes[0]-1)  {
			print OUTFILE "$DerivedCountsHash{$Pop1Count}\t";
		}
		
		print OUTFILE "$DerivedCountsHash{$DiploidSizes[0]}";
		
		
		system "rm TempFiles/*";
		system "rmdir TempFiles";
	
		print "Unfolded site frequency spectrum has been generated and printed to Formatting directory with name _DAFpop0.obs.\n";
	}
}	
