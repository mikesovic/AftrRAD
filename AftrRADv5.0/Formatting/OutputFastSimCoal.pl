#! /usr/bin/perl
use warnings;
use strict;

#Before running this script, edit the file "Genotypes/Output/SNPMatrix_X.Y.txt" by reordering individuals into populations and keep track of the number of individuals per population.
#Also, an outgroup individual must be the first individual in this file.

#Built this script to produce multidimensional SFS from the March Bug Fix script.

print "\nInformation for OutputFastSimCoal.pl\n";
print "\nBefore running this script, make sure...\n";
print "1) samples in the file \"Genotypes/Output/SNPMatrix_X.Y.txt\" are ordered by population\n";
print "2) if creating unfolded SFS, the first sample in the SNPMatrix file should be the outgroup individual\n";
print "3) you know the number of individuals in each population, in the order they occur in the SNPMatrix file\n";
print "4) you have data for more than one population.  If you have data from a single population, run the script OutputFastSimCoal_SinglePop.pl instead.\n\n";


#Set parameters for the run

# -resamp Flag indicating whether to create resampled datasets.
# -num  Number of resampled datasets to create (default 1).
# -prop Percent of loci sampled to create the resampled datasets (1-99).
# -rep  Flag indicating whether to sample with replacement when generating resampled datasets.  Default is 0 (sample loci without replacement).
# -unlinked Flag indicating whether to include only unlinked SNPs in the site frequency spectrum.  Default is 0 (all SNPs are included).
# -multi  Flag indicating that a single multidimensional SFS will be produced.  Default is 0 (pairwise joint spectra will be output).
# -subsamp  Individuals to include in the SFS.  Indicate with a comma-delimited vector such as 3,4,9:11. This would include the third, fourth, ninth, tenth, and eleventh samples in the SNPMatrix file.
# -folded Indicates whether to output a folded or unfolded SFS.  Default is an unfolded SFS, which requires that an outgroup sample is the first sample in the SNPMatrix file.

my %RunArguments = ();

#Defaults
$RunArguments{resamp} = 0;
$RunArguments{num} = 1;
$RunArguments{pct} = 50;
$RunArguments{rep} = 0;
$RunArguments{unlinked} = 0;
$RunArguments{multi} = 0;
$RunArguments{folded} = 0;
$RunArguments{MonoScaled} = 0;

#read in command line arguments

for my $argument (@ARGV)  {
	my @TempArgs = split(/-/,$argument);
	
	#print help information
	if (($TempArgs[1] =~ /^h$/)|| ($TempArgs[1] =~ /^help$/)) {
		
		print "\nInformation for OutputFastSimCoal_JointSFS.pl...\n\n";
		print "This script creates site frequency spectra that can be used as input for analyses in FastSimCoal\n";  
		print "It requires a SNPMatrix file with no missing data (SNPMatrix_100.X.txt) in the Output/Genotypes folder. The samples in this file must be grouped by population.\n\n";
		
		print "\n\nCommand line arguments available for this script...\n\n";
		print "unlinked\nFlag indicating whether to include only unlinked SNPs in the site frequency spectrum.\n";
		print "Default is 0 (all SNPs are included)\n\n";
		print "resamp\nFlag indicating whether to create resampled datasets (resample loci).\n";
		print "Set this to 1 to perform resampling.\n";
		print "Resampled SNP datasets and associated SFS will be printed to folders \"ResampledDatasets\" and \"ResampledSFS\", respectively.\n\n";
		print "num\nNumber of resampled datasets to create (Default is 1).\n\n";
		print "pct\nPercent of loci sampled to create each resampled dataset (1-99, default is 50).\n\n";
		print "rep\nFlag indicating whether to sample with replacement when generating resampled datasets.  Default is 0 (sample loci without replacement).\n\n";
		print "multi\nFlag indicating that a single multidimensional SFS will be produced.  Default is 0 (pairwise joint spectra will be output).\n";
		print "Note that this is necessary if performing model choice (i.e. AIC).\n\n";
		print "folded\nFlag indicating whether to output a folded or unfolded SFS.  Default is an unfolded SFS,\n";
		print "which requires that an outgroup sample is the first sample in the SNPMatrix file.\n Set to '1' for folded SFS\n\n";
		print "MonoScaled\nFlag indicating whether to scale the number of monomorphic sites based on the proportion of SNPs retained after removing linked SNPs.\n";
		print "Default is 0, meaning all monomorphic sites are counted.  If set to '1', the proportion of total SNPs that are unlinked is calculated, and the total\n";
		print "number of monomorphic sites is scaled by this proportion\n\n.";
		
		exit;
		
	}
	
	#update default values with entered parameters as appropriate and continue with run
	elsif (exists($RunArguments{$TempArgs[0]}))  {
		$RunArguments{$TempArgs[0]} = $TempArgs[1];
	}
}


#Set command line arguments for the run

my $resamp = $RunArguments{resamp};
my $NumResampReps = $RunArguments{num};
my $PctLoci = $RunArguments{pct};
my $Replacement = $RunArguments{rep};
my $Unlinked = $RunArguments{unlinked};
my $multidim = $RunArguments{multi};
#my $SubsampleIndividuals = $RunArguments{subsamp};
my $Folded = $RunArguments{folded};
my $ScaleMonos = $RunArguments{MonoScaled};

my $FileName;

#read the SNPMatrix files available for the run

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
		print "\nEnter the name of the SNPMatrix file you want to use to create the FastSimCoal infile.\n";
	
		$FileName = <STDIN>;
		chomp($FileName);
		
}




print "\nEnter the number of individuals in each population, in the order the populations occur in the SNPMatrix_X.Y.txt file.  Separate each with a tab.\n";
print "If you're creating an unfolded SFS, do not include the outgroup sample here.\n";
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
	print "OutputFastSimCoal_SingleSFS.pl should be used to create SFS from a single population.";
	exit;
}

else {
	print "Recognized $NumPops populations.\n";  
}	
	

print "\n\n";







##############################################################################################################################################
##############################################################################################################################################
#If the resamp argument is set to 1, need to generate the resampled datasets.

if ($multidim == 1) {

	
	if ($Folded == 1) {
	
	
		if ($resamp == 1) {	#Multidim, Folded, with resamp
	
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
							if ($_ =~ /Sequence/) {
								next;
							}
							
							else {
								if ($_ =~ /[ATGC]/) {
									chomp($_);
									my @Temp = split(/\t/, $_);
									my $TempSeq = $Temp[0];
									$RetainedReadLength = length($TempSeq);
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
					my $ScaledMonomorphicLoci = int(($PctLoci/100)*$NumMonomorphicLoci);   #If resampled dataset includes less than 100% of the sites (based on maxSNP argument), need to scale the monomorphic sites.
					my $MonomorphicSites = 	$ScaledMonomorphicLoci*$RetainedReadLength;
					
					my $TotalLengthOfPolymorphicLoci = $NumPolymorphicLoci*$RetainedReadLength;
					my $NumMonomorphicSitesInPolymorphicLoci = $TotalLengthOfPolymorphicLoci-$TotalNumSNPs;
				
					$TotalMonomorphicSites = $MonomorphicSites+$NumMonomorphicSitesInPolymorphicLoci;
				
			
					if ($ScaleMonos == 1)  {
						my $ProportionUnlinkedSNPs = $NumUnlinkedBiallelicLoci/$TotalNumSNPs;
						$TotalMonomorphicSites = int($TotalMonomorphicSites*$ProportionUnlinkedSNPs);
					
					}
			
			
				
			##############################################################################################################################################
					
					#Define a hash that contains keys that have all combinations of counts of derived alleles in the form xxxyyyzzz (note that no population can have more than 999 samples).
					my %DerivedCountsHash = ();  
				
					
					#create all of the combinations
				
					my @StartArray = ();
					my @EndArray = ();
				
					my $FirstSampSize = $DiploidSizes[0];
					
					foreach my $value (0..$FirstSampSize) {
						if ($value < 10)  {
							$value = '00'.$value;
						}
						
						elsif ($value < 100) {
							$value = '0'.$value;
						}	
						
						push (@StartArray, $value);
					}	
					
					foreach my $sampsize (@DiploidSizes[1..$NumPops-1])  {
						
						my $CurrentMax = $sampsize;
						
						foreach my $value (@StartArray)  {
							
							foreach my $value2 (0..$CurrentMax) {
								if ($value2 < 10)  {
									$value2 = '00'.$value2;
								}
								
								elsif ($value2 < 100) {
									$value2 = '0'.$value2;
								}
								
								my $TempCombination = $value.$value2;
								
								push (@EndArray, $TempCombination);
							}
						}
						
						@StartArray = @EndArray;
						
						@EndArray = ();
					}
				
					
					#populate the hash with keys (combinations)
					
					foreach my $combination (@StartArray) {
						$DerivedCountsHash{$combination} = 0;
					}
				
				
					my $NumCombinations = keys(%DerivedCountsHash);
					
					print "\nNumber of combinations is $NumCombinations\n\n";
					
					
			
			
					
					
					
					#Go through file "TempFiles/BiallelicSNPs_SFS.txt" to get locus names
			
					my @UnlinkedBiallelicsNames = ();
	
					open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
				
					
					while(<TOTALBIALLELICS>) {
						@UnlinkedBiallelicsNames = split(/\t/, $_);
						shift(@UnlinkedBiallelicsNames);
						last;
					}
	
					close TOTALBIALLELICS;
			
			
			
		
					
					my @AncestralAllelePopCounts = ();
					my @DerivedAllelePopCounts = ();
					my @EqualFreqsFlags = ();		#This array keeps a flag (0 or 1) for each locus.  1's indicate that the locus had equal frequencies of the 
										#two alleles, and so ancestral/derived couldn't be determined.  
					
							
							
					#Define ancestral and derived alleles at each locus, as the spectrum is folded
							
							
							
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
							
							
						#At this point, we've counted the number of first alleles and second alleles at the current locus.
						#These are $NumFirstAlleles and $NumSecondAlleles
						
						my $CurrentDerivedAllele;
						my $CurrentAncestralAllele;
						my $EqualFreqFlag = 0;
								
						if ($NumFirstAlleles > $NumSecondAlleles)  {
							$CurrentDerivedAllele = $CurrentSecondAllele;
							$CurrentAncestralAllele = $CurrentFirstAllele;
								
						}
								
						elsif ($NumSecondAlleles > $NumFirstAlleles)  {
							$CurrentDerivedAllele = $CurrentFirstAllele;
							$CurrentAncestralAllele = $CurrentSecondAllele;
						}
								
						elsif ($NumSecondAlleles == $NumFirstAlleles) {
							$EqualFreqFlag = 1;
							$CurrentDerivedAllele = $CurrentSecondAllele;
							$CurrentAncestralAllele = $CurrentFirstAllele;
							
						}
								
								
						close TOTALBIALLELICS;
								
								
						#Reopen the TempFiles/BiallelicSNPs_SFS.txt file and find the current locus.
						#Then go through one population at a time, counting the number of derived alleles.
								
								
						open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
							
						$TempCounter = 0;
						my $CurrentPopEndLineNumber = $DiploidSizes[0];
						my $CurrentDerivedCountsValue = 0;
						my $CurrentAncestralCountsValue = 0;
						my $ConcatenatedDerivedCountsVector;
						my $ConcatenatedAncestralCountsVector;
						my $PopulationNumber = 0;
								
						while (<TOTALBIALLELICS>)  {
								
							if ($TempCounter == 0)  {
								$TempCounter++;
								next;	
							}
									
							my @TempArray = split(/\t/,$_);
							my $CurrentRead = $TempArray[$CurrentSiteInTotalBiallelics];
									
							if ($CurrentRead eq $CurrentDerivedAllele)  {
								$CurrentDerivedCountsValue++;
							}	
							
							else {
								$CurrentAncestralCountsValue++;	#This is only relevant if the EqualFreqsFlag is 1.
							}	
								
								
							if ($TempCounter == $CurrentPopEndLineNumber) {		#on the last line of the current pop
										
								if ($CurrentDerivedCountsValue<10) {
									$CurrentDerivedCountsValue = "00".$CurrentDerivedCountsValue;	
								}	
										
								elsif ($CurrentDerivedCountsValue<100) {
									$CurrentDerivedCountsValue = "0".$CurrentDerivedCountsValue;	
								}
								
								
								#This if statement only matters if the EqualFreqsFlag is 1.
								if ($CurrentAncestralCountsValue<10) {
									$CurrentAncestralCountsValue = "00".$CurrentAncestralCountsValue;	
								}	
										
								elsif ($CurrentAncestralCountsValue<100) {
									$CurrentAncestralCountsValue = "0".$CurrentAncestralCountsValue;	
								}
								
								
										
								if ($ConcatenatedDerivedCountsVector) {
									$ConcatenatedDerivedCountsVector = $ConcatenatedDerivedCountsVector.$CurrentDerivedCountsValue;
								}
										
								else {
									$ConcatenatedDerivedCountsVector = $CurrentDerivedCountsValue;
								}
								
								
								#This if statement only matters if the EqualFreqsFlag is 1.
								if ($ConcatenatedAncestralCountsVector) {
									$ConcatenatedAncestralCountsVector = $ConcatenatedAncestralCountsVector.$CurrentAncestralCountsValue;
								}
										
								else {
									$ConcatenatedAncestralCountsVector = $CurrentAncestralCountsValue;
								}
								
								
								
								
								$CurrentDerivedCountsValue = 0;
								$CurrentAncestralCountsValue = 0;
								$PopulationNumber++;
										
								if ($PopulationNumber == $NumPops) {
									last;
								}
										
								$CurrentPopEndLineNumber = $CurrentPopEndLineNumber+$DiploidSizes[$PopulationNumber];
										
										
							}
									
							$TempCounter++;
						}
										
						
						
						close TOTALBIALLELICS;
						
						#Add the current locus to the derived counts hash.
								
						if ($EqualFreqFlag == 0) {
							$DerivedCountsHash{$ConcatenatedDerivedCountsVector}++;
						}		
								
						else {	#EqualFreqFlag is 1, so can't determine ancestral/derived allele
							my $CurrentDerivedValue = $DerivedCountsHash{$ConcatenatedDerivedCountsVector};
							my $UpdatedDerivedValue = $CurrentDerivedValue+0.5;
							$DerivedCountsHash{$ConcatenatedDerivedCountsVector} = $UpdatedDerivedValue;
							
							my $CurrentAncestralValue = $DerivedCountsHash{$ConcatenatedDerivedCountsVector};
							my $UpdatedAncestralValue = $CurrentAncestralValue+0.5;
							$DerivedCountsHash{$ConcatenatedDerivedCountsVector} = $UpdatedAncestralValue;
						}		
								
					}
					
					
						
				
						mkdir "ResampledSFS/Resample$RepDataset";
						
						
						my @SortedHashValues = sort { $a<=>$b } keys %DerivedCountsHash;
			
						my $MonomorphicKey = 0 x (3*$NumPops);
						
						my $CurrentMonoValue = $DerivedCountsHash{$MonomorphicKey};
						
						$DerivedCountsHash{$MonomorphicKey} = $CurrentMonoValue+$TotalMonomorphicSites;
						
						
						
						open OUTFILE, ">ResampledSFS/Resample$RepDataset/_DSFS.obs" or die$!;
						
						print OUTFILE "1 observations\n";
						
						print OUTFILE "$NumPops\t";
						
						foreach my $samplesize (@DiploidSizes[0..$NumPops-2]) {
							print OUTFILE "$samplesize\t";
						}
						
						print OUTFILE "$DiploidSizes[$NumPops-1]\n";
						
						
						foreach my $combination (@SortedHashValues[0..$NumCombinations-2])  {
							print OUTFILE "$DerivedCountsHash{$combination}\t";
						}
						
						print OUTFILE "$DerivedCountsHash{$SortedHashValues[$NumCombinations-1]}";
					
						
						
					
				}	
					
				system "rm TempFiles/*";
				system "rmdir TempFiles";
			
		print "Resampled folded site frequency spectra have been generated and printed to folder ResampledSFS in Formatting directory\n";	
			
		}
		

		
	
		
		
		
		
		
		else {	#Multidim, Folded without resample	
			
			
			print "\n\nWorking on creating multidimensional site frequency spectrum\n\n";
				
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
				
				
				
				
			#Run R script "OutputBiallelicSingleSNPs".  This outputs two files: AllBiallelicSNPsRaw.txt and UnlinkedBiallelicSNPs_Raw.txt. A maximum of one SNP is output for each locus in UnlinkedBiallelicSNPs_Raw.txt.
					
			system "R --vanilla --slave < ../RScripts/OutputBiallelicSingleSNPs.R";
				
				
			my $NumUnlinkedBiallelicLoci;
				
				
				
			#Each line in the R output files ends with \t\n and has quotes.  Remove these in the appropriate file.
				
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
				
				
			print "\nNumber of biallelic loci for site frequency spectrum is $NumUnlinkedBiallelicLoci\n";
				
			
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
					if ($_ =~ /Sequence/) {
						next;
					}
					
					else {
						if ($_ =~ /[ATGC]/) {
							chomp($_);
							my @Temp = split(/\t/, $_);
							my $TempSeq = $Temp[0];
							$RetainedReadLength = length($TempSeq);
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
				
			#Define a hash that contains keys that have all combinations of counts of derived alleles in the form xxxyyyzzz (note that no population can have more than 999 samples).
			my %DerivedCountsHash = ();  
				
					
			#create all of the combinations
				
			my @StartArray = ();
			my @EndArray = ();
				
			my $FirstSampSize = $DiploidSizes[0];
					
			foreach my $value (0..$FirstSampSize) {
				if ($value < 10)  {
					$value = '00'.$value;
				}
						
				elsif ($value < 100) {
					$value = '0'.$value;
				}	
						
				push (@StartArray, $value);
			}	
					
			foreach my $sampsize (@DiploidSizes[1..$NumPops-1])  {
					
				my $CurrentMax = $sampsize;
						
				foreach my $value (@StartArray)  {
							
					foreach my $value2 (0..$CurrentMax) {
						if ($value2 < 10)  {
							$value2 = '00'.$value2;
						}
								
						elsif ($value2 < 100) {
							$value2 = '0'.$value2;
						}
								
						my $TempCombination = $value.$value2;
								
						push (@EndArray, $TempCombination);
					}
				}
						
				@StartArray = @EndArray;
						
				@EndArray = ();
			}
				
					
			#populate the hash with keys (combinations)
					
			foreach my $combination (@StartArray) {
				$DerivedCountsHash{$combination} = 0;
			}
				
				
			my $NumCombinations = keys(%DerivedCountsHash);
				
			print "\nNumber of combinations is $NumCombinations\n\n";
					
					
			
			
					
					
					
			#Go through file "TempFiles/BiallelicSNPs_SFS.txt" to get locus names
			
			open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die#!;
			
			my @UnlinkedBiallelicsNames = ();
	
			while(<TOTALBIALLELICS>) {
				@UnlinkedBiallelicsNames = split(/\t/, $_);
				shift(@UnlinkedBiallelicsNames);
				last;
			}
	
			close TOTALBIALLELICS;
			
			
			
			my @AncestralAllelePopCounts = ();
			my @DerivedAllelePopCounts = ();
			my @EqualFreqsFlags = ();		#This array keeps a flag (0 or 1) for each locus.  1's indicate that the locus had equal frequencies of the 
								#two alleles, and so ancestral/derived couldn't be determined.  
			
					
					
			#Define ancestral and derived alleles at each locus, as the spectrum is folded
					
					
					
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
					
					
				#At this point, we've counted the number of first alleles and second alleles at the current locus.
				#These are $NumFirstAlleles and $NumSecondAlleles
				
				my $CurrentDerivedAllele;
				my $CurrentAncestralAllele;
				my $EqualFreqFlag = 0;
						
				if ($NumFirstAlleles > $NumSecondAlleles)  {
					$CurrentDerivedAllele = $CurrentSecondAllele;
					$CurrentAncestralAllele = $CurrentFirstAllele;
						
				}
						
				elsif ($NumSecondAlleles > $NumFirstAlleles)  {
					$CurrentDerivedAllele = $CurrentFirstAllele;
					$CurrentAncestralAllele = $CurrentSecondAllele;
					
				}
						
				elsif ($NumSecondAlleles == $NumFirstAlleles) {
					$EqualFreqFlag = 1;
					$CurrentDerivedAllele = $CurrentSecondAllele;
					$CurrentAncestralAllele = $CurrentFirstAllele;
					
				}
						
						
				close TOTALBIALLELICS;
						
						
				#Reopen the TempFiles/BiallelicSNPs_SFS.txt file and find the current locus.
				#Then go through one population at a time, counting the number of derived alleles.
						
						
				open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
					
				$TempCounter = 0;
				my $CurrentPopEndLineNumber = $DiploidSizes[0];
				my $CurrentDerivedCountsValue = 0;
				my $CurrentAncestralCountsValue = 0;
				my $ConcatenatedDerivedCountsVector;
				my $ConcatenatedAncestralCountsVector;
				my $PopulationNumber = 0;
						
				while (<TOTALBIALLELICS>)  {
						
					if ($TempCounter == 0)  {
						$TempCounter++;
						next;	
					}
							
					my @TempArray = split(/\t/,$_);
					my $CurrentRead = $TempArray[$CurrentSiteInTotalBiallelics];
							
					if ($CurrentRead eq $CurrentDerivedAllele)  {
						$CurrentDerivedCountsValue++;
					}	
					
					else {
						$CurrentAncestralCountsValue++;	#This is only relevant if the EqualFreqsFlag is 1.
					}	
						
						
					if ($TempCounter == $CurrentPopEndLineNumber) {		#on the last line of the current pop
								
						if ($CurrentDerivedCountsValue<10) {
							$CurrentDerivedCountsValue = "00".$CurrentDerivedCountsValue;	
						}	
								
						elsif ($CurrentDerivedCountsValue<100) {
							$CurrentDerivedCountsValue = "0".$CurrentDerivedCountsValue;	
						}
						
						
						#This if statement only matters if the EqualFreqsFlag is 1.
						if ($CurrentAncestralCountsValue<10) {
							$CurrentAncestralCountsValue = "00".$CurrentAncestralCountsValue;	
						}	
								
						elsif ($CurrentAncestralCountsValue<100) {
							$CurrentAncestralCountsValue = "0".$CurrentAncestralCountsValue;	
						}
						
						
								
						if ($ConcatenatedDerivedCountsVector) {
							$ConcatenatedDerivedCountsVector = $ConcatenatedDerivedCountsVector.$CurrentDerivedCountsValue;
						}
								
						else {
							$ConcatenatedDerivedCountsVector = $CurrentDerivedCountsValue;
						}
						
						
						#This if statement only matters if the EqualFreqsFlag is 1.
						if ($ConcatenatedAncestralCountsVector) {
							$ConcatenatedAncestralCountsVector = $ConcatenatedAncestralCountsVector.$CurrentAncestralCountsValue;
						}
								
						else {
							$ConcatenatedAncestralCountsVector = $CurrentAncestralCountsValue;
						}
						
						
						
						
						$CurrentDerivedCountsValue = 0;
						$CurrentAncestralCountsValue = 0;
						$PopulationNumber++;
								
						if ($PopulationNumber == $NumPops) {
							last;
						}
								
						$CurrentPopEndLineNumber = $CurrentPopEndLineNumber+$DiploidSizes[$PopulationNumber];
								
								
					}
							
					$TempCounter++;
				}
								
						
				
				close TOTALBIALLELICS;
				
				#Add the current locus to the derived counts hash.
						
				if ($EqualFreqFlag == 0) {
					$DerivedCountsHash{$ConcatenatedDerivedCountsVector}++;
				}		
						
				else {	#EqualFreqFlag is 1, so can't determine ancestral/derived allele
					my $CurrentDerivedValue = $DerivedCountsHash{$ConcatenatedDerivedCountsVector};
					my $UpdatedDerivedValue = $CurrentDerivedValue+0.5;
					$DerivedCountsHash{$ConcatenatedDerivedCountsVector} = $UpdatedDerivedValue;
					
					my $CurrentAncestralValue = $DerivedCountsHash{$ConcatenatedDerivedCountsVector};
					my $UpdatedAncestralValue = $CurrentAncestralValue+0.5;
					$DerivedCountsHash{$ConcatenatedDerivedCountsVector} = $UpdatedAncestralValue;
				}		
						
			}
				
				
			my @SortedHashValues = sort { $a<=>$b } keys %DerivedCountsHash;
			
			my $MonomorphicKey = 0 x (3*$NumPops);
						
			my $CurrentMonoValue = $DerivedCountsHash{$MonomorphicKey};
						
			$DerivedCountsHash{$MonomorphicKey} = $CurrentMonoValue+$TotalMonomorphicSites;
						
						
			open OUTFILE, ">_DSFS.obs" or die$!;
			print OUTFILE "1 observations\n";
						
			print OUTFILE "$NumPops\t";
						
			foreach my $samplesize (@DiploidSizes[0..$NumPops-2]) {
				print OUTFILE "$samplesize\t";
			}
						
			print OUTFILE "$DiploidSizes[$NumPops-1]\n";
						
			
			##########################################################
			##########################################################
			
			foreach my $combination (@SortedHashValues[0..$NumCombinations-2])  {
				print OUTFILE "$DerivedCountsHash{$combination}\t";
			}
						
			print OUTFILE "$DerivedCountsHash{$SortedHashValues[$NumCombinations-1]}";
		
			
			print "Folded multidimensional site frequency spectrum has been generated and printed to Formatting directory with name _DSFS.obs\n";	
		
			
		}	
							
							
			
			
	
			
	}
	
	
	
	
	
	else {	#unfolded
	
	
	
		if ($resamp == 1)  {	#Multidim, Unfolded with resample
			
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
			
				#Get number of monomorphic loci
			
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
						if ($_ =~ /Sequence/) {
							next;
						}
						
						else {
							if ($_ =~ /[ATGC]/) {
								chomp($_);
								my @Temp = split(/\t/, $_);
								my $TempSeq = $Temp[0];
								$RetainedReadLength = length($TempSeq);
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
			
			
					
				my %DerivedCountsHash = ();  #Define a hash that contains keys that have all combinations of counts of derived alleles in the form xxyyzz (note that no population can have more than 99 samples.
				
				#create all of the combinations and add them to hash
				
				my @StartArray = ();
				my @EndArray = ();
				
				my $FirstSampSize = $DiploidSizes[0];
					
				foreach my $value (0..$FirstSampSize) {
					if ($value < 10)  {
						$value = '00'.$value;
					}
						
					elsif ($value < 100) {
						$value = '0'.$value;
					}	
						
					push (@StartArray, $value);
				}	
					
				foreach my $sampsize (@DiploidSizes[1..$NumPops-1])  {
						
					my $CurrentMax = $sampsize;
						
					foreach my $value (@StartArray)  {
							
						foreach my $value2 (0..$CurrentMax) {
							if ($value2 < 10)  {
								$value2 = '00'.$value2;
							}
								
							elsif ($value2 < 100) {
								$value2 = '0'.$value2;
							}
								
							my $TempCombination = $value.$value2;
								
							push (@EndArray, $TempCombination);
						}
					}
						
					@StartArray = @EndArray;
						
					@EndArray = ();
				}
				
					
				#populate the hash with keys
					
				foreach my $combination (@StartArray) {
					$DerivedCountsHash{$combination} = 0;
				}
				
					#for keys(%DerivedCountsHash) {
					#	print $_\t;
					#}
				
				my $NumCombinations = keys(%DerivedCountsHash);
					
				print "\nNumber of combinations is $NumCombinations\n\n";
					
					
					
					
					
					
					
			
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
							
							
						my $TempDerivedCombinationForHash;
								
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
									
							if ($DerivedAlleleCount < 10)  {
								$DerivedAlleleCount = '00'.$DerivedAlleleCount;
							}	
									
							elsif ($DerivedAlleleCount < 100) {
								$DerivedAlleleCount = '0'.$DerivedAlleleCount;
							}
									
									
							if ($TempDerivedCombinationForHash) {
								$TempDerivedCombinationForHash = $TempDerivedCombinationForHash.$DerivedAlleleCount;
							}	
									
							else {
								$TempDerivedCombinationForHash = $DerivedAlleleCount;
							}
									
									
							$PopulationStarter = $PopulationStarter+$SampleSizeInCurrentPop;
						}	
							
						$DerivedCountsHash{$TempDerivedCombinationForHash}++;
						
						$TempDerivedCombinationForHash = ();
						
						
							
					}
				
				}			
							
						
						
						
						
						
				
				mkdir "ResampledSFS/Resample$RepDataset";
						
						
				my @SortedHashValues = sort { $a<=>$b } keys %DerivedCountsHash;
			
				my $MonomorphicKey = 0 x (3*$NumPops);
						
				my $CurrentMonoValue = $DerivedCountsHash{$MonomorphicKey};
						
				$DerivedCountsHash{$MonomorphicKey} = $CurrentMonoValue+$TotalMonomorphicSites;
						
						
						
				open OUTFILE, ">ResampledSFS/Resample$RepDataset/_DSFS.obs" or die$!;
						
				print OUTFILE "1 observations\n";
						
				print OUTFILE "$NumPops\t";
						
				foreach my $samplesize (@DiploidSizes[0..$NumPops-2]) {
					print OUTFILE "$samplesize\t";
				}
						
				print OUTFILE "$DiploidSizes[$NumPops-1]\n";
						
						
				foreach my $combination (@SortedHashValues[0..$NumCombinations-2])  {
					print OUTFILE "$DerivedCountsHash{$combination}\t";
				}
						
				print OUTFILE "$DerivedCountsHash{$SortedHashValues[$NumCombinations-1]}";
						
			}	
							
							
				system "rm TempFiles/*";
				system "rmdir TempFiles";
						
				print "Resampled unfolded site frequency spectra have been generated and printed to folder ResampledSFS in Formatting directory\n";	
					
				
			
		}
	
	
	
			
			
		
		
		else {  #Multidim, Unfolded, no resample
			
			
			
			print "\n\nWorking on creating multidimensional site frequency spectrum\n\n";
				
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
				
				
				
				
			#Run R script "OutputBiallelicSingleSNPs".  This outputs two files: AllBiallelicSNPsRaw.txt and UnlinkedBiallelicSNPs_Raw.txt. A maximum of one SNP is output for each locus in UnlinkedBiallelicSNPs_Raw.txt.
					
			system "R --vanilla --slave < ../RScripts/OutputBiallelicSingleSNPs.R";
				
				
			my $NumUnlinkedBiallelicLoci;
				
				
				
			#Each line in the R output files ends with \t\n and has quotes.  Remove these in the appropriate file.
				
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
				
				
			print "\nNumber of biallelic loci for site frequency spectrum is $NumUnlinkedBiallelicLoci\n";
			
			##############################################################################################################################################
			#Go through file BiallelicSNPs_SFS.txt and count ancestral and derived alleles at each locus.
				
				
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
				
				
			my %DerivedCountsHash = ();  #Define a hash that contains keys that have all combinations of counts of derived alleles in the form xxyyzz (note that no population can have more than 99 samples.
				
			#create all of the combinations and add them to hash
				
			my @StartArray = ();
			my @EndArray = ();
			
			my $FirstSampSize = $DiploidSizes[0];
				
			foreach my $value (0..$FirstSampSize) {
				if ($value < 10)  {
					$value = '00'.$value;
				}
					
				elsif ($value < 100) {
					$value = '0'.$value;
				}	
					
				push (@StartArray, $value);
			}	
				
			foreach my $sampsize (@DiploidSizes[1..$NumPops-1])  {
					
				my $CurrentMax = $sampsize;
					
				foreach my $value (@StartArray)  {
						
					foreach my $value2 (0..$CurrentMax) {
						if ($value2 < 10)  {
							$value2 = '00'.$value2;
						}
							
						elsif ($value2 < 100) {
							$value2 = '0'.$value2;
						}
							
						my $TempCombination = $value.$value2;
							
						push (@EndArray, $TempCombination);
					}
				}
					
				@StartArray = @EndArray;
					
				@EndArray = ();
			}
			
				
			#populate the hash with keys
				
			foreach my $combination (@StartArray) {
				$DerivedCountsHash{$combination} = 0;
			}
			
				#for keys(%DerivedCountsHash) {
				#	print $_\t;
				#}
			
			my $NumCombinations = keys(%DerivedCountsHash);
				
			print "\nNumber of combinations is $NumCombinations\n\n";
						
						
				
				
				
				
				
				
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
								#This sample gives us the ancestral allele for the current locus.
								
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
						push(@CurrentSiteIndividualIDs, $TempCounter);			#keep track of what sample (individual) we're working on
						push(@CurrentSiteAlleles, $TempArray[$CurrentSiteInTotalBiallelics]);	#push the current allele to @CurrentSiteAlleles - this is paired with the array CurrentSiteIndividualIDs.
								
						if ($TempArray[$CurrentSiteInTotalBiallelics] !~ /$CurrentAncestralAllele/)	{	#Define the derived allele at the current locus.
							$CurrentDerivedAllele = $TempArray[$CurrentSiteInTotalBiallelics];
						}
								
						$TempCounter++;
								
					}
				}
						
				close TOTALBIALLELICS;
						
						
				if ($CurrentDerivedAllele)  {	#this should not be relevant anymore, as all loci should be biallelic and therefore have a derived allele, so this will always be true.
							
					my $TempDerivedCombinationForHash;
						
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
								
						if ($DerivedAlleleCount < 10)  {
							$DerivedAlleleCount = '00'.$DerivedAlleleCount;
						}	
								
						elsif ($DerivedAlleleCount < 100) {
							$DerivedAlleleCount = '0'.$DerivedAlleleCount;
						}
								
								
						if ($TempDerivedCombinationForHash) {
							$TempDerivedCombinationForHash = $TempDerivedCombinationForHash.$DerivedAlleleCount;
						}	
								
						else {
							$TempDerivedCombinationForHash = $DerivedAlleleCount;
						}	
								
						$PopulationStarter = $PopulationStarter+$SampleSizeInCurrentPop;
					}	
							
						
					$DerivedCountsHash{$TempDerivedCombinationForHash}++;
						
					$TempDerivedCombinationForHash = ();
						
							
				}
				
			}			
							
				
				
			
			
			
			
			
			my @SortedHashValues = sort { $a<=>$b } keys %DerivedCountsHash;
				
			#Add monomorphic sites to appropriate hash key.
				
				
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
					if ($_ =~ /Sequence/) {
						next;
					}
					
					else {
						if ($_ =~ /[ATGC]/) {
							chomp($_);
							my @Temp = split(/\t/, $_);
							my $TempSeq = $Temp[0];
							$RetainedReadLength = length($TempSeq);
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
				
			close MONOMORPHICINALL;
				
			my $NumMonomorphicLoci = @Monomorphics;
				
			my $MonomorphicSites = 	$NumMonomorphicLoci*$RetainedReadLength;
				
			my $TotalLengthOfPolymorphicLoci = $NumPolymorphicLoci*$RetainedReadLength;
			my $NumMonomorphicSitesInPolymorphicLoci = $TotalLengthOfPolymorphicLoci-$TotalNumSNPs;
					
			$TotalMonomorphicSites = $MonomorphicSites+$NumMonomorphicSitesInPolymorphicLoci;
				
			if ($ScaleMonos == 1)  {
				my $ProportionUnlinkedSNPs = $NumUnlinkedBiallelicLoci/$TotalNumSNPs;
				$TotalMonomorphicSites = int($TotalMonomorphicSites*$ProportionUnlinkedSNPs);
					
			}
			
			
			my $MonomorphicKey = 0 x (3*$NumPops);
				
			my $CurrentMonoValue = $DerivedCountsHash{$MonomorphicKey};
				
			$DerivedCountsHash{$MonomorphicKey} = $CurrentMonoValue+$TotalMonomorphicSites;
				
				
				
			open OUTFILE, ">_DSFS.obs" or die$!;
				
			print OUTFILE "1 observations\n";
				
			print OUTFILE "$NumPops\t";
				
			foreach my $samplesize (@DiploidSizes[0..$NumPops-2]) {
				print OUTFILE "$samplesize\t";
			}
				
			print OUTFILE "$DiploidSizes[$NumPops-1]\n";
				
				
			foreach my $combination (@SortedHashValues[0..$NumCombinations-2])  {
				print OUTFILE "$DerivedCountsHash{$combination}\t";
			}
				
			print OUTFILE "$DerivedCountsHash{$SortedHashValues[$NumCombinations-1]}";
				
				
					
					
			system "rm TempFiles/*";
			system "rmdir TempFiles";
			
			print "Multidimensional unfolded site frequency spectrum has been generated and printed to Formatting directory with name _DSFS.obs\n";	
		
		}	


	}

}











else {	#outputing joint SFS's
	
	
	if ($Folded == 1)  {

		if ($resamp == 1)  {	#Joint, folded, with resample
			
			
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
				#First, get the percent of loci scored and retained read length
				
				my @FirstSplit = split(/_/, $FileName);
				my $SecondElement = $FirstSplit[1];
				my @SecondSplit = split(/\./, $SecondElement);
				my $PctLociScored = $SecondSplit[0];
				$RetainedReadLength = $SecondSplit[1];
					
				if ($RetainedReadLength =~ /All/) {
					open MONOMORPHICINALL, "../Output/Genotypes/Monomorphics_$PctLociScored.txt" or die$!;
					while(<MONOMORPHICINALL>)  {
						if ($_ =~ /Sequence/) {
							next;
						}
						
						else {
							if ($_ =~ /[ATGC]/) {
								chomp($_);
								my @Temp = split(/\t/, $_);
								my $TempSeq = $Temp[0];
								$RetainedReadLength = length($TempSeq);
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
				
				
				
				#Go through file "TempFiles/BiallelicSNPs_SFS.txt" to get locus names
			
				
				my @UnlinkedBiallelicsNames = ();
		
				open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
				
				while(<TOTALBIALLELICS>) {
					@UnlinkedBiallelicsNames = split(/\t/, $_);
					shift(@UnlinkedBiallelicsNames);
					last;
				}
		
				close TOTALBIALLELICS;
				
						
				
				
				my @AncestralAllelePopCounts = ();
				my @DerivedAllelePopCounts = ();
				my @EqualFreqsFlags = ();		#This array keeps a flag (0 or 1) for each locus.  1's indicate that the locus had equal frequencies of the 
									#two alleles, and so ancestral/derived couldn't be determined.  
			
				
						
				#Define ancestral and derived alleles at each locus, as the spectrum is folded
					
						
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
						
						
					#At this point, we've counted the number of first alleles and second alleles at the current locus.
					#These are $NumFirstAlleles and $NumSecondAlleles
						
					my $CurrentDerivedAllele;
					my $CurrentAncestralAllele;
					my $EqualFreqFlag = 0;
							
					if ($NumFirstAlleles > $NumSecondAlleles)  {
						$CurrentDerivedAllele = $CurrentSecondAllele;
						$CurrentAncestralAllele = $CurrentFirstAllele;
							
					}
							
					elsif ($NumSecondAlleles > $NumFirstAlleles)  {
						$CurrentDerivedAllele = $CurrentFirstAllele;
						$CurrentAncestralAllele = $CurrentSecondAllele;
						
					}
							
					elsif ($NumSecondAlleles == $NumFirstAlleles) {
						$EqualFreqFlag = 1;
						$CurrentDerivedAllele = $CurrentSecondAllele;
						$CurrentAncestralAllele = $CurrentFirstAllele;
						
					}
							
							
					close TOTALBIALLELICS;
							
							
					#Reopen the TempFiles/BiallelicSNPs_SFS.txt file and find the current locus.
					#Then go through one population at a time, counting the number of derived alleles.
							
							
					open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
						
					$TempCounter = 0;
					my $CurrentPopEndLineNumber = $DiploidSizes[0];
					my $CurrentDerivedCounts = 0;
					my $CurrentAncestralCounts = 0;
			
					my $PopulationNumber = 0;
							
					while (<TOTALBIALLELICS>)  {
							
						if ($TempCounter == 0)  {
							$TempCounter++;
							next;	
						}
			
						
						my @TempArray = split(/\t/,$_);
						my $CurrentRead = $TempArray[$CurrentSiteInTotalBiallelics];
								
					
						
						if ($CurrentRead eq $CurrentDerivedAllele)  {
							$CurrentDerivedCounts++;
						}	
						
						
						else {
							$CurrentAncestralCounts++;
						}	
							
						if ($TempCounter == $CurrentPopEndLineNumber) {		#on the last line of the current pop
							
							push(@AncestralAllelePopCounts, $CurrentAncestralCounts);
							push(@DerivedAllelePopCounts, $CurrentDerivedCounts);
							push(@EqualFreqsFlags, $EqualFreqFlag);
							
							$CurrentDerivedCounts = 0;
							$CurrentAncestralCounts = 0;
							$PopulationNumber++;
									
							if ($PopulationNumber == $NumPops) {
								last;
							}
									
							$CurrentPopEndLineNumber = $CurrentPopEndLineNumber+$DiploidSizes[$PopulationNumber];
									
									
						}
								
						$TempCounter++;
					}
			
					close TOTALBIALLELICS;
			
				}
	
			
				
				
				mkdir "ResampledSFS/Resample$RepDataset";
			
				foreach my $pop (0..$NumPops-2) {			#Do pairwaise comparisons for each pair of populations.
						
					my $MonomorphicSites2 = 0;
						
					foreach my $comparison ($pop+1..$NumPops-1)  {
							
							
						my $FirstPop = '_'.$pop;
						my $SecondPop = $comparison;
						my $CurrentComparison = $SecondPop . $FirstPop;		
							
						open OUTFILE, ">ResampledSFS/Resample$RepDataset/_jointMAFpop$CurrentComparison.obs" or die$!;		#Create file for current pair of populations to store DFS.
							
						
						
						my $NumSamplesInFirstPop = $DiploidSizes[$pop];			#Get number of alleles sampled in each population for the current pair of pops.
						my $NumSamplesInSecondPop = $DiploidSizes[$comparison];
							
						print OUTFILE "1 observations\n";
						print OUTFILE "\t";
							
						foreach my $Pop1Count (0..$NumSamplesInFirstPop-1)  {		#Print the header line in the outfile
							print OUTFILE "d$pop\_$Pop1Count\t";
								
						}	
							
						print OUTFILE "d$pop\_$NumSamplesInFirstPop\n";
							
							
							
						#Get number of monomorphic sites in current pair of populations
							
						my $NumSitesInPolymorphicReads = $NumPolymorphicLoci*$RetainedReadLength;
					
						my $NumMonomorphicSitesInPolymorphicLoci = $NumSitesInPolymorphicReads-$TotalNumSNPs;
					
						$MonomorphicSites2 = $MonomorphicSites+$NumMonomorphicSitesInPolymorphicLoci;
												
							
						my %DerivedCountsHash = ();
								
						foreach my $samples1 (0..$NumSamplesInFirstPop)  {		#Populate the hash with keys that account for all of the possibilities, and assign each a value of zero.
								
								
							foreach my $samples2 (0..$NumSamplesInSecondPop)  {
									
									
								my $combin = $samples1.'a'.$samples2;
									
								$DerivedCountsHash{$combin} = 0;
							}
								
								
								
								
						}	
								
								
								
								
									
						#Go through array to get paired counts of derived allele for each locus for the current pop comparison.
								
						foreach my $currentlocusnumber (0..$NumUnlinkedBiallelicLoci-1) {	#For each locus
									
							my $FirstPosition = $pop+($NumPops*$currentlocusnumber);
							my $SecondPosition = $comparison+($NumPops*$currentlocusnumber);
	
							if ($EqualFreqsFlags[$FirstPosition] == 1)  {	#This means the derived allele couldn't be identified at the locus
								
								#Add 0.5 to appropriate hash value for "derived" allele
								my $FirstCount = $DerivedAllelePopCounts[$FirstPosition];
								my $SecondCount = $DerivedAllelePopCounts[$SecondPosition];
								my $CurrentHashLookup = $FirstCount.'a'.$SecondCount;
								my $CurrentHashValue = $DerivedCountsHash{$CurrentHashLookup};
								my $NewValue = $CurrentHashValue+0.5;
								$DerivedCountsHash{$CurrentHashLookup} = $NewValue;
								
								#Add 0.5 to appropriate hash value for "ancestral" allele
								
								$FirstCount = $AncestralAllelePopCounts[$FirstPosition];
								$SecondCount = $AncestralAllelePopCounts[$SecondPosition];
								$CurrentHashLookup = $FirstCount.'a'.$SecondCount;
								$CurrentHashValue = $DerivedCountsHash{$CurrentHashLookup};
								$NewValue = $CurrentHashValue+0.5;
								$DerivedCountsHash{$CurrentHashLookup} = $NewValue;
								
								
							}	
						
							
							else {	
								
								
								my $FirstCount = $DerivedAllelePopCounts[$FirstPosition];
								my $SecondCount = $DerivedAllelePopCounts[$SecondPosition];
									
								my $CurrentHashLookup = $FirstCount.'a'.$SecondCount;
									
								$DerivedCountsHash{$CurrentHashLookup}++;
							}
						}
								
									
							
						#Add monomorphic sites from monomorphic loci to monomorphic sites counted above - these monomrphic sites that were counted
						#above can only occur when more than two populations are in the dataset.  For example, if there are 3 populations, and
						#one of the three has a singleton, then when the other two are compared, they will not have any derived sites.
							
							
						my $CurrentMonomorphicSitesSum = $DerivedCountsHash{'0a0'};
						my $TotalMonomorphicSites2 = $CurrentMonomorphicSitesSum+$TotalMonomorphicSites;
					
						$DerivedCountsHash{'0a0'} = $TotalMonomorphicSites2;
							
							
							
						foreach my $Pop2Count (0..$NumSamplesInSecondPop)  {			#For each row in the DAFS
								
							print OUTFILE "d$comparison\_$Pop2Count\t";
								
							foreach my $pop1potentialcount (0..$NumSamplesInFirstPop-1)  {
									
								my $combin2 = $pop1potentialcount.'a'.$Pop2Count;
									
								print OUTFILE "$DerivedCountsHash{$combin2}\t";
									
							}
								
							my $combin3 = $NumSamplesInFirstPop.'a'.$Pop2Count;
								
							print OUTFILE "$DerivedCountsHash{$combin3}\n";	
							
						}	
							
									
								
						
					}
				
				
				}
				
				
				
			}
		
		
		print "Resampled folded joint site frequency spectra have been generated and printed to folder ResampledSFS in Formatting directory\n";	
		
		
		}
		
		
	
		
		
		
		
		else {	#Joint, folded, no resample
		
		
			print "\n\nWorking on creating joint site frequency spectra\n\n";
			
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
			
			
			
			
			#Run R script "OutputBiallelicSingleSNPs".  This outputs two files: All BiallelicSNPsRaw.txt and UnlinkedBiallelicSNPs_Raw.txt. A maximum of one SNP is output for each locus.
				
			system "R --vanilla --slave < ../RScripts/OutputBiallelicSingleSNPs.R";
			
			
			my $NumUnlinkedBiallelicLoci;
			
			
			
			#Each line in the R output files ends with \t\n and has quotes.  Remove these in the appropriate file.
			
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
			
			
			print "\nNumber of biallelic loci for site frequency spectrum is $NumUnlinkedBiallelicLoci\n";
			
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
					if ($_ =~ /Sequence/) {
						next;
					}
					
					else {
						if ($_ =~ /[ATGC]/) {
							chomp($_);
							my @Temp = split(/\t/, $_);
							my $TempSeq = $Temp[0];
							$RetainedReadLength = length($TempSeq);
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
			
	
			
			
			#Go through file "TempFiles/BiallelicSNPs_SFS.txt" to get locus names
			
			my @UnlinkedBiallelicsNames = ();
	
			open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
			
			while(<TOTALBIALLELICS>) {
				@UnlinkedBiallelicsNames = split(/\t/, $_);
				shift(@UnlinkedBiallelicsNames);
				last;
			}
	
			close TOTALBIALLELICS;
			
					
			
			
			my @AncestralAllelePopCounts = ();
			my @DerivedAllelePopCounts = ();
			my @EqualFreqsFlags = ();		#This array keeps a flag (0 or 1) for each locus.  1's indicate that the locus had equal frequencies of the 
								#two alleles, and so ancestral/derived couldn't be determined.  
			
			
					
			#Define ancestral and derived alleles at each locus, as the spectrum is folded
				
					
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
					
					
				#At this point, we've counted the number of first alleles and second alleles at the current locus.
				#These are $NumFirstAlleles and $NumSecondAlleles
				
				my $CurrentDerivedAllele;
				my $CurrentAncestralAllele;
				
				my $EqualFreqFlag = 0;
						
				if ($NumFirstAlleles > $NumSecondAlleles)  {
					$CurrentDerivedAllele = $CurrentSecondAllele;
					$CurrentAncestralAllele = $CurrentFirstAllele;
						
				}
						
				elsif ($NumSecondAlleles > $NumFirstAlleles)  {
					$CurrentDerivedAllele = $CurrentFirstAllele;
					$CurrentAncestralAllele = $CurrentSecondAllele;
					
				}
						
				elsif ($NumSecondAlleles == $NumFirstAlleles) {
					$EqualFreqFlag = 1;
					$CurrentDerivedAllele = $CurrentSecondAllele;
					$CurrentAncestralAllele = $CurrentFirstAllele;
					
				}
						
						
				close TOTALBIALLELICS;
						
						
				#Reopen the TempFiles/BiallelicSNPs_SFS.txt file and find the current locus.
				#Then go through one population at a time, counting the number of derived alleles.
						
						
				open TOTALBIALLELICS, "TempFiles/BiallelicSNPs_SFS.txt" or die$!;
					
				$TempCounter = 0;
				my $CurrentPopEndLineNumber = $DiploidSizes[0];
				my $CurrentDerivedCounts = 0;
				my $CurrentAncestralCounts = 0;
				my $ConcatenatedDerivedCountsVector;
				my $PopulationNumber = 0;
				
				while (<TOTALBIALLELICS>)  {
						
					if ($TempCounter == 0)  {
						$TempCounter++;
						next;	
					}
		
					
					
					my @TempArray = split(/\t/,$_);
					my $CurrentRead = $TempArray[$CurrentSiteInTotalBiallelics];
							
					
					if ($CurrentRead eq $CurrentDerivedAllele)  {
						$CurrentDerivedCounts++;
					}	
					
					
					else {
						$CurrentAncestralCounts++;
					}
					
						
					if ($TempCounter == $CurrentPopEndLineNumber) {		#on the last line of the current pop for the current locus
					
						
						push(@AncestralAllelePopCounts, $CurrentAncestralCounts);
						push(@DerivedAllelePopCounts, $CurrentDerivedCounts);
						push(@EqualFreqsFlags, $EqualFreqFlag);
						
						$CurrentDerivedCounts = 0;
						$CurrentAncestralCounts = 0;
						$PopulationNumber++;
								
						if ($PopulationNumber == $NumPops) {
							last;
						}
						
						$CurrentPopEndLineNumber = $CurrentPopEndLineNumber+$DiploidSizes[$PopulationNumber];
								
								
					}
							
					$TempCounter++;
				}
			
			
				close TOTALBIALLELICS;
				
			}
			
			
			foreach my $pop (0..$NumPops-2) {			#Do pairwaise comparisons for each pair of populations.
				
				my $MonomorphicSites2 = 0;
				
				foreach my $comparison ($pop+1..$NumPops-1)  {
					
					
					my $FirstPop = '_'.$pop;
					my $SecondPop = $comparison;
					my $CurrentComparison = $SecondPop . $FirstPop;		
					
					open OUTFILE, ">_jointMAFpop$CurrentComparison.obs" or die$!;		#Create file for current pair of populations to store DFS.
					
				
				
					my $NumSamplesInFirstPop = $DiploidSizes[$pop];			#Get number of alleles sampled in each population for the current pair of pops.
					my $NumSamplesInSecondPop = $DiploidSizes[$comparison];
					
					print OUTFILE "1 observations\n";
					print OUTFILE "\t";
					
					foreach my $Pop1Count (0..$NumSamplesInFirstPop-1)  {		#Print the header line in the outfile
						print OUTFILE "d$pop\_$Pop1Count\t";
						
					}	
					
					print OUTFILE "d$pop\_$NumSamplesInFirstPop\n";
					
					
					
					#Get number of monomorphic sites in current pair of populations
					
					my $NumSitesInPolymorphicReads = $NumPolymorphicLoci*$RetainedReadLength;
			
					my $NumMonomorphicSitesInPolymorphicLoci = $NumSitesInPolymorphicReads-$TotalNumSNPs;
					$MonomorphicSites2 = $MonomorphicSites+$NumMonomorphicSitesInPolymorphicLoci;
					
					
					
					
					
					
					
					
					my %DerivedCountsHash = ();
						
					foreach my $samples1 (0..$NumSamplesInFirstPop)  {		#Populate the hash with keys that account for all of the possibilities, and assign each a value of zero.
						
						
						foreach my $samples2 (0..$NumSamplesInSecondPop)  {
							
							
							my $combin = $samples1.'a'.$samples2;
							
							$DerivedCountsHash{$combin} = 0;
						}
						
						
						
						
					}	
						
						
						
						
							
					#Go through array to get paired counts of derived allele for each locus for the current pop comparison.
						
					foreach my $currentlocusnumber (0..$NumUnlinkedBiallelicLoci-1) {	#For each locus
						
						my $FirstPosition = $pop+($NumPops*$currentlocusnumber);
						my $SecondPosition = $comparison+($NumPops*$currentlocusnumber);

						if ($EqualFreqsFlags[$FirstPosition] == 1)  {	#This means the derived allele couldn't be identified at the locus
							
							#Add 0.5 to appropriate hash value for "derived" allele
							my $FirstCount = $DerivedAllelePopCounts[$FirstPosition];
							my $SecondCount = $DerivedAllelePopCounts[$SecondPosition];
							my $CurrentHashLookup = $FirstCount.'a'.$SecondCount;
							my $CurrentHashValue = $DerivedCountsHash{$CurrentHashLookup};
							my $NewValue = $CurrentHashValue+0.5;
							$DerivedCountsHash{$CurrentHashLookup} = $NewValue;
							
							
							#Add 0.5 to appropriate hash value for "ancestral" allele
							
							$FirstCount = $AncestralAllelePopCounts[$FirstPosition];
							$SecondCount = $AncestralAllelePopCounts[$SecondPosition];
							$CurrentHashLookup = $FirstCount.'a'.$SecondCount;
							$CurrentHashValue = $DerivedCountsHash{$CurrentHashLookup};
							$NewValue = $CurrentHashValue+0.5;
							$DerivedCountsHash{$CurrentHashLookup} = $NewValue;
							
						}	
					
						
						else {	
							
							
							my $FirstCount = $DerivedAllelePopCounts[$FirstPosition];
							my $SecondCount = $DerivedAllelePopCounts[$SecondPosition];
								
							my $CurrentHashLookup = $FirstCount.'a'.$SecondCount;
								
							$DerivedCountsHash{$CurrentHashLookup}++;
						}	
				
					}	
							
					
					#Add monomorphic sites from monomorphic loci to monomorphic sites counted above - these monomrphic sites that were counted
					#above can only occur when more than two populations are in the dataset.  For example, if there are 3 populations, and
					#one of the three has a singleton, then when the other two are compared, they will not have any derived sites.
					
					
					my $CurrentMonomorphicSitesSum = $DerivedCountsHash{'0a0'};
					my $TotalMonomorphicSites2 = $CurrentMonomorphicSitesSum+$TotalMonomorphicSites;
					
					$DerivedCountsHash{'0a0'} = $TotalMonomorphicSites2;
					
					
					
					foreach my $Pop2Count (0..$NumSamplesInSecondPop)  {			#For each row in the DAFS
						
						print OUTFILE "d$comparison\_$Pop2Count\t";
						
						foreach my $pop1potentialcount (0..$NumSamplesInFirstPop-1)  {
							
							my $combin2 = $pop1potentialcount.'a'.$Pop2Count;
							
							print OUTFILE "$DerivedCountsHash{$combin2}\t";
							
						}
						
						my $combin3 = $NumSamplesInFirstPop.'a'.$Pop2Count;
						
						print OUTFILE "$DerivedCountsHash{$combin3}\n";	
					
					}	
					
							
						
				
				}
				
				
				
			}
			
			
			
		
			
			
			
			
		print "Joint folded site frequency spectra have been generated and printed to Formatting directory with names _jointMAFpopX_Y.obs.\n";	
			
	
			
		}	
		
		
	}	
	
	
	
	
	
	
	
	
	else {
	
		if ($resamp == 1)  {	#Joint, unfolded, with resample
	
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
			
				#Create the current resampled dataset - each is stored in a folder named ResampledDatasets
		
				print "Creating resampled dataset $RepDataset of $NumResampReps\n";
				system "R --vanilla --slave < ../RScripts/Resample_Edit.R";
		
				#Rename and move the output from R
				my $TempFileName = $FileName;
				$TempFileName =~ s/.txt//;
				system "mv ResampledDatasets/SNPMatrix_resamp.txt ResampledDatasets/$TempFileName.$RepDataset.txt";
				
				
				
				#Create SFS from current resampled dataset
				open FILE, "ResampledDatasets/$TempFileName.$RepDataset.txt" or die$!;
				open OUTFILE, ">TempFiles/SNPMatrix_Edit.txt" or die$!;
		
				my @LocusNames = ();
				my $LineNumber = 0;
		
				
				#Clean up locus names from R - create file TempFiles/SNPMatrix_Edit.txt which will be used going forward
				
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
		
				#The hash gets rid of duplicate locus names, so this gives us the actual number of polymorphic loci
				foreach my $name (@LocusNames)  {
					$AllPolymorphicLoci{$name}=1;
				}
		
				my $NumPolymorphicLoci = keys %AllPolymorphicLoci;
				my $TotalNumSNPs = @LocusNames;
		
				close FILE;
				close OUTFILE;
		
		
		
		
		
		
		
				#Replace any NA's in file TempFiles/SNPMatrix_Edit.txt with "N".
		
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
		
		
		
		
				#Run R script "OutputBiallelicSingleSNPs".  This outputs two files: "UnlinkedBiallelicSNPs_Raw.txt" and "AllBiallelicSNPsRaw.txt". A maximum of one SNP is output for each locus in the unlinked file.
			
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
		
				#The variable "$NumUnlinkedBiallelicLoci" may not necessarily represent unlinked loci - if $Unlinked is set to zero, this variable will represent all loci.
				print "\n\nNumber of biallelic loci for site frequency spectrum is $NumUnlinkedBiallelicLoci\n\n";
		
				#Get number of monomorphic loci
		
				
				
				
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
						if ($_ =~ /Sequence/) {
							next;
						}
						
						else {
							if ($_ =~ /[ATGC]/) {
								chomp($_);
								my @Temp = split(/\t/, $_);
								my $TempSeq = $Temp[0];
								$RetainedReadLength = length($TempSeq);
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
				my $ScaledMonomorphicLoci = int(($PctLoci/100)*$NumMonomorphicLoci);	 #If resampled dataset includes less than 100% of the loci, need to scale the monomorphic sites.
				my $MonomorphicSites = 	$ScaledMonomorphicLoci*$RetainedReadLength;
						
					
				my $TotalLengthOfPolymorphicLoci = $NumPolymorphicLoci*$RetainedReadLength;
				my $NumMonomorphicSitesInPolymorphicLoci = $TotalLengthOfPolymorphicLoci-$TotalNumSNPs;
					
				$TotalMonomorphicSites = $MonomorphicSites+$NumMonomorphicSitesInPolymorphicLoci;
					
				
				if ($ScaleMonos == 1)  {
					my $ProportionUnlinkedSNPs = $NumUnlinkedBiallelicLoci/$TotalNumSNPs;
					$TotalMonomorphicSites = int($TotalMonomorphicSites*$ProportionUnlinkedSNPs);
					
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
						
			
					mkdir "ResampledSFS/Resample$RepDataset";
			
					foreach my $pop (0..$NumPops-2) {			#Do pairwaise comparisons for each pair of populations.
						
						my $MonomorphicSites2 = 0;
						
						foreach my $comparison ($pop+1..$NumPops-1)  {
							
							
							my $FirstPop = '_'.$pop;
							my $SecondPop = $comparison;
							my $CurrentComparison = $SecondPop . $FirstPop;		
							
							open OUTFILE, ">ResampledSFS/Resample$RepDataset/_jointDAFpop$CurrentComparison.obs" or die$!;		#Create file for current pair of populations to store DFS.
							
						
						
							my $NumSamplesInFirstPop = $DiploidSizes[$pop];			#Get number of alleles sampled in each population for the current pair of pops.
							my $NumSamplesInSecondPop = $DiploidSizes[$comparison];
							
							print OUTFILE "1 observations\n";
							print OUTFILE "\t";
							
							foreach my $Pop1Count (0..$NumSamplesInFirstPop-1)  {		#Print the header line in the outfile
								print OUTFILE "d$pop\_$Pop1Count\t";
								
							}	
							
							print OUTFILE "d$pop\_$NumSamplesInFirstPop\n";
							
							$MonomorphicSites2 = $TotalMonomorphicSites;
							
							
							
							
							
							
							my %DerivedCountsHash = ();
								
							foreach my $samples1 (0..$NumSamplesInFirstPop)  {		#Populate the hash with keys that account for all of the possibilities, and assign each a value of zero.
								
								
								foreach my $samples2 (0..$NumSamplesInSecondPop)  {
									
									
									my $combin = $samples1.'a'.$samples2;
									
									$DerivedCountsHash{$combin} = 0;
								}
								
								
								
								
							}	
								
								
								
								
									
								#Go through array to get paired counts of derived allele for each locus for the current pop comparison.
								
								foreach my $currentlocusnumber (0..$NumUnlinkedBiallelicLoci-1) {	#For each locus
									
									my $FirstPosition = $pop+($NumPops*$currentlocusnumber);
									my $SecondPosition = $comparison+($NumPops*$currentlocusnumber);
									
									my $FirstCount = $DerivedAllelePopCounts[$FirstPosition];
									my $SecondCount = $DerivedAllelePopCounts[$SecondPosition];
									
									my $CurrentHashLookup = $FirstCount.'a'.$SecondCount;
									
									$DerivedCountsHash{$CurrentHashLookup}++;
								}
								
									
							
							#Add monomorphic sites from monomorphic loci to monomorphic sites counted above - these monomrphic sites that were counted
							#above can only occur when more than two populations are in the dataset.  For example, if there are 3 populations, and
							#one of the three has a singleton, then when the other two are compared, they will not have any derived sites.
							
							
							my $CurrentMonomorphicSitesSum = $DerivedCountsHash{'0a0'};
							my $TotalMonomorphicSites2 = $CurrentMonomorphicSitesSum+$TotalMonomorphicSites;
					
							$DerivedCountsHash{'0a0'} = $TotalMonomorphicSites2;
							
							
							
							foreach my $Pop2Count (0..$NumSamplesInSecondPop)  {			#For each row in the DAFS
								
								print OUTFILE "d$comparison\_$Pop2Count\t";
								
								foreach my $pop1potentialcount (0..$NumSamplesInFirstPop-1)  {
									
									my $combin2 = $pop1potentialcount.'a'.$Pop2Count;
									
									print OUTFILE "$DerivedCountsHash{$combin2}\t";
									
								}
								
								my $combin3 = $NumSamplesInFirstPop.'a'.$Pop2Count;
								
								print OUTFILE "$DerivedCountsHash{$combin3}\n";	
							
							}	
							
									
								
						
						}
						
						
						
					}
			
			
			
			}
	
	
			system "rm TempFiles/*";
			system "rmdir TempFiles";
			
			print "Resampled unfolded joint site frequency spectra have been generated and printed to folder ResampledSFS in Formatting directory\n";	
			
		}
	
	
	
	
	
	
		else {  #Joint, unfolded, no resampling
		
		
		
			print "\n\nWorking on creating joint site frequency spectra\n\n";
			
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
			
			
			
			
			#Run R script "OutputBiallelicSingleSNPs".  This outputs two files: All BiallelicSNPsRaw.txt and UnlinkedBiallelicSNPs_Raw.txt. A maximum of one SNP is output for each locus.
				
			system "R --vanilla --slave < ../RScripts/OutputBiallelicSingleSNPs.R";
			
			
			my $NumUnlinkedBiallelicLoci;
			
			
			
			#Each line in the R output files ends with \t\n and has quotes.  Remove these in the appropriate file.
			
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
			
			
			print "\nNumber of biallelic loci for site frequency spectrum is $NumUnlinkedBiallelicLoci\n";
			
			
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
					if ($_ =~ /Sequence/) {
						next;
					}
					
					else {
						if ($_ =~ /[ATGC]/) {
							chomp($_);
							my @Temp = split(/\t/, $_);
							my $TempSeq = $Temp[0];
							$RetainedReadLength = length($TempSeq);
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
			
			#Go through file BiallelicSNPs_SFS.txt and count ancestral and derived alleles at each locus.
			
			
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
										#This sample gives us the ancestral allele for the current locus.
							
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
						
			
			
			
			
			
			
			foreach my $pop (0..$NumPops-2) {			#Do pairwaise comparisons for each pair of populations.
				
				my $MonomorphicSites2 = 0;
				
				foreach my $comparison ($pop+1..$NumPops-1)  {
					
					
					my $FirstPop = '_'.$pop;
					my $SecondPop = $comparison;
					my $CurrentComparison = $SecondPop . $FirstPop;		
					
					open OUTFILE, ">_jointDAFpop$CurrentComparison.obs" or die$!;		#Create file for current pair of populations to store DFS.
					
				
				
					my $NumSamplesInFirstPop = $DiploidSizes[$pop];			#Get number of alleles sampled in each population for the current pair of pops.
					my $NumSamplesInSecondPop = $DiploidSizes[$comparison];
					
					print OUTFILE "1 observations\n";
					print OUTFILE "\t";
					
					foreach my $Pop1Count (0..$NumSamplesInFirstPop-1)  {		#Print the header line in the outfile
						print OUTFILE "d$pop\_$Pop1Count\t";
						
					}	
					
					print OUTFILE "d$pop\_$NumSamplesInFirstPop\n";
					
					
					
					$MonomorphicSites2 = $TotalMonomorphicSites;
					
					
					
					
					
					
					my %DerivedCountsHash = ();
						
					foreach my $samples1 (0..$NumSamplesInFirstPop)  {		#Populate the hash with keys that account for all of the possibilities, and assign each a value of zero.
						
						
						foreach my $samples2 (0..$NumSamplesInSecondPop)  {
							
							
							my $combin = $samples1.'a'.$samples2;
							
							$DerivedCountsHash{$combin} = 0;
						}
						
						
						
						
					}	
						
						
						
						
							
						#Go through array to get paired counts of derived allele for each locus for the current pop comparison.
						
						foreach my $currentlocusnumber (0..$NumUnlinkedBiallelicLoci-1) {	#For each locus
							
							my $FirstPosition = $pop+($NumPops*$currentlocusnumber);
							my $SecondPosition = $comparison+($NumPops*$currentlocusnumber);
							
							my $FirstCount = $DerivedAllelePopCounts[$FirstPosition];
							my $SecondCount = $DerivedAllelePopCounts[$SecondPosition];
							
							my $CurrentHashLookup = $FirstCount.'a'.$SecondCount;
							
							$DerivedCountsHash{$CurrentHashLookup}++;
						}
						
							
					
					#Add monomorphic sites from monomorphic loci to monomorphic sites counted above - these monomrphic sites that were counted
					#above can only occur when more than two populations are in the dataset.  For example, if there are 3 populations, and
					#one of the three has a singleton, then when the other two are compared, they will not have any derived sites.
					
					
					my $CurrentMonomorphicSitesSum = $DerivedCountsHash{'0a0'};
					my $TotalMonomorphicSites2 = $CurrentMonomorphicSitesSum+$TotalMonomorphicSites;
					
					$DerivedCountsHash{'0a0'} = $TotalMonomorphicSites2;
					
					
					
					foreach my $Pop2Count (0..$NumSamplesInSecondPop)  {			#For each row in the DAFS
						
						print OUTFILE "d$comparison\_$Pop2Count\t";
						
						foreach my $pop1potentialcount (0..$NumSamplesInFirstPop-1)  {
							
							my $combin2 = $pop1potentialcount.'a'.$Pop2Count;
							
							print OUTFILE "$DerivedCountsHash{$combin2}\t";
							
						}
						
						my $combin3 = $NumSamplesInFirstPop.'a'.$Pop2Count;
						
						print OUTFILE "$DerivedCountsHash{$combin3}\n";	
					
					}	
					
							
						
				
				}
				
				
				
			}
		
		
			system "rm TempFiles/*";
			system "rmdir TempFiles";
			
			print "Joint unfolded site frequency spectra have been generated and printed to Formatting directory with names _jointDAFpopX_Y.obs.\n";	
		}	

	}
}	
