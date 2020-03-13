#! /usr/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;


#Set parameters for the run

# -MinReads  Minimum coverage required at a locus in an individual to apply a binomial test and call a genotype.
# -pvalLow   p-value applied for scoring genotypes at loci that have less than pvalThresh total reads.
# -pvalHigh  p-value applied for scoring genotypes at loci that have more than than pvalThresh total reads.
# -pvalThresh  Threshold number of reads at each locus that determines which p-value to use for the binomial test to call the genotype.
# -subset  Only a subset of the samples will be analyzed.  

my %RunArguments = ();

$RunArguments{MinReads} = 10;
$RunArguments{pvalLow} = 0.0001;
$RunArguments{pvalHigh} = 0.00001;
$RunArguments{pvalThresh} = 100;
$RunArguments{Help} = 0;
$RunArguments{subset} = 0;
$RunArguments{subsetfile} = 'SamplesForSubset.txt';
$RunArguments{maxProcesses} = 1;

for my $argument (@ARGV)  {
	my @TempArgs = split(/-/,$argument);
	
	#print help information
	if (($TempArgs[1] =~ /^h$/)|| ($TempArgs[1] =~ /^help$/)) {
		print "Information for Genotype.pl...\n\n";
		print "Genotype.pl should be run after the AftrRAD.pl script has completed.\n";
		print "Genotype.pl performs a number of operations that include calling genotypes at each locus, calculating the proportion of missing data in each sample, removing bad samples from further analyses, and identifying SNP locations along the reads.\n";
		print "\n\n Command line arguments available in Genotype.pl...\n\n";
		print "MinReads\nThe total number of reads at each locus required to call a genotype in an individual.\n";
		print "Default is 10.\n\n";
		print "pvalLow\np-value applied for scoring genotypes at loci that have less than pvalThresh total reads.\n";
		print "Default is 0.0001\n\n";
		print "pvalHigh\np-value applied for scoring genotypes at loci that have more than than pvalThresh total reads.\n";
		print "Default is 0.00001\n\n";
		print "pvalThresh\nThreshold number of reads at each locus that determines the p-value for the binomial test when calling genotypes.\n";
		print "Default is 100\n\n";
		print "subset\nOption to include only a subset of the individuals in genotyping.\n";
		print "Default (0) is to include all samples.  If set to '1', a text file must be provided with the names of the samples to include.\n";
		print "The default name for this file is 'SamplesForSubset.txt', and it should be located in the main AftrRAD directory (the same place as Genotypes.pl).\n";
		print "This file should contain one sample name per line.\n\n";
		print "subsetfile\nFile name/path containing sample names to include in genotyping, if different from the default above.\n";
		print "Only applies if 'subset' argument is set to '1'.\n\n";
		print "maxProcesses\nThe maximum number of processors to use for parallel runs.\n";
		print "Default is 1.  Currently the initial filtering and the ACANA alignment step are set up to run in parallel.\n\n";
		
		
		exit;
	}	
	
	#update default values with entered parameters as appropriate and continue with run
	elsif (exists($RunArguments{$TempArgs[0]}))  {
		$RunArguments{$TempArgs[0]} = $TempArgs[1];
	}
}

my $MinReads = $RunArguments{MinReads};
my $pvalLow = $RunArguments{pvalLow};
my $pvalHigh = $RunArguments{pvalHigh};
my $pvalThresh = $RunArguments{pvalThresh};
my $Subset = $RunArguments{subset};
my $SubFileName = $RunArguments{subsetfile};
my $MaxProcesses = $RunArguments{maxProcesses};


print "\n\n";


print "Arguments entered are...\n";
for (keys %RunArguments) {
	print "$_\t$RunArguments{$_}\n";
}	


print "\n\nRunning Genotype.pl...\n\n";

#############################################################################################################
#Get names of all samples in the analysis, and store them in array AllSampleNames
#Omit any samples that did not have any reads assigned to them.

my @OmittedSamples = ();  #These are samples that are in the MasterBarcodeFile, but did not have any unique reads assigned to them.
my @AllSampleNames = ();  #All samples with reads assigned to them.

my $TotalNumberOfIndividuals = 0;

if ($Subset == 1)  {
	
	open CURRENTBARCODEFILE, "$SubFileName" or die$!;

	while (<CURRENTBARCODEFILE>)  {
		chomp ($_);
		my @TabSeparatedArray = split (/\t/, "$_");
				
		if (-e "out/TempFiles/ForBinomialTest$TabSeparatedArray[0].txt")  {
			push (@AllSampleNames, $TabSeparatedArray[0]);
			$TotalNumberOfIndividuals++;
		}	
		
		else {
			push(@OmittedSamples, $TabSeparatedArray[1]);
		}	
	}
}



else {
	open CURRENTBARCODEFILE, "out/TempFiles/MasterBarcodeFile.txt" or die$!;
	
	while (<CURRENTBARCODEFILE>)  {
		chomp ($_);
		my @TabSeparatedArray = split (/\t/, "$_");
				
		if (-e "out/TempFiles/ForBinomialTest$TabSeparatedArray[1].txt")  {
			push (@AllSampleNames, $TabSeparatedArray[1]);
			$TotalNumberOfIndividuals++;
		}	
		
		else {
			push(@OmittedSamples, $TabSeparatedArray[1]);
		}	
	}
}



close CURRENTBARCODEFILE;




############################################################################################################
############################################################################################################
#For each individual, run R script on each pair of alleles at each locus and use binomial test to genotype.
#Whether a given locus is genotyped in a specific individual depends on whether the number of reads for that individual
#at that locus exceeds a certain threshold (MinReads above).  The default threshold is 10, but that can be changed with a command line argument.
#Print scored genotypes to file Genotypes.txt
############################################################################################################
############################################################################################################

#First, update the Genotype.R template script to set the desired coverage threshold and p-values.
 

open GENOTYPERSCRIPT, "RScripts/Genotype.R" or die$!;
open GENOTYPERSCRIPTEDIT, ">RScripts/Genotype_Edit.R" or die$!;
 
 
while(<GENOTYPERSCRIPT>)  {
	 if ($_ =~ /NumTrials<5/) {
 	 	 print GENOTYPERSCRIPTEDIT "if (NumTrials<$MinReads) {\n";
 	 	 next;
 	 }
 	 
 	 elsif ($_ =~ /NumTrials<100/)  {
 	 	print GENOTYPERSCRIPTEDIT "if (NumTrials<$pvalThresh) {\n";
 	 	next;
 	 }
 	 
 	 elsif ($_ =~ /homozygoteLow/)  {
 	 	 print GENOTYPERSCRIPTEDIT "if (CumProb<$pvalLow) { #Have a homozygoteLow\n";
 	 	 next;
 	 }
 	 
 	 elsif ($_ =~ /homozygoteHigh/)  {
 	 	 print GENOTYPERSCRIPTEDIT "if (CumProb<$pvalHigh) { #Have a homozygoteHigh\n";
 	 	 next;
 	 }
	 
 	 else {
 	 	 print GENOTYPERSCRIPTEDIT "$_";
 	 }
}	 
 


############################################################################################################
print "Scoring genotypes at all nonparalogous loci.\n\n";
print "Currently genotyping sample...";
 

#Perform binomial tests for genotyping.


if ($MaxProcesses > 1)  {
	
	#print "\nGenotyping sample...\n";
	
	#get number of subdirectories for parallel run
	
	my $NumSubdirectories = $MaxProcesses;
	my @SubprocessDirectories = ();
	
	#create subdirectories in TempFiles folder
	
	my $LengthAllSamples = @AllSampleNames;
	
	if ($MaxProcesses > $LengthAllSamples)  {
		
		foreach my $number (0..$LengthAllSamples-1)  {
			mkdir "out/TempFiles/SubGenotypes_$number";
			push (@SubprocessDirectories, "SubGenotypes_$number");
			open GENS, ">out/TempFiles/SubGenotypes_$number/SubGenotypes.txt" or die$!;
			open NAMES, ">out/TempFiles/SubGenotypes_$number/SubNames.txt" or die$!;
			close GENS;
			close NAMES;
		}
	}	
		
	else {
		foreach my $number (0..$NumSubdirectories-1)  {
			mkdir "out/TempFiles/SubGenotypes_$number";
			push (@SubprocessDirectories, "SubGenotypes_$number");
			open GENS, ">out/TempFiles/SubGenotypes_$number/SubGenotypes.txt" or die$!;
			open NAMES, ">out/TempFiles/SubGenotypes_$number/SubNames.txt" or die$!;
			close GENS;
			close NAMES;
		}
	}	
	
	
	#split up sample names into subdirs for parallel run
	
	open CURRENTBARCODEFILE, "out/TempFiles/MasterBarcodeFile.txt" or die$!;
	
	my $LineCounter = 0;
	
	while (<CURRENTBARCODEFILE>)  {
		chomp ($_);
		my @TabSeparatedArray = split (/\t/, "$_");
				
		if (-e "out/TempFiles/ForBinomialTest$TabSeparatedArray[1].txt")  {
			open NAMES, ">>out/TempFiles/SubGenotypes_$LineCounter/SubNames.txt" or die$!;
			print NAMES "$TabSeparatedArray[1]\n";
			close NAMES;
		}
		
		if ($LineCounter == $MaxProcesses-1)  {
			$LineCounter =0;
		}
		
		else {
			$LineCounter++;
		}
		
	
	}


	my $GenoManager = new Parallel::ForkManager( $MaxProcesses );	
		
	foreach my $subdir (@SubprocessDirectories)  {	
				
		$GenoManager->start and next;
	
		open NAMES, "out/TempFiles/$subdir/SubNames.txt" or die $!;
		
		open GENOTYPESFILE, ">>out/TempFiles/$subdir/SubGenotypes.txt" or die$!;
		
		open TEMPGEN, ">out/TempFiles/$subdir/TempGenotypes.txt" or die$!;
		
		while(<NAMES>) {
			
			
			if ($~ =~ /[A-Za-z0-9]/) {
				
				chomp($_);
				my $CurrentName = $_;
				
				print "$CurrentName\t";
				
				open TEMPGEN, ">out/TempFiles/$subdir/TempGenotypes.txt" or die$!;
				close TEMPGEN;
				
				open LOCIFILE, "out/TempFiles/ForBinomialTest$CurrentName.txt" or die$!;
				open TEMPLOCIFILE, ">out/TempFiles/$subdir/TempBiallelicLociForGenotyping.txt" or die$!;
		 
				while (<LOCIFILE>)  {
					print TEMPLOCIFILE "$_";
				}
		 
				close LOCIFILE;
				close TEMPLOCIFILE;
				
				system "cp RScripts/Genotype_Edit.R out/TempFiles/$subdir";
				
				
				
				open GENOTYPERSCRIPT, "out/TempFiles/$subdir/Genotype_Edit.R" or die$!;
				open GENOTYPERSCRIPTEDIT, ">out/TempFiles/$subdir/Genotype_SubEdit.R" or die$!;
				 
				 
				while(<GENOTYPERSCRIPT>)  {
					 if ($_ =~ /Genotypes<-read/) {
						 print GENOTYPERSCRIPTEDIT "Genotypes<-read.table(file=\"out/TempFiles/$subdir/TempBiallelicLociForGenotyping.txt\")\n";
						 next;
					 }
					 
					 elsif ($_ =~ /TempGenotypes/)  {
						my $SubstituteText = 'out/TempFiles/'.$subdir.'/TempGenotypes';
					 	$_ =~ s/TempGenotypes/$SubstituteText/;
					 	print GENOTYPERSCRIPTEDIT "$_";
						next;
					 }
					 
					 else {
					 	 print GENOTYPERSCRIPTEDIT "$_";
					 }	 
					 
					 
				}
				
				close GENOTYPERSCRIPT;
				close GENOTYPERSCRIPTEDIT;
				
				print "Running R\n";
				
				system "R --vanilla --slave < out/TempFiles/$subdir/Genotype_SubEdit.R";
				
				print "Finished R\n";
	
				open TEMPGENOTYPES, "out/TempFiles/$subdir/TempGenotypes.txt" or die$!;
				my $LineCounter = 0;
				
				
				
				
				print GENOTYPESFILE "Individual$CurrentName\t";
			 
				while (<TEMPGENOTYPES>)  {	#print the first allele at each locus for the current individual to file Genotypes
					$LineCounter++;
			   
					if ($LineCounter == 2)  {
						chomp ($_);
						print GENOTYPESFILE "$_\t";
					}
			   
					if ($LineCounter == 3)  {
						$LineCounter = 0;
					}
				}
				
				print GENOTYPESFILE "\n";
				close TEMPGENOTYPES;
			 
			 
				open TEMPGENOTYPES, "out/TempFiles/$subdir/TempGenotypes.txt" or die$!;
				$LineCounter = 0;
				print GENOTYPESFILE "Individual$CurrentName\t";
			       
				while (<TEMPGENOTYPES>)  {	#print the second allele at each locus for the current individual to file Genotypes
					$LineCounter++;     
			   
					if ($LineCounter == 3)  {
						chomp ($_);
						print GENOTYPESFILE "$_\t";
						$LineCounter = 0;
					}
				}
				
				print GENOTYPESFILE "\n";
		 
				system "rm out/TempFiles/$subdir/TempGenotypes.txt";
				
			}	
	
	
	
		}	
	
		$GenoManager->finish;
	}
	
	$GenoManager->wait_all_children;
	
	
	open GENOTYPESFILE, ">out/TempFiles/GenotypesLocusNamesForCat.txt" or die$!;
	
	my $LociPrinted = 0;
	
	foreach my $f (@AllSampleNames)  {
	     
		 if (-e "out/TempFiles/ForBinomialTest$f.txt")  {	
		 	 
		 	 open LOCIFILE, "out/TempFiles/ForBinomialTest$f.txt" or die$!;
		 	 
		 	 print GENOTYPESFILE "\t";
		 	 
		 	 while (<LOCIFILE>)  {
				if ($_ =~ /Locus/)  {
					my @TempArray = split(/\s/, $_);
					print GENOTYPESFILE "$TempArray[0]\t";
				}
			}
			print GENOTYPESFILE "\n";
			
			$LociPrinted = 1;
		}	
		 	 
		if ($LociPrinted == 1)  {
			last;
		}
	}	
		 	 
	close GENOTYPESFILE;	 	 
		
	
	my $FilesToCat;
	
	
	for my $filenumber (0..$NumSubdirectories-1) {		
		$FilesToCat = $FilesToCat."out/TempFiles/SubGenotypes_$filenumber/SubGenotypes.txt ";		
	}
		
	system "cat out/TempFiles/GenotypesLocusNamesForCat.txt $FilesToCat > out/TempFiles/Genotypes.txt";
	
}	









else {
	
	#print "\nGenotyping sample...\n";
	
	open GENOTYPESFILE, ">out/TempFiles/Genotypes.txt" or die$!;
	
	my $FirstIndividual = 1;
	
	foreach my $f (@AllSampleNames)  {
	     
		 if (-e "out/TempFiles/ForBinomialTest$f.txt")  {
		 
			 print "$f\t";
		 
			 open LOCIFILE, "out/TempFiles/ForBinomialTest$f.txt" or die$!;
			 open TEMPLOCIFILE, ">out/TempFiles/TempBiallelicLociForGenotyping.txt" or die$!;
		 
			 while (<LOCIFILE>)  {
				 print TEMPLOCIFILE "$_";
			 }
		 
			close LOCIFILE;
			close TEMPLOCIFILE;
		 
			system "R --vanilla --slave < RScripts/Genotype_Edit.R";  	
		 
		 
		 
			if ($FirstIndividual == 1)  {		#If first individual, print the locus names to file TempGenotypes.txt before moving to the genotypes for this individual.
				open TEMPGENOTYPES, "TempGenotypes.txt" or die$!;
				print GENOTYPESFILE "\t";
		   
				while (<TEMPGENOTYPES>)  {
					if ($_ =~ /Locus/)  {
						chomp ($_);
						print GENOTYPESFILE "$_\t";
					}
				}  
				close TEMPGENOTYPES;
				$FirstIndividual = 0;
				print GENOTYPESFILE "\n";
			} 
		  
		 
		 
			open TEMPGENOTYPES, "TempGenotypes.txt" or die$!;
			my $LineCounter = 0;
			print GENOTYPESFILE "Individual$f\t";
		 
			while (<TEMPGENOTYPES>)  {	#print the first allele at each locus for the current individual to file Genotypes
				$LineCounter++;
		   
				if ($LineCounter == 2)  {
					chomp ($_);
					print GENOTYPESFILE "$_\t";
				}
		   
				if ($LineCounter == 3)  {
					$LineCounter = 0;
				}
			}
			
			print GENOTYPESFILE "\n";
			close TEMPGENOTYPES;
		 
		 
			open TEMPGENOTYPES, "TempGenotypes.txt" or die$!;
			$LineCounter = 0;
			print GENOTYPESFILE "Individual$f\t";
		       
			while (<TEMPGENOTYPES>)  {	#print the second allele at each locus for the current individual to file Genotypes
				$LineCounter++;     
		   
				if ($LineCounter == 3)  {
					chomp ($_);
					print GENOTYPESFILE "$_\t";
					$LineCounter = 0;
				}
			}
			
			print GENOTYPESFILE "\n";
	 
			system "rm TempGenotypes.txt";
	 
		 } 
	 }

	 close GENOTYPESFILE;

}	 
 
 
 


###########################################################################################################
###########################################################################################################
##Evaluate for potentially bad individuals (# of NA's in SNPMatrixAll is > 2 stdev above mean).
###########################################################################################################
###########################################################################################################
 
#First, create SNPMatrixAll, which has every polymorphic locus not identified as paralogous - doesn't matter how many individuals it was scored in.

print "Creating SNPMatrixAll from R Script\n";

system "R --vanilla --slave < RScripts/OutputSNPMatrixAll.R";


###########################################################################################################
#Get the number of NA's for each individual and store it in array NACounts
print "\n\nChecking for individuals with large amounts of missing data (have >2 StDev more missing data than the average)\n";


open SNPMATRIX, "out/TempFiles/SNPMatrixAll.txt" or die$!;

my $Starter = 0;
my @SampleNames = ();
my @NACounts = ();
my $OddNumberLine = 1;
my %AllLocusNames = ();
my @AllLocusNamesArray = ();
my $NumOfLoci = 0;

while (<SNPMATRIX>)  {
	
	#Get count of total number of unique polymorphic loci in the dataset.
	
	
	if ($Starter == 0)  {	#on the first line (locus names)
		$Starter++;
		@AllLocusNamesArray = split(/\t/,$_);
		
		for my $name (@AllLocusNamesArray)  {
			if (exists($AllLocusNames{$name}))  {
				next;
			}
			
			else {
				$NumOfLoci++;
				$AllLocusNames{$name} = 1;
			}	
		
		}
	}
	
	
	

	else {			#On a line with genotypes for one of the samples.  Count number of loci with NA's.
	
	    my %LociWithNAs = ();
		
	    if ($OddNumberLine == 1)  {	
		
	    	my $TempNACount = 0;
		chomp($_);
		my @CurrentGenotypes = split (/\t/, $_);
		push (@SampleNames, $CurrentGenotypes[0]);
		
		my $CurrentLocationInArray = 1;
		my $CurrentGenotypesLength = @CurrentGenotypes;
		
		foreach my $genotype (@CurrentGenotypes[1..$CurrentGenotypesLength-1]) {
			
			if ($genotype =~ /NA/)  {
				
				my $CurrentLocus = $AllLocusNamesArray[$CurrentLocationInArray];
				
				if (exists($LociWithNAs{$CurrentLocus}))  {	#make sure that loci that have multiple SNPs only get counted once if they have NA.
					$CurrentLocationInArray++;
					next;
				}
				
				else {
					$LociWithNAs{$CurrentLocus} = 1;
					$TempNACount++;
					$CurrentLocationInArray++;
				}
			}
		}
		
		push (@NACounts, $TempNACount);
		$TempNACount = 0;
		
		$OddNumberLine = 0;
		
	    }
		
	    else {
			$OddNumberLine = 1;
	    }
	}
}




###########################################################################################################
#Get the proportion of missing data for each individual and print it to file MissingData.txt.

my @NAProportions;
my $SampleElementCounter = 0;

for my $count (@NACounts)  {
	
	my $TempProportion = $count/$NumOfLoci;
	push (@NAProportions, $TempProportion);
}

open NAPROPORTIONSFILE, ">out/OutputRunInfo/MissingDataProportions.txt" or die$!;

for my $name (@SampleNames)  {
	
	print NAPROPORTIONSFILE "$name\t";
	print NAPROPORTIONSFILE "$NAProportions[$SampleElementCounter]\n";
	$SampleElementCounter++;
	
}	


close NAPROPORTIONSFILE;

##############################################################################################################################################################################		
#Get the mean and stdev of the NACount values

my $NumOfIndividuals = @NACounts;

my $TotalNACount = 0;

foreach my $count (@NACounts)  {
	$TotalNACount = $TotalNACount + $count;
}

my $MeanCount = $TotalNACount/$NumOfIndividuals;

my @SquaredDifferences = ();

foreach my $count (@NACounts) {
	my $TempDifference = $count - $MeanCount;
	my $SquaredDifference = $TempDifference*$TempDifference;
	push (@SquaredDifferences, $SquaredDifference);
}

my $SumSquaredDifferences = 0;

foreach my $SquaredDifference (@SquaredDifferences)  {
	$SumSquaredDifferences = $SumSquaredDifferences + $SquaredDifference;
}

my $MeanSquaredDiff = $SumSquaredDifferences/$NumOfIndividuals;

my $StDev = sqrt($MeanSquaredDiff);

close SNPMATRIX;


#############################################################################################################################################################################		
#Identify individuals with >2 stdev more missing data than the mean.

print "Identifying samples with >2 stdev missing data\n";

open FLAGGEDSAMPLES, ">out/TempFiles/FlaggedSamples.txt" or die$!;
open GOODSAMPLES, ">out/TempFiles/GoodSamples.txt" or die$!;

my $BadSample = 0;
my $TwoStDev = 2*$StDev;
my $NAThresholdHigh = $MeanCount+$TwoStDev;
my $CurrentSample = 0;
my @BadSamples = ();
my $CurrentIndividualNumber = 0;

open NACOUNTSTOPLOT, ">out/TempFiles/NACountsToPlot.txt" or die$!;
open SAMPLENAMESTOPLOT, ">out/TempFiles/NASampleNamesToPlot.txt" or die$!;
open NACOUNTS, ">out/TempFiles/NACounts.txt" or die$!;

foreach my $name (@SampleNames)  {
	print SAMPLENAMESTOPLOT "$name\t";
	print NACOUNTS "$name\t";
	print NACOUNTS "$NACounts[$CurrentIndividualNumber]\n";
	$CurrentIndividualNumber++;
}	

foreach my $count (@NACounts)  {
	print NACOUNTSTOPLOT "$count\t";
	
	$BadSample = 0;
	
	if ($count > $NAThresholdHigh)  {	#Current sample has >2 stdev more missing data than the mean.
		$BadSample = 1;
	}

	if ($BadSample == 1)  {
		push (@BadSamples, $SampleNames[$CurrentSample]);
		my $TempSample = $SampleNames[$CurrentSample];
		$TempSample =~ s/"//g;
		print FLAGGEDSAMPLES "$TempSample\n";
	}

	$BadSample = 0;
	$CurrentSample++;
}	

close NACOUNTSTOPLOT;
close SAMPLENAMESTOPLOT;
close FLAGGEDSAMPLES;

system "R --vanilla --slave < RScripts/PlotNACounts.R";






############################################################################################################################
#Determine whether any samples will be left out of further analyses.  
#The script makes suggestions about potential bad samples, but the user makes the final decision.

open REMOVEDSAMPLES, ">out/TempFiles/MasterReport.txt" or die$!;

if ($OmittedSamples[0])  {
	print REMOVEDSAMPLES "Samples @OmittedSamples were omitted from the analysis because no unique reads were assigned to them.\n\n";
}

print REMOVEDSAMPLES "A total of $TotalNumberOfIndividuals individuals were analyzed.\n\n";



my %SamplesToRemove = ();

if (@BadSamples)  {
	print "\n\nSamples identified as potentially bad samples are @BadSamples.\n\n";
	
	print "Do you want to remove all(a), some(s), or none(n) of these potentially bad samples from the analysis? (a/s/n)\n";
	my $Remove = <STDIN>;
	chomp($Remove);
	
	my %BadNamesHash = ();
	my @BadSamplesUpdate = ();

	if ($Remove =~ /s/)  {
		
		 print "Enter the names of the samples you want to remove. Separate each with a tab\n";
		 my $RemoveString = <STDIN>;
		 chomp($RemoveString);
		
		 if ($RemoveString =~ /\t/)  {
		 
		 	 my @SamplesToRemove = split(/\t/, $RemoveString);
		  
		
		 	 foreach my $name (@SamplesToRemove)  {
		 	 	 $name =~ s/Individual//;
		 	 	 $SamplesToRemove{$name} = 1;
		 	 }
		 }
		 
		 else {
		 	 $RemoveString =~ s/Individual//;
		 	 $SamplesToRemove{$RemoveString} = 1;
		 }	 
		
	}
		 
		 
	elsif ($Remove =~ /a/)  {	 
		 
		foreach my $badsample (@BadSamples)  {
			
			$SamplesToRemove{$badsample} = 1;
		}	
	}
	
	
	elsif ($Remove =~ /n/)  {
		print "All flagged samples will be retained for analyses.\n";
	}

	else {
		print "Unrecognized input (a/s/n).  Samples identified as bad will be omitted from further analyses.\n"; 
	
		foreach my $badsample (@BadSamples)  {
			
			$SamplesToRemove{$badsample} = 1;
		}
	}	 
		 
}


if (@BadSamples)  {
	
	print "Are there any additional samples you would like to remove from the analysis? (y/n)\n";
	
	my $AdditionalRemove = <STDIN>;
	chomp($AdditionalRemove);
	
	my $AdditionalRemoveString;
	
	if ($AdditionalRemove =~ /y/i)  {
		print "Enter the names of the samples you want to remove. Separate each with a tab\n";
		$AdditionalRemoveString = <STDIN>;
		chomp($AdditionalRemoveString);
		
	        my @AdditionalRemoves = split(/\t/, $AdditionalRemoveString);

	        foreach my $AddRemove (@AdditionalRemoves)  {
	           	
	           $AddRemove =~ s/Individual//;	
		   $SamplesToRemove{$AddRemove} = 1;
	        }
	}
}

else {
	
	print "No samples were flagged as potentially bad samples.  Would like to remove any samples from further analyses anyway? (y/n)\n";
	
	my $AdditionalRemove = <STDIN>;
	chomp($AdditionalRemove);
	
	my $AdditionalRemoveString;
	
	if ($AdditionalRemove =~ /y/i)  {
		print "Enter the names of the samples you want to remove. Separate each with a tab\n";
		$AdditionalRemoveString = <STDIN>;
		chomp($AdditionalRemoveString);
		
	        my @AdditionalRemoves = split(/\t/, $AdditionalRemoveString);

	        foreach my $AddRemove (@AdditionalRemoves)  {
	           	
	           $AddRemove =~ s/Individual//;	
		   $SamplesToRemove{$AddRemove} = 1;
	        }
	}
	
	
	elsif ($AdditionalRemove =~ /n/) {
		
		print "All individuals will be analyzed.\n";
		print REMOVEDSAMPLES "No bad samples identified in the dataset. All were analyzed.\n";
	
	}	

}
	




############################################################################################################################
#Get rid of samples identified as bad by the user.
#Go through array AllSampleNames and remove the names that are in SamplesToRemove hash.
#Retained individuals will be stored in array called @GoodSampleNames.  These are also printed to file out/TempFiles/GoodSamples.txt.

my @GoodSampleNames;
my %SamplesToRemoveUpdate;
my $NumberSamplesRemoved = 0;
open GOODSAMPLES, ">out/TempFiles/GoodSampleNames.txt" or die$!;

if (%SamplesToRemove)  {
	
	for (keys %SamplesToRemove)  {
		$_ =~ s/"//g;
		$_ =~ s/Individual//;
		
		$SamplesToRemoveUpdate{$_} = 1;
	}	
		

	foreach my $name (@AllSampleNames)  {
		  if (exists($SamplesToRemoveUpdate{$name})) {
			  $NumberSamplesRemoved++;	
		  	  next;
		  }
		   
		  else {	
			push (@GoodSampleNames, $name);
			print GOODSAMPLES "$name\n";
		  }
	}

	
	print "The samples that have been eliminated from further analyses are...\n";
	
	print REMOVEDSAMPLES "You chose to remove $NumberSamplesRemoved of the $TotalNumberOfIndividuals samples prior to genotyping.\n ";

	print REMOVEDSAMPLES "The samples removed from analyses were...\n";
  
	for (keys %SamplesToRemoveUpdate)  {
		print "$_\n";
		print REMOVEDSAMPLES "$_\n";
	}

}


else  {
	
	foreach my $name (@AllSampleNames)  {
		push (@GoodSampleNames, $name);
		print GOODSAMPLES "$name\n";
	}
	
	print "All samples were included in analyses.\n";
	print REMOVEDSAMPLES "No samples were removed from the analyses prior to genotyping.\n";

			
}
	


close REMOVEDSAMPLES;
close GOODSAMPLES;




#######################################################################################################################
#######################################################################################################################
#Based on the retained samples, identify monomorphic loci.
#These are all of the reads stored in the ErrorTestOut file that are not part of a polymorphic locus, or a paralogous locus.
#######################################################################################################################
#######################################################################################################################

#Get all unique locus names that are polymorphic after binomial tests and not paralogous - put them in array UniqueLoci
#Will use these to identify the reads that are part of a polymorphic locus.

print "Reading SNPMatrixAll file\n";

open SNPMATRIX, "out/TempFiles/SNPMatrixAll.txt" or die$!;

my $Counter = 0;
my @UniqueLoci = ();
my %UniqueLocusNames = ();

while (<SNPMATRIX>)  {
	
	if ($Counter == 0)  {	#We're on the first line (locus names) of the file

		$_ =~ s/"//g;
		chomp($_);
		my @TotalLocusNames = split (/\t/, "$_");
		shift (@TotalLocusNames);	#Get rid of blank element at beginning
		
		
		foreach my $name (@TotalLocusNames)  {	#populate the hash UniqueLocusNames with the locus names and counts of each
			
			if (exists($UniqueLocusNames{$name})) {
			
				$UniqueLocusNames{$name}++;
			}
			
			else {			#It's a new locus name - push the location to NumbersWanted and the name to UniqueLoci.
				$UniqueLocusNames{$name} = 1;
				push (@UniqueLoci, $name);
			}
		}
		
		$Counter++;
		
	}

	
	
	else {	
		last;

	}
}

#close NODUPSOUT;
close SNPMATRIX;




##############################################################################################
#Check each locus name in the first RawReadCounts_NonParalogous file and see if it's in the array UniqueLoci.  If so, print that locus
#with it's reads to RawReadCounts_NonParalogousAllPoly file
#This is necessary because some loci in the RawReadCounts_NonParalogous file may actually be monomorphic after the genotyping was done.  
#For example, say an error read had a read depth of 10 in a single individual, at an otherwise monomorphic locus.  The read depth for
#the true sequence was 150 for this individual.  It has two alleles in the RawReadCounts_NonParalogous file, but after the binomial test,
#the error read is removed, and now it's a monomorphic locus.  
#The RawReadCounts_NonParalogousAllPoly has these loci eliminated, because they don't exist in the SNPMatrixAll file, and therefore, don't get in the array UniqueLoci.


if (-e "out/TempFiles/RawReadCountFiles/RawReadCounts_NonParalogous$GoodSampleNames[0].txt") {

	open RAWGENOSREAD, "out/TempFiles/RawReadCountFiles/RawReadCounts_NonParalogous$GoodSampleNames[0].txt" or die$!;
	open RAWGENOSWRITE, ">out/TempFiles/RawReadCounts_NonParalogousAllPoly$GoodSampleNames[0].txt" or die$!;
}

elsif (-e "out/TempFiles/RawReadCountFiles/RawReadCounts_NonParalogous$GoodSampleNames[1].txt") {	#in case the first one doesn't exist

	open RAWGENOSREAD, "out/TempFiles/RawReadCountFiles/RawReadCounts_NonParalogous$GoodSampleNames[1].txt" or die$!;
	open RAWGENOSWRITE, ">out/TempFiles/RawReadCounts_NonParalogousAllPoly$GoodSampleNames[1].txt" or die$!;
}

else {				#in case the first two don't exist - this is probably unnecessary now that samples with no data are automatically omitted from the dataset.
	open RAWGENOSREAD, "TemwpFiles/RawReadCountFiles/RawReadCounts_NonParalogous$GoodSampleNames[2].txt" or die$!;
	open RAWGENOSWRITE, ">out/TempFiles/RawReadCounts_NonParalogousAllPoly$GoodSampleNames[2].txt" or die$!;
}
	

$Starter = 0;
my $PreviousLocus;
my $LastPrinted = 0;
my @CurrentSeqsAndCounts;


while(<RAWGENOSREAD>)  {
	
	chomp($_);
	
	if ($Starter == 0)  {
		
		$PreviousLocus = $_;
		$Starter++;

	}	

	
	else  {
		
		if ($_ =~ /Locus/)  {
		
			foreach my $LocusName (@UniqueLoci)  {
				
				if ($PreviousLocus eq $LocusName)  {
					
					print RAWGENOSWRITE "$PreviousLocus\n";
					
					foreach my $line (@CurrentSeqsAndCounts)  {
						print RAWGENOSWRITE "$line\n";
					}
					
					last;
				}
				
				
			}
			
			$PreviousLocus = $_;
			$LastPrinted = 1;
			@CurrentSeqsAndCounts = ();
			
		}
		
		else {
			push (@CurrentSeqsAndCounts, $_);
			$LastPrinted = 0;
			
		}	
		
	}	

}


if ($LastPrinted == 0)  {
	
	
		foreach my $LocusName (@UniqueLoci)  {
				
				if ($PreviousLocus eq $LocusName)  {
					
					foreach my $line (@CurrentSeqsAndCounts)  {
						print RAWGENOSWRITE "$line\n";
					}
					
					last;
				}
				
				
			}
			
			$PreviousLocus = $_;
			$LastPrinted = 1;
			@CurrentSeqsAndCounts = ();


}


close RAWGENOSREAD;
close RAWGENOSWRITE;






#######################################################################################################################
#######################################################################################################################
#Identify monomorphic loci.  Create a hash from the IndividualRawGenotypeFileAllPoly file.  Keys will be the sequences.
#Search all of the reads in ErrorTestOut.txt against the keys of this hash.  
#ErrorTestOut contains all of the reads across the entire dataset that were retained as non-error reads.
#If the read don't exist in the hash, print to file MonomorphicLoci.txt.

print "Identifying monomorphic loci\n";

open MONOMORPHICLOCI, ">out/TempFiles/MonomorphicLoci.txt" or die$!;
open ERROROUT, "out/TempFiles/ErrorReadTest/ErrorTestOut.txt" or die$!;

my $FirstName = $GoodSampleNames[0];

open RAWGENOTYPES, "out/TempFiles/RawReadCounts_NonParalogousAllPoly$FirstName.txt" or die$!;

my @SeqsAndCountsForHash;

while (<RAWGENOTYPES>)  {
	
	chomp($_);
	
	if ($_ !~ /[ATGC]/)   {
		next;
	}
	
	else {
		my @TempArray = split(/\t/, $_);
		
		my $SeqToRemoveDashes = $TempArray[0];
		$SeqToRemoveDashes =~ s/-//g;
		my $SeqWODashes = $SeqToRemoveDashes;
		push (@SeqsAndCountsForHash, $SeqWODashes);
		push (@SeqsAndCountsForHash, $TempArray[1]);	
	}
}

close RAWGENOTYPES;




my %PolymorphicSeqs = @SeqsAndCountsForHash;
my $keyscounter = 0;

for my $key (keys(%PolymorphicSeqs)) {
     $keyscounter++;
}

while(<ERROROUT>) {
	chomp($_);
	if (exists $PolymorphicSeqs{$_}) {
		next;
	}
	
	else {
		print MONOMORPHICLOCI "$_\n";
	}
}


close MONOMORPHICLOCI;
close ERROROUT;


########################################################################################
#MonomorphicLoci file currently contains all reads that were part of loci identified as paralogous.  Remove all of these.

#First, put all of the paralogous reads in a hash without any dashes.

open PARALOGOUSLOCI, "out/TempFiles/ParalogousLoci.txt" or die$!;

my %Paralogs;

while (<PARALOGOUSLOCI>)  {
	
	chomp($_);
	
	if ($_ =~ /[ATGC]/)  {
		my @TempArray = split(/\t/, $_);
		my $TempRead = $TempArray[0];
		$TempRead =~ s/-//g;
		
		$Paralogs{$TempRead} = 1;
	}
}

close PARALOGOUSLOCI;

#Now, check reads in MonomorphicLoci to make sure they don't occur in Paralogs hash.  If they don't, print them to file MonomorphicsNoParalogs.

open MONOMORPHICS, "out/TempFiles/MonomorphicLoci.txt" or die$!;
open MONOMORPHICSNOPARALOGS, ">out/TempFiles/MonomorphicsNoParalogs.txt" or die$!;

while (<MONOMORPHICS>)  {
	
	chomp($_);
	
	if (exists $Paralogs{$_})  {
		
		next;
	}

	else {
		print MONOMORPHICSNOPARALOGS "$_\n";
	}	
}


close MONOMORPHICS;
close MONOMORPHICSNOPARALOGS;



########################################################################################
#Add counts for each individual to monomorphic loci file - new file called "MonomorphicLociWithCounts.txt"

my $Number = 0;
my $IndividualCounter = 0;

my $FinalNumber;

foreach my $Individual (@GoodSampleNames) {
	
	if (-e "out/TempFiles/UniqueWithCountsIndividual$Individual.txt")  {	
		$Number++;
		my $NumberPlus1 = $Number+1;
	
		if ($IndividualCounter == 0) {
		
			$IndividualCounter++;
			open CURRENTINDIVIDUAL, "out/TempFiles/UniqueWithCountsIndividual$Individual.txt" or die$!;
			my @CurrentIndividualArray = ();
		
			while (<CURRENTINDIVIDUAL>)  {
		    
				if ($_ =~ /[a-zA-Z]/) {	
					chomp($_);
					$_ =~ s/^\s*//;
					my @TempArray = split (/\s/, $_);
					push (@CurrentIndividualArray, $TempArray[1]);
					push (@CurrentIndividualArray, $TempArray[0]);
				}
			}	
	
			my %TempHash = @CurrentIndividualArray;
	
			open MONOMORPHICSEQS, "out/TempFiles/MonomorphicsNoParalogs.txt" or die$!;
			open MONOMORPHICUPDATE, ">out/TempFiles/MonomorphicUpdate$Number.txt" or die$!;
	
			while (<MONOMORPHICSEQS>) {
				chomp($_);
				my @TempArray = split (/\s+/, $_);
				my $CurrentSeq = $TempArray[0];
				my $Value = 0;
		
				if (exists($TempHash{$CurrentSeq})) {
					my $CurrentCount = $TempHash{$CurrentSeq};
					print MONOMORPHICUPDATE "@TempArray\t$CurrentCount\n";
				}
		   
				else {
					print MONOMORPHICUPDATE "@TempArray\t0\n";
				}	
			}
	
			close MONOMORPHICSEQS;
			close MONOMORPHICUPDATE;
		close CURRENTINDIVIDUAL;
	
		}	

		else {
		
			open CURRENTINDIVIDUAL, "out/TempFiles/UniqueWithCountsIndividual$Individual.txt" or die$!;
			my @CurrentIndividualArray = ();
		
			while (<CURRENTINDIVIDUAL>)  {
				if ($_ =~ /[a-zA-Z]/) {	
					chomp($_);
					$_ =~ s/^\s*//;
					my @TempArray = split (/\s/, $_);
					push (@CurrentIndividualArray, $TempArray[1]);
			    push (@CurrentIndividualArray, $TempArray[0]);
			    	}
			}	
	
			my %TempHash = ();
			%TempHash = @CurrentIndividualArray;
 
			my $NumberMinus1 = $Number-1;	
		
			open MONOMORPHICUPDATE, "out/TempFiles/MonomorphicUpdate$NumberMinus1.txt" or die$!;
			open MONOMORPHICWRITE, ">out/TempFiles/MonomorphicUpdate$Number.txt" or die$!;
			$FinalNumber = $Number;
		
			while (<MONOMORPHICUPDATE>) {
				chomp($_);
				$_ =~ s/\s/\t/;
				my @TempArray = split (/\s+/, $_);
				my $CurrentSeq = $TempArray[0];
				my $Value = 0;
		
				if (exists($TempHash{$CurrentSeq})) {
					my $CurrentCount = $TempHash{$CurrentSeq};
					print MONOMORPHICWRITE "@TempArray\t$CurrentCount\n";
				}
		
				else {
					print MONOMORPHICWRITE "@TempArray\t0\n";
				}	
			}
	
			close MONOMORPHICUPDATE;
			close MONOMORPHICWRITE;
			close CURRENTINDIVIDUAL;
	
		}
	
	}	
}




system "mv out/TempFiles/MonomorphicUpdate$FinalNumber.txt out/TempFiles/MonomorphicLociWithCounts.txt";
system "rm out/TempFiles/MonomorphicUp*";



#######################################################################################################

#Get summary of numbers of loci


my $Monomorphics;	#Line count from out/TempFiles/MonomorphicLoci.txt
my $TotalPolymorphicCountNonParalogous;		#From unique locus names in SNPMatrixAll.txt
my $ParalogousLoci;			#Line count from ParalogousLoci.txt


open MASTER, ">>out/TempFiles/MasterReport.txt" or die$!;
print MASTER "\n\n";


open PARALOGOUS, "out/TempFiles/ParalogousLoci.txt" or die$!;
my $NumParalogous = 0;

while (<PARALOGOUS>)  {
	if ($_ =~ /^Locus/) {
		$NumParalogous++;
	}
}	


print MASTER "Total number of paralogous loci identified and removed is\t$NumParalogous\n";

close PARALOGOUS;


print MASTER "Total number of monomorphic loci identified is\t";
system "wc -l out/TempFiles/MonomorphicsNoParalogs.txt | tee -a out/TempFiles/MasterReport.txt";
print MASTER "\n";



open SNPMATRIXALL, "out/TempFiles/SNPMatrixAll.txt" or die$!;

my $SNPMatrixAllCounter = 1;
my %PolymorphicNamesUnique;

while(<SNPMATRIXALL>)  {
	
	if ($SNPMatrixAllCounter==1)  {
		
		$_ =~ s/"//g;
		
		my @TempArray = split(/\t/,$_);
		
		shift(@TempArray);
		
		foreach my $name (@TempArray)  {
			
			$PolymorphicNamesUnique{$name} = 1;
		}
		
		last;
	}
	
}

my $HashSize = scalar keys %PolymorphicNamesUnique;


print MASTER "Total number of polymorphic loci identified is\t$HashSize";
print MASTER "\n";

close MASTER;
close SNPMATRIXALL;









####################################################################################################


#Update the Genotypes.txt file to include only samples in the GoodSampleNames file.

print "Updating Genotypes.txt file\n";

open GOODSAMPLES, "out/TempFiles/GoodSampleNames.txt" or die$!;

my @TotalFileSampleNamesUpdate;

while (<GOODSAMPLES>)  {
	chomp($_);
	push (@TotalFileSampleNamesUpdate, $_);
}	

close GOODSAMPLES;

my %GoodNamesHash;

foreach my $GoodName (@TotalFileSampleNamesUpdate) {
	$GoodNamesHash{$GoodName} = 1;
}	

open GENOTYPES, "out/TempFiles/Genotypes.txt" or die$!;
open GENOTYPESUPDATE, ">out/TempFiles/GenotypesUpdate.txt" or die$!;
	
while (<GENOTYPES>)  {
	chomp($_);
	
	if ($_ =~ /Locus/)  {
		print GENOTYPESUPDATE "$_\n";
	}

	else  {	
	
	   my @TempArray = split("\t", $_);
	   my $FirstElement = $TempArray[0];
	   $FirstElement =~ s/Individual//;

	   if (exists($GoodNamesHash{$FirstElement}))  {
		
		print GENOTYPESUPDATE "$_\n";
	   }
	
	}
	
}	

close GENOTYPES;
close GENOTYPESUPDATE;





#####################################################################################################
####################################################################################################
#Pull out individual SNPs and get their locations along the read.
#####################################################################################################
#####################################################################################################

print "\nIdentifying SNPs and determining the position of each SNP along the read.\n";

system "R --vanilla --slave < RScripts/OutputSNPMatrix.R";	#Outputs TempSNPMatrix.txt which has all SNPs in it and also SNPLocations.txt file, which has the position of each SNP along the read.



#####################################################################################################
#Plot frequency of SNP positions along the reads.  Give the option to remove SNPs after a certain position.


my @SNPLocationNames = ();
my @SNPLocationPositions = ();

open SNPLOCATIONS, "out/TempFiles/SNPLocations.txt" or die$!;		#This file is created by OutputSNPMatrix.R
open TEMP, ">out/TempFiles/TempLocationsToPlot.txt" or die$!;


while(<SNPLOCATIONS>)  {
	
	if ($_ =~ /^"V/)  {
		next;
	}

	elsif ($_ =~ /^"Locus"/)  {
		$_ =~ s/^"Locus" //;
		$_ =~ s/"//g;
		my @TempArray = split(/\s/,$_);
		foreach my $name (@TempArray)  {
			push (@SNPLocationNames, $name);
		}
	}

	else {
		$_ =~ s/^"SNPLocation" //;
		$_ =~ s/"//g;
		my @TempArray = split(/\s/,$_);
		foreach my $value (@TempArray)  {
			push (@SNPLocationPositions, $value);
			print TEMP "$value\t";
		}		
	}	

}


close TEMP;
close SNPLOCATIONS;



open TEMP, "out/TempFiles/TempLocationsToPlot.txt" or die$!;
open LOCATIONSTOPLOT, ">out/TempFiles/SNPLocationsToPlot.txt" or die$!;

my $LocationsCounter = 0;

while (<TEMP>)  {
	if ($LocationsCounter==0)  {
		$_ =~ s/\t$//;
		$LocationsCounter++;
		print LOCATIONSTOPLOT "$_\n";
	}

	else {	
		$_ =~ s/\t$//;
		print LOCATIONSTOPLOT "$_";
	}	
}

close TEMP;
close LOCATIONSTOPLOT;





#Now have edited the SNPLocations.txt file so that it can be read in to R and plotted.

print "Plotting SNP locations\n";

system "R --vanilla --slave < RScripts/PlotSNPLocations.R";





#####################################################################################################

print "\nLoci have all been genotyped, and scored genotypes are stored in the file 'Genotypes.txt' in the TempFiles folder.\n\n";
print "The locations of each SNP are plotted in the file SNPLocations.pdf in the out/OutputRunInfo folder.  Use this plot to choose a cut-off value for omitting SNPs due to artifacts toward the ends of reads.\n\n";
print "Next, run 'FilterSNPs.pl' to filter and output SNPs and monomorphic loci based on 1.) the proportion of samples genotyped at each locus and 2.) SNP positions along the reads.\n";



