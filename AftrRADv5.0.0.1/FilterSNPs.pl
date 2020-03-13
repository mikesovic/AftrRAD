#! /usr/bin/perl
use warnings;
use strict;



#Set parameters for the run

# -pctScored  Percent of individuals that must be genotyped in order to retain the locus (1-100).
# -maxSNP   The maximum location along the reads to score SNPs.  Use the plot in RunInfo/SNPLocations.pdf to choose this value.
# -MinReads  The minimum number of reads required to score a monomorphic locus in an individual.

my %RunArguments = ();
my $MinReadsFromGenotypeEdit;

open IN, "RScripts/Genotype_Edit.R" or die$!;

while(<IN>)  {
	if ($_ =~ /NumTrials</)  {
		$_ =~ s/</q/;
		$_ =~ s/\)/q/;
	
	

	
		my @TempArray = split(/q/, $_);
	
		$MinReadsFromGenotypeEdit = $TempArray[1];
		last;
	}
}	

close IN;


$RunArguments{pctScored} = 100;
$RunArguments{maxSNP} = 0;
$RunArguments{MinReads} = $MinReadsFromGenotypeEdit;
$RunArguments{Help} = 0;

for my $argument (@ARGV)  {
	my @TempArgs = split(/-/,$argument);
	
	#print help information
	if (($TempArgs[1] =~ /^h$/)|| ($TempArgs[1] =~ /^help$/)) {
		print "Information for FilterSNPs.pl...\n\n";
		print "FilterSNPs.pl should be run after the Genotype.pl script has completed.\n";
		print "FilterSNPs.pl outputs three major files each time it is run:  Monomorphics_X.txt, SNPMatrix_X.Y.txt and Haplotypes_X.Y.txt\n"; 
		print "The X and Y in the output filenames refer to the percent of samples that must be genotyped to retain a locus (pctScored), and the maximum location along the read that SNPs are retained (maxSNPS), respectively.\n\n";
		
		print "\n\n Command line arguments available in FilterSNPs.pl...\n\n";
		
		print "pctScored\nPercent of individuals that must be genotyped in order to retain the locus (1-100).\n";
		print "Default is 100 (all samples must be genotyped to retain the locus.\n\n";
		print "maxSNP\nThe maximum location along the reads to score SNPs.\n";
		print "Use the plot in RunInfo/SNPLocations.pdf to choose this value.\n";
		print "Default is '0', which prints all SNPs.\n\n";
		print "MinReads\nThe minimum number of reads necessary to score a monomorphic locus in an individual.\n";
		print "Default is to use the same value used in the most recent Genotypes.pl run for genotyping polymorphic loci.\n\n";
		
		print "To change defaults, enter the appropriate argument and value, separated by a dash, when executing FilterSNPs.pl.\n";
		print "For example, the command 'perl FilterSNPs.pl pctScored-90 maxSNP-75' would output loci (monomorphic and polymorphic) scored in at least 90\% of the samples and would only output SNPs located in the first 75 bases of the reads.\n\n"; 
		exit;
		
	}	
	
	elsif (exists($RunArguments{$TempArgs[0]}))  {
		$RunArguments{$TempArgs[0]} = $TempArgs[1];
	}
}

my $MinPercent = $RunArguments{pctScored};
my $MaxLocation = $RunArguments{maxSNP};
my $MinReads = $RunArguments{MinReads};


print "\n\n";

print "Arguments entered are...\n";
for (keys %RunArguments) {
	print "$_\t$RunArguments{$_}\n";
}


print "\n\nRunning FilterSNPs.pl...\n\n";


if ($MaxLocation == 0)  {
	print "Identifying all SNPs that are scored in $MinPercent\% of the samples\n";
}

else {
	print "Identifying SNPs that are located in the first $MaxLocation positions along the read and are scored in $MinPercent\% of the samples.\n"; 
}




mkdir "out/Output/Genotypes" unless(-d "out/Output/Genotypes");

#####################################################################################################

#Get SNPs that occur before the threshold position.

my @SNPLocationNames = ();
my @SNPLocationPositions = ();

open SNPLOCATIONS, "out/TempFiles/SNPLocations.txt" or die$!;		#This file is created by OutputSNPMatrix.R


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
			#print TEMP "$value\t";
		}		
	}	

}


close SNPLOCATIONS;



my $MinProportion = $MinPercent/100;


if ($MaxLocation > 0)  {
	
	
	open RSCRIPT, "RScripts/MinNumber.R" or die$!;
	open RSCRIPTOUT, ">RScripts/MinNumberEdit.R" or die$!;

	while (<RSCRIPT>)  {
		#replace the line in this script with the entered threshold
		chomp($_);
	
		if ($_ =~ /^MinNumberAllelesGenotypedToRetainLocus/)  {
			print RSCRIPTOUT "MinNumberAllelesGenotypedToRetainLocus<-NumRowsInMatrix*$MinProportion\n";
		}
	
		elsif ($_ =~ /write.table\(SNPMatrixWithLociWith/)  {
				print RSCRIPTOUT "write.table(SNPMatrixWithLociWithMinGenotypes, file = \"out/Output/Genotypes/SNPMatrix_$MinPercent.$MaxLocation.Temp1.txt\", sep = \"\\t\",row.names=FALSE, col.names=FALSE)";
		}
		
		
	
		else {
			print RSCRIPTOUT "$_\n";
		}	
	}	

	close RSCRIPT;
	close RSCRIPTOUT;
	
	
	
	
	
	

	my @SitesToRetain = ();
	my $CurrentElementCounter = 0;

	foreach my $element (@SNPLocationPositions)  {
		$CurrentElementCounter++;
		if ($element < $MaxLocation)  {
			push (@SitesToRetain, $CurrentElementCounter)
		}
	}

	open SITESTOKEEP, ">out/TempFiles/SitesToRetain.txt" or die$!;

	foreach my $site (@SitesToRetain)  {
		print SITESTOKEEP "$site\t";
	}

	close SITESTOKEEP;

	system "R --vanilla --slave < RScripts/OutputSNPGoodLocations.R";
}

else {
	system "cp out/TempFiles/TempSNPMatrix.txt out/TempFiles/SNPMatrix_GoodLocations.txt";
	
	open RSCRIPT, "RScripts/MinNumber.R" or die$!;
	open RSCRIPTOUT, ">RScripts/MinNumberEdit.R" or die$!;

	while (<RSCRIPT>)  {
		#replace the line in this script with the entered threshold
		chomp($_);
	
		if ($_ =~ /^MinNumberAllelesGenotypedToRetainLocus/)  {
			print RSCRIPTOUT "MinNumberAllelesGenotypedToRetainLocus<-NumRowsInMatrix*$MinProportion\n";
		}
	
		elsif ($_ =~ /write.table\(SNPMatrixWithLociWith/)  {
				print RSCRIPTOUT "write.table(SNPMatrixWithLociWithMinGenotypes, file = \"out/Output/Genotypes/SNPMatrix_$MinPercent.All.Temp1.txt\", sep = \"\\t\",row.names=FALSE, col.names=FALSE)";
		}
	
		else {
			print RSCRIPTOUT "$_\n";
		}	
	}	

	close RSCRIPT;
	close RSCRIPTOUT;
}	

#####################################################################################################
#####################################################################################################

#Apply filter based on proportion of individuals that must be genotyped at each locus.

system "R --vanilla --slave < RScripts/MinNumberEdit.R";


###########################################################################################
#Remove dups in the SNPMatrix file and create haplotypes.

if ($MaxLocation == 0)  {
	open SNPMATRIXTEMP, "out/Output/Genotypes/SNPMatrix_$MinPercent.All.Temp1.txt" or die$!;
	
	open SNPMATRIX, ">out/Output/Genotypes/SNPMatrix_$MinPercent.All.Temp.txt" or die$!;
	
	my $Counter = 0;
	
	while(<SNPMATRIXTEMP>)  {
		if ($Counter == 0)  {
			my @TempArray = split(/\t/,$_);
			my $LocusNamesString = "\t";
			foreach my $name (@TempArray)  {
				if ($name !~ /L/)  {
					next;
				}
				
				else {
					$name =~ s/\..*/"/;
					$LocusNamesString = $LocusNamesString.$name;
					$LocusNamesString = $LocusNamesString."\t";
				}
			}
			$LocusNamesString =~ s/\t$//;
			print SNPMATRIX "$LocusNamesString";
			$Counter++;
		}
		
		else {
			print SNPMATRIX "$_";
		}
	}	
close SNPMATRIXTEMP;
close SNPMATRIX;
system "rm out/Output/Genotypes/SNPMatrix_$MinPercent.All.Temp1.txt";

}




else {
	open SNPMATRIXTEMP, "out/Output/Genotypes/SNPMatrix_$MinPercent.$MaxLocation.Temp1.txt" or die$!;
	
	open SNPMATRIX, ">out/Output/Genotypes/SNPMatrix_$MinPercent.$MaxLocation.Temp.txt" or die$!;
	
	my $Counter = 0;
	
	while(<SNPMATRIXTEMP>)  {
		if ($Counter == 0)  {
			my @TempArray = split(/\t/,$_);
			my $LocusNamesString = "\t";
			foreach my $name (@TempArray)  {
				if ($name !~ /L/)  {
					next;
				}
				
				else {
					$name =~ s/\..*/"/;
					$LocusNamesString = $LocusNamesString.$name;
					$LocusNamesString = $LocusNamesString."\t";
				}
			}
			$LocusNamesString =~ s/\t$//;
			print SNPMATRIX "$LocusNamesString";
			$Counter++;
		}
		
		else {
			print SNPMATRIX "$_";
		}
	}
	
close SNPMATRIXTEMP;
close SNPMATRIX;
system "rm out/Output/Genotypes/SNPMatrix_$MinPercent.$MaxLocation.Temp1.txt";

}






open NODUPSOUT, ">out/TempFiles/NoDupsOutHaplotypes.txt" or die$!;

if ($MaxLocation == 0)  {
	
	open SNPMATRIX, "out/Output/Genotypes/SNPMatrix_$MinPercent.All.Temp.txt" or die$!;
}

else {
	
	open SNPMATRIX, "out/Output/Genotypes/SNPMatrix_$MinPercent.$MaxLocation.Temp.txt" or die$!;
}	
	

my $Counter = 0;
my @NumbersWanted = ();		#Positions of the first instance of each unique locus
my @UniqueLoci = ();
my %UniqueLocusNames = ();	#Contains counts of the number of times each unique is repeated.

my @TotalLocusNames = ();

while (<SNPMATRIX>)  {
	
	if ($Counter == 0)  {	#We're on the first line (locus names) of the file

		my $CurrentElement = 0;
		$_ =~ s/"//g;
		chomp($_);
		@TotalLocusNames = split (/\t/, "$_");
		shift (@TotalLocusNames);	#Get rid of blank element at beginning
		
		foreach my $name (@TotalLocusNames)  {	#populate the hash UniqueLocusNames with the locus names and counts of each
			
			if (exists($UniqueLocusNames{$name})) {
			
				$UniqueLocusNames{$name}++;
			}
			
			else {			#It's a new locus name - push the location to NumbersWanted and the name to UniqueLoci.
				$UniqueLocusNames{$name} = 1;
				push (@NumbersWanted, $CurrentElement);
				push (@UniqueLoci, $name);
			}
			$CurrentElement++;
		}
		$Counter++;
		
	


		for my $key (keys(%UniqueLocusNames))  {	#subtract 1 from each value in the hash, to have the number of times the unique locus is repeated
	
			my $NewCount = $UniqueLocusNames{$key}-1;
	
			$UniqueLocusNames{$key} = $NewCount;
	
		}
		
		foreach my $LocusName (@UniqueLoci)  {	#Create the header line (names) in NoDupsOutHaplotypesForMonoTest.txt
			
			print NODUPSOUT "\t$LocusName";
		}
		
		my $NumInNumbersWanted = @NumbersWanted;
		my $NumInUniqueLoci = @UniqueLoci;
		
	}

	
	
	else {	 #We're on a line of SNP matrix that has genotypes (not the header line)
		chomp($_);
		$_ =~ s/"//g;
		my @TempArray = split (/\t/, "$_");
		
		print NODUPSOUT "\n$TempArray[0]";	#print the sample name
		
		foreach my $Number (@NumbersWanted)  {
			
			my $CurrentUniqueLocusName = $TotalLocusNames[$Number];
			
			if ($UniqueLocusNames{$CurrentUniqueLocusName} == 0)  {			#if the count associated with that number is zero, it's a single - print it
			
			   print NODUPSOUT "\t$TempArray[$Number+1]";
			}
			
			else {
				print NODUPSOUT "\t$TempArray[$Number+1]";
				my $TempNumberWanted = $UniqueLocusNames{$CurrentUniqueLocusName};
				
				foreach my $number (1..$TempNumberWanted)  {
					print NODUPSOUT "$TempArray[$Number+1+$number]";
				}
			}	
		}
	}
}


close NODUPSOUT;
close SNPMATRIX;


open NODUPS, "out/TempFiles/NoDupsOutHaplotypes.txt" or die$!;

if ($MaxLocation == 0)  {
	open NODUPSUPDATED, ">out/Output/Genotypes/Haplotypes_$MinPercent.All.txt" or die$!;
}

else {
	open NODUPSUPDATED, ">out/Output/Genotypes/Haplotypes_$MinPercent.$MaxLocation.txt" or die$!;
}	

while (<NODUPS>)  {

	$_ =~ s/(NA)+/NA/g;
	
print NODUPSUPDATED "$_";	
}	


close NODUPS;
close NODUPSUPDATED;




if ($MaxLocation == 0)  {
	
	open SNPMATRIX, "out/Output/Genotypes/SNPMatrix_$MinPercent.All.Temp.txt" or die$!;
	open OUT, ">out/Output/Genotypes/SNPMatrix_$MinPercent.All.txt" or die$!;
	
		while (<SNPMATRIX>)  {
			$_ =~ s/"//g;
			print OUT "$_";
		}	
}

else {
	
	open SNPMATRIX, "out/Output/Genotypes/SNPMatrix_$MinPercent.$MaxLocation.Temp.txt" or die$!;
	open OUT, ">out/Output/Genotypes/SNPMatrix_$MinPercent.$MaxLocation.txt" or die$!;
	
		while (<SNPMATRIX>)  {
			$_ =~ s/"//g;
			print OUT "$_";
		}	

}

close SNPMATRIX;
close OUT;


system "rm out/Output/Genotypes/SNPMatrix_$MinPercent*Tem*.txt";


#Output monomorphic loci scored in the threshold proportion of samples.

open MONOMORPHICSWITHCOUNTS, "out/TempFiles/MonomorphicLociWithCounts.txt" or die$!;
open MONOMORPHICOUT, ">out/Output/Genotypes/Monomorphics_$MinPercent.txt" or die$!;

#Print the sample names as column headings in the file

my @AllGoodSampleNames;

open GOODSAMPLES, "out/TempFiles/GoodSampleNames.txt" or die$!;

while(<GOODSAMPLES>)  {
	if ($_ =~ /[A-Za-z0-9]/) {
		chomp($_);
		push(@AllGoodSampleNames, $_);	
	}
}

close GOODSAMPLES;

print MONOMORPHICOUT "Sequence";

for my $name (@AllGoodSampleNames)  {
	print MONOMORPHICOUT "\t$name";
}

print MONOMORPHICOUT "\n";


my $ThresholdToCall;

while (<MONOMORPHICSWITHCOUNTS>)  {
	
	my $NumScored = 0;
	
	chomp($_);
	my @TempArray = split (/\s+/, $_);
	my $CurrentNumElements = @TempArray;
	my $NumSamples = $CurrentNumElements-1;
	
	$ThresholdToCall = int($NumSamples*$MinPercent/100);
	
	foreach my $value (@TempArray[1..$CurrentNumElements-1])  {
		
		if ($value >= $MinReads)  {
			$NumScored++;
		}
	}
	
	if ($NumScored >= $ThresholdToCall)  {
		print MONOMORPHICOUT "$_\n";
	}	
	
}	

close MONOMORPHICOUT;
close MONOMORPHICSWITHCOUNTS;


my $FirstName;

open GOODSAMPLES, "out/TempFiles/GoodSampleNames.txt" or die$!;

while(<GOODSAMPLES>)  {
	chomp($_);
	$FirstName = $_;
	last;
}	



#Get mean/median read counts per individual/locus.


#First, create hash (AllSeqs) that contains all reads from one of the NonParalogousRawGenotypes files, then add all monomorphic seq reads scored to this hash.


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

my %AllSeqs = %PolymorphicSeqs;



open MONOMORPHICS, "out/Output/Genotypes/Monomorphics_$MinPercent.txt" or die$!;


while (<MONOMORPHICS>)  {
	chomp($_);
	
	if ($_ =~ /[ATGC]/) {
	    
		my @TempArray = split(/\t/, $_);
		my $CurrentSeq = $TempArray[0];
		
		$AllSeqs{$CurrentSeq} = 1;
	}	
	 
		

}	
close MONOMORPHICS;





#For each line in AllReadsAndDepths, see if the seq is in the AllSeqs hash (which doesn't include paralogous loci).  If so, push all of the non-zero counts associated with the seq to an array.
#There are a couple of problems here - first, it's not accounting for individuals that were removed from the analysis, as they are still in the AllReadsAndDepths file.  Maybe print a new file that doesn't include these?
#Also, this counts loci that might not have been genotyped at all, because the read counts were too low.  Could reference against the locus names of the Haplotypes file (maybe should do this?)


open CURRENTBARCODEFILE, "out/TempFiles/MasterBarcodeFile.txt" or die$!;

my $TotalNumberOfIndividuals = 0;
my @AllSampleNames;
my @OmittedSamples;


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

close CURRENTBARCODEFILE;




my @TotalCounts = ();
my $SampleSize = 0;
my $Sum = 0;

open VARFILE, "out/TempFiles/ErrorReadTest/AllReadsAndDepths.txt" or die$!;

while (<VARFILE>)  {
	chomp($_);
	my @TempArray = split(/\t/, $_);
	my $CurrentSeq = $TempArray[0];
	#print "$CurrentSeq.55";
	
	if (exists($AllSeqs{$CurrentSeq}))  {
		
		foreach my $a (@TempArray[1..$TotalNumberOfIndividuals])  {
			if ($a == 0)  {
				next;
			}
			
			else {
				push (@TotalCounts, $a);
				$SampleSize++;
			}
		}
	}
}


#Sum the values in the array and divide by the total length of the array.
#Also, print values to file to plot in R.

open TEMP, ">out/TempFiles/TempDepthsToPlot.txt" or die$!;

foreach my $count (@TotalCounts)  {
	
	print TEMP "$count\t";
	$Sum = $Sum+$count;
}

close TEMP;

open TEMP, "out/TempFiles/TempDepthsToPlot.txt" or die$!;
open DEPTHS, ">out/TempFiles/DepthsToPlot.txt" or die$!;

while (<TEMP>)  {
	$_ =~ s/\t$//;
	print DEPTHS "$_";
}

close TEMP;
close DEPTHS;
	

my $average = $Sum/$SampleSize;


#Print results to MasterReport file.

system "cp out/TempFiles/MasterReport.txt out/Output/RunInfo/MasterReportTemp.txt";

open TEMP, "out/Output/RunInfo/MasterReportTemp.txt" or die$!;
open OUT, ">out/Output/RunInfo/MasterReport.txt" or die$!;


while(<TEMP>)  {
	
	$_ =~ s/TempFiles\/MonomorphicsNoParalogs.txt//;
	
	if ($_ =~ /depth/)  {
		next;
	}

	elsif ($_ =~ /paralogous/) {
		print OUT "\n$_";
	
		
	}
	
	elsif ($_ =~ /[A-Za-z0-9]/) {
		
		print OUT "$_";
	}	
	
}	

close TEMP;


print OUT "\nThe average read depth across all polymorphic and monomorphic loci is $average\n";

my @SortedCounts = sort {$a <=> $b} @TotalCounts;

my $Median;

if ($SampleSize%2) {
	$Median=$SortedCounts[int($SampleSize/2)];
}

else {
	$Median = ($SortedCounts[($SampleSize/2)-1] + $SortedCounts[($SampleSize/2)-1])/2;
}	


print OUT "The median read depth across all polymorphic and monomorphic loci is $Median\n\n";



my $MinCount = $SortedCounts[0];
my $MaxCount = $SortedCounts[-1];

#print MASTER "The minimum read count at a locus was $MinCount and the maximum count was $MaxCount\n"; 

close OUT;
close VARFILE;


#Plot the read counts

#Edit PlotSNPLocations.R script to set x-axis

open RFILE, "RScripts/PlotReadCounts.R" or die$!;
open OUTFILE, ">RScripts/PlotReadCountsEdit.R" or die$!;

while(<RFILE>) {
	if ($_ =~ /^MaxLength/)  {
		print OUTFILE "MaxLength<-$MaxCount\n";
	}
	
	else {
		print OUTFILE "$_";
	}
}

close RFILE;
close OUTFILE;

system "rm out/Output/RunInfo/MasterReportTemp.txt";


#system "R --vanilla --slave < RScripts/PlotReadCountsEdit.R";

print "Individual SNPs, haplotypes, and monomorphic loci have been printed to files located in out/Output/Genotypes.\n";
print "These data can be further formatted using the scripts in the 'Formatting' directory\n";
