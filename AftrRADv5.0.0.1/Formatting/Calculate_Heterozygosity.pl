#! /usr/bin/perl
use warnings;
use strict;



#Set parameters for the run

# -MinReads  The minimum number of reads required to score a monomorphic locus in an individual.
# -PopSizes  An ordered vector of population sizes for calculation of population heterozygosities.  Populations must be ordered in the haplotypes file.  Separate population sizes with a comma.  

my %RunArguments = ();
my $MinReadsFromGenotypeEdit;

open IN, "../RScripts/Genotype_Edit.R" or die$!;

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

$RunArguments{MinReads} = $MinReadsFromGenotypeEdit;
$RunArguments{PopSizes} = 0;
$RunArguments{Help} = 0;

for my $argument (@ARGV)  {
	my @TempArgs = split(/-/,$argument);
	
	#print help information
	if (($TempArgs[1] =~ /^h$/)|| ($TempArgs[1] =~ /^help$/)) {
		print "Information for Calculate_Heterozygosity.pl...\n\n";
		print "\n\n Command line arguments available in Calculate_Heterozygosity.pl...\n\n";
		
		print "MinReads\nThe minimum number of reads necessary to score a monomorphic locus in an individual.\n";
		print "Default is to use the same value used in the most recent Genotypes.pl run for genotyping polymorphic loci.\n\n";
		
		print "PopSizes\nAn ordered vector of population sizes for calculation of population heterozygosities.\n";
		print "Populations must be ordered in the haplotypes file.  Separate population sizes with a comma.\n";
		print "If this argument is given, the script will print heterozygosities for each population in addition to the defaults\n";
		print "of individual heterozygosities and overall heterozygosity\n";
		
		print "To change defaults, enter the appropriate argument and value, separated by a dash, when executing Calculate_Heteroaygosity.pl.\n";
		print "For example, the command 'perl Calculate_Heterozygosity.pl PopSizes-3,5,6' would include population heterozygosity values in the output for three populations with 3, 5, and 6 diploid samples, respectively.\n\n"; 
		exit;
		
	}	
	
	elsif (exists($RunArguments{$TempArgs[0]}))  {
		$RunArguments{$TempArgs[0]} = $TempArgs[1];
	}
}

my $MinReads = $RunArguments{MinReads};
my $OrderedPopSizes = $RunArguments{PopSizes};

my @OrderedPopSizes;


my $NumSamplesInVector = 0;

if ($RunArguments{PopSizes} ne "0")  {
	my @TempPopSizeArray = split(/,/, $OrderedPopSizes);
	for my $size (@TempPopSizeArray)  {
		if ($size =~ /[0-9]/) {
			push (@OrderedPopSizes, $size);
			$NumSamplesInVector = $NumSamplesInVector+$size;
		}
	}
}	


print "\n\n";

print "Arguments entered are...\n";
for (keys %RunArguments) {
	print "$_\t$RunArguments{$_}\n";
}

		
	
#Update the R script dealing with monomorphic loci with the MinReads parameter.

open IN, "../RScripts/Calc_Het.R" or die$!;
open OUT, ">../RScripts/Calc_Het_Edit.R" or die$!;

while(<IN>) {
	if ($_ =~ /MinReads<-edit/) {
		print OUT "MinReads<-$MinReads\n";
	}

	else {
		print OUT "$_";
	}
}

close IN;
close OUT;




#Find the correct Haplotypes file to use

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
	
	
my $LineCounter = 0;

if ($RunArguments{PopSizes} ne "0")  {
	my $NumPopsToScore = @OrderedPopSizes;
	print "Recognized $NumPopsToScore populations with a total of $NumSamplesInVector samples.\n"; 
	
	my $TotalChromosomesInHaplotypes = 0;
	
	

	open INFILE, "../Output/Genotypes/$FileName" or die$!;

	while(<INFILE>)  {

		if ($LineCounter == 0)  {
			$LineCounter++;
			next;
		}
		
		else {
			if ($_ =~ /[a-zA-Z0-9]/) {
				$TotalChromosomesInHaplotypes++;
			}
		}
	}	
			
	my $TotalSamplesInHaplotypes = $TotalChromosomesInHaplotypes/2;
	
	close INFILE;
	
	if ($TotalSamplesInHaplotypes != $NumSamplesInVector) {
		print "Warning...Number of samples accounted for in command line argument ($NumSamplesInVector) does not equal total number of samples in Haplotypes file ($TotalSamplesInHaplotypes)\n";
	}

}



	

print "\n\nRunning Calculate_Heterozygosity.pl...\n\n";



#Count heterozygous sites in the Haplotypes file

my @FirstAlleleArray = ();
my @SecondAlleleArray = ();

$LineCounter = 0;
my $NumHetSites = 0;
my $TotalScoredSites = 0;
my @SampleNames = ();

my @NumHetSitesPerSample = ();
my @TotalSitesScoredPerSample = ();

open INFILE, "../Output/Genotypes/$FileName" or die$!;

while(<INFILE>)  {

	if ($LineCounter == 0)  {
		$LineCounter++;
		next;
	}
	
	else {
		if ($LineCounter == 1) {
			@FirstAlleleArray = split(/\t/, $_);
			push(@SampleNames, $FirstAlleleArray[0]);
			shift(@FirstAlleleArray);
			my $length = @FirstAlleleArray;
			$LineCounter++;
			next;
		}	
	
		elsif ($LineCounter == 2)  {
			my @SecondAlleleArray = split(/\t/, $_);
			shift(@SecondAlleleArray);
			$LineCounter = 1;
			
			my $NumLoci = @SecondAlleleArray;
			
			for my $locus (0..$NumLoci-1) {
				
				if ($FirstAlleleArray[$locus] eq "NA") {
					next;
				}	
				
				else {
					$TotalScoredSites++;
					
					if (($FirstAlleleArray[$locus]) ne ($SecondAlleleArray[$locus]))  {
						$NumHetSites++;				
					}
				}	
			}
			
			push(@NumHetSitesPerSample, $NumHetSites);
			push(@TotalSitesScoredPerSample, $TotalScoredSites);
			$NumHetSites = 0;
			$TotalScoredSites = 0;
			
		}
	
	}
}	

close INFILE;

#print "@NumHetSitesPerSample\n";
#print "@TotalSitesScoredPerSample\n";



#Deal with monomorphic sites


my @FirstSplit = split(/_/, $FileName);
my $SecondElement = $FirstSplit[1];
my @SecondSplit = split(/\./, $SecondElement);
my $PctLociScored = $SecondSplit[0];

open MONOMORPHICINALL, "../Output/Genotypes/Monomorphics_$PctLociScored.txt" or die$!;	
open RINFILE, ">R_Mono_Infile.txt" or die$!;	

my @Monomorphics = ();

while(<MONOMORPHICINALL>)  {
	
	if ($_ =~ /Sequence/)  {
		next;
	}	
	
	if ($_ =~ /[ATGC]/)  {
		my @TempArray = split (/\s/, $_);
		shift(@TempArray);
		print RINFILE "@TempArray\n";
		
	}
}

close MONOMORPHICINALL;
close RINFILE;


system "R --vanilla --slave < ../RScripts/Calc_Het_Edit.R";

open MONOS, "MonosOut.txt" or die$!;

my @MonomorphicSitesScoredPerSample = ();

while(<MONOS>)  {
	if ($_ =~ /[0-9]/) {
		chomp($_);
		push(@MonomorphicSitesScoredPerSample, $_);
	}
}

close MONOS;

#Now have 3 arrays: @NumHetSitesPerSample, @TotalSitesScoredPerSample, and @MonomorphicSitesScoredPerSample
#The TotalSitesScoredPerSample array is only polymorphic loci, and doesn't include monomorphics.
#Each array has an element for each sample in the dataset.
#Also have sample names stored in array @SampleNames.

my $NumSamples = @MonomorphicSitesScoredPerSample;

open OUT, ">Heterozygosity_Results.txt" or die$!;

print OUT "Individual Heterozygosities...\n";

foreach my $i (0..$NumSamples-1) {
	print OUT "$SampleNames[$i]\t";
	
	my $HetValue = $NumHetSitesPerSample[$i]/($TotalSitesScoredPerSample[$i]+$MonomorphicSitesScoredPerSample[$i]);
	print OUT "$HetValue\n";
}

print OUT "\n\n";

if ($RunArguments{PopSizes} ne "0")  {
	print OUT "Population heterozygosities...\n";
	
	my $StartPosition = 0;
	my $PopNumber = 1;
	
	for my $size (@OrderedPopSizes) {
		
		my $EndPosition = $StartPosition+$size-1;
		
		my $CurrentPopHetSites = 0;
		my $CurrentPopTotalPolySitesScored = 0;
		my $CurrentPopMonoSitesScored = 0;
		
		for my $het (@NumHetSitesPerSample[$StartPosition..$EndPosition]) {
			$CurrentPopHetSites = $CurrentPopHetSites+$het;
		}
		
		for my $poly (@TotalSitesScoredPerSample[$StartPosition..$EndPosition]) {
			$CurrentPopTotalPolySitesScored = $CurrentPopTotalPolySitesScored+$poly;
		}
		
		for my $mono (@MonomorphicSitesScoredPerSample[$StartPosition..$EndPosition]) {
			$CurrentPopMonoSitesScored = $CurrentPopMonoSitesScored+$mono;
		}
		
		my $CurrentPopHetValue = $CurrentPopHetSites/($CurrentPopTotalPolySitesScored+$CurrentPopMonoSitesScored);
		
		print OUT "Population $PopNumber\t$SampleNames[$StartPosition]\tn=$size\t$CurrentPopHetValue\n";
		
		$StartPosition = $EndPosition+1;
		$PopNumber++;
		
	}

	print OUT "\n\n";
	
}

print OUT "Overall heterozygosity...\n";

my $TotalHetSites = 0;
my $TotalPolySites = 0;
my $TotalMonoSites = 0;

for my $het (@NumHetSitesPerSample)  {
	$TotalHetSites = $TotalHetSites+$het;
}

for my $poly (@TotalSitesScoredPerSample) {
	$TotalPolySites = $TotalPolySites+$poly;
}

for my $mono (@MonomorphicSitesScoredPerSample) {
	$TotalMonoSites = $TotalMonoSites+$mono;
}

my $OverallHet = $TotalHetSites/($TotalPolySites+$TotalMonoSites);

print OUT "$OverallHet";

close OUT;

print "Results printed to file Heterozygosity_Results.txt in Formatting directory\n";
