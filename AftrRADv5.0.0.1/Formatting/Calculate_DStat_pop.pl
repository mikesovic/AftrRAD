#! /usr/bin/perl
use warnings;
use strict;



#Set parameters for the run

my %RunArguments = ();

$RunArguments{bootstrap} = 100;


for my $argument (@ARGV)  {
	my @TempArgs = split(/-/,$argument);
	
	#print help information
	if (($TempArgs[1] =~ /^h$/)|| ($TempArgs[1] =~ /^help$/)) {
		
		print "Information for Calculate_DStat_pop.pl...\n\n";
		print "This script calculates a D statistic from population (allele frequency) data, as described in Equation 2 of Durand et al 2011.\n";  
			
		print "Prior to running this script, make sure...\n";
		print "1.) Samples from 4 populatins are grouped by population in the SNPMatrix file you will use\n";
		print "2.) The order of the populations in this file is P1, P2, P3, P4 (P4 = outgroup).\n";
		print "2.) You know the number of chromosomes sampled in each population\n";
		print "3.) The R library 'seqinr' is installed and available to load\n";
		
		print "\n\nCommand line arguments available for this script...\n\n";
		print "bootstrap\nNumber of bootstrap replicates used to generate a p-value.\n";
		print "Set this to 0 for no bootstrapping.\n";
		print "Default is 100.\n";
		
		exit;
		
	}
	
	#update default values with entered parameters as appropriate and continue with run
	elsif (exists($RunArguments{$TempArgs[0]}))  {
		$RunArguments{$TempArgs[0]} = $TempArgs[1];
	}
}





my $NumBootstraps = $RunArguments{bootstrap};

my $FileName;

opendir GENOS, "../Output/Genotypes/";
my @AllFileNames = grep { $_ ne '.' && $_ ne '..' && $_ ne '.DS_Store' } readdir(GENOS);
close GENOS;

my @HaplotypeFileNames = ();

for my $name (@AllFileNames)  {
	if ($name =~ /SNPMatrix/) {
		push (@HaplotypeFileNames, $name);
	}
}

my $NumHaplotypeFiles = @HaplotypeFileNames;

if ($NumHaplotypeFiles == 1)  {
	$FileName = $HaplotypeFileNames[0];
}

else {
	print "The following SNPMatrix files are available in the Output/Genotypes directory...\n";
	for my $file (@HaplotypeFileNames) {
		print "$file\n";
	}
	print "\nEnter the name of the SNPMatrix file you want to use for the ABBA_BABA test\n";

	$FileName = <STDIN>;
	chomp($FileName);
}


print "Enter the number of samples (chromosomes) from each population, in the order the populations occur in the SNPMatrix file. Separate with tabs.\n";
my $SampSizes = <STDIN>;
chomp($SampSizes);

my @SampleSizes = split(/\t/,$SampSizes);




mkdir "TempFiles" unless (-d "TempFiles");
	

#open INFILE, "../Output/Genotypes/$FileName" or die$!;
#open OUTFILE, ">TempFiles/FASTA_In_Temp.txt" or die$!;


	
open RSCRIPT, "../RScripts/ID_Unlinked_SNPs.R" or die$!;
open ROUT, ">../RScripts/ID_Unlinked_SNPs_Edit.R" or die$!;
		
while(<RSCRIPT>)  {
	if ($_ =~ /DataTableFromFile<-read/)  {
		print ROUT "DataTableFromFile<-read.table(\"../Output/Genotypes/$FileName\", header=TRUE, sep = \"\\t\")\n";
	}
			
	else {
		print ROUT "$_";
	}
}
		
close RSCRIPT;
close ROUT;
		
system "R --vanilla --slave < ../RScripts/ID_Unlinked_SNPs_Edit.R";
		
		
open SNPFILE, "TempFiles/UnlinkedSNPsRaw.txt" or die$!;
open SNPEDIT, ">TempFiles/UnlinkedSNPsEdit.txt" or die$!;
		
while(<SNPFILE>)  {
	$_ =~ s/"//g;
	print SNPEDIT "$_";
}
		
close SNPFILE;
close SNPEDIT;
		
		
open SNPFILE, "TempFiles/UnlinkedSNPsEdit.txt" or die$!;
		
open SNPOUT, ">TempFiles/FASTARaw.txt" or die$!;
		
		
	
my $SNPFileCounter = 0;
my $LineCounter = 1;
			
while(<SNPFILE>)  {
	chomp($_);
				
	if ($SNPFileCounter == 0)  {
		$SNPFileCounter++;
		next;
	}
				
	elsif  ($LineCounter == 1) {
					
		my @TempArray = split(/\t/, $_);
		print SNPOUT ">$TempArray[0]\n";
		$_ =~ s/NA/N/g;
		@TempArray = split(/\t/, $_);
		shift(@TempArray);	#Get rid of sample name at beginning of array.
		my $SNPs = join("",@TempArray);	#join all SNPs together
					
		print SNPOUT "$SNPs\n";
		$LineCounter++;
	}
				
	elsif ($LineCounter == 2) {
					
		my @TempArray = split(/\t/, $_);
		print SNPOUT ">$TempArray[0]_2\n";
		$_ =~ s/NA/N/g;
		@TempArray = split(/\t/, $_);
		shift(@TempArray);	#Get rid of sample name at beginning of array.
		my $SNPs = join("",@TempArray);	#join all SNPs together
					
		print SNPOUT "$SNPs\n";
		$LineCounter=1;
	}
}
			
close SNPFILE;
close SNPOUT;
		


	
open INFILE, "TempFiles/FASTARaw.txt" or die$!;
$FileName =~ s/.txt$//;
open OUTFINAL, ">ABBA_BABA_$FileName.fasta" or die$!;
		
while (<INFILE>)  {
	if ($_ =~ /[A-Za-z0-9]/)  {
		$_ =~ s/"//g;
		$_ =~ s/\s//g;
		print OUTFINAL "$_\n";
	}
}
		
close INFILE;
close OUTFINAL;
		
system "rm -r TempFiles";		
		
		
		
print "Successfully created the fasta input file for population ABBA_BABA test.\n\n";


#Edit the R script with the file name and the sample sizes.


open RORIGINAL, "../RScripts/Pop_DStat_Resamp_Binomial.R" or die$!;
open REDIT, ">../RScripts/Pop_DStat_Resamp_Binomial_Edit.R" or die$!;

while(<RORIGINAL>) {
	if ($_ =~ /^SampleSizesVector<-c/) {
		print REDIT "SampleSizesVector<-c(";
		
		for my $pop (0..2) {
			print REDIT "$SampleSizes[$pop],";
		}
		
		print REDIT "$SampleSizes[3])\n";
	}

	elsif ($_ =~ /read.alignment\(file=/) {
		print REDIT "alignment <- read.alignment(file=\"ABBA_BABA_$FileName.fasta\", format = \"fasta\")\n";
	}
	
	elsif ($_ =~ /for \(k in 1:100\)/) {
		print REDIT "for(k in 1:$NumBootstraps){\n";
	}

	elsif ($_ =~ /^NumBootstrapsFromPerl/) {
		print REDIT "NumBootstrapsFromPerl<-$NumBootstraps\n";
		
	}	

	else {
		print REDIT "$_";
	}
}


close RORIGINAL;
close REDIT;

system "R --vanilla --slave < ../RScripts/Pop_DStat_Resamp_Binomial_Edit.R";

print "\nResults printed to file DStat_pop_out.txt in Formatting directory.\n";
