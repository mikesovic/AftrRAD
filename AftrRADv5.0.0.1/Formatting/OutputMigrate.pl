#! /usr/bin/perl
use warnings;
use strict;
use List::MoreUtils qw(uniq);
use Data::Dumper qw(Dumper);


#Prior to running OutputMigrate.pl, make sure you have a text file named "PopulationsForMigrate.txt" in the Formatting directory.
#This file should contain the word 'Population' on the first line, followed by a population designation for each sample in the order the samples appear in the file "TempFiles/GenotypesUpdate.txt".
#In the case of diploid samples, each individual should be represented by two lines in the 'PopulationsForMigrate.txt' file.

if (! -e "PopulationsForMigrate.txt") {
	print "Warning: File 'PopulationsForMigrate' not found.  This file is necessary for running OutputMigrate.pl\n";
	die;
}	

mkdir "../out/formatted_files" unless(-d "../out/formatted_files");

#This script was modified from the OutputFASTA.pl script.

my %RunArguments = ();

#Defaults
$RunArguments{SNPsOnly} = 0;
$RunArguments{unlinked} = 0;
$RunArguments{ambig} = 0;



for my $argument (@ARGV)  {
	my @TempArgs = split(/-/,$argument);
	
	#print help information
	if (($TempArgs[1] =~ /^h$/)|| ($TempArgs[1] =~ /^help$/)) {
		
		print "Information for OutputFASTA.pl...\n\n";
		print "This script creates a FASTA file, with the option of including all sites, or just variable sites.\n";  
		
		print "\n\nCommand line arguments available for this script...\n\n";
		print "SNPsOnly\nFlag to print only SNPs (variable sites) to the fasta file.\n";
		print "Set this to 1 to print SNPs only.\n";
		print "Default is 0 (all sites are included)\n";
		print "unlinked\nFlag to print only unlinked SNPs.\n";
		print "Only applies if SNPsOnly flag is set to 1.\n";
		print "Default is 0 (all variable sites are included)\n";
		print "ambig\nFlag to print only line per sample, with ambiguity codes used for heterozygous sites.\n";
		print "In the case of sites heterozygous for an indel, the base is printed, and the indel is ignored.\n";
		print "Default is 0 (two lines are printed for each sample in the fasta file.\n";
		
		exit;
		
	}
	
	#update default values with entered parameters as appropriate and continue with run
	elsif (exists($RunArguments{$TempArgs[0]}))  {
		$RunArguments{$TempArgs[0]} = $TempArgs[1];
	}
}





my $SNPsOnly = $RunArguments{SNPsOnly};
my $unlinked = $RunArguments{unlinked};
my $ambig = $RunArguments{ambig};









my $FileName;

opendir GENOS, "../out/Output/Genotypes/";
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
	print "\nEnter the name of the SNPMatrix file you want to use to create the Migrate infile\n";

	$FileName = <STDIN>;
	chomp($FileName);
}



mkdir "TempFiles" unless (-d "TempFiles");
	

open INFILE, "../out/Output/Genotypes/$FileName" or die$!;


if ($SNPsOnly == 1)  {
	
	if ($unlinked == 1)  {
		
		open RSCRIPT, "../RScripts/ID_Unlinked_SNPs.R" or die$!;
		open ROUT, ">../RScripts/ID_Unlinked_SNPs_Edit.R" or die$!;
		
		while(<RSCRIPT>)  {
			if ($_ =~ /DataTableFromFile<-read/)  {
				print ROUT "DataTableFromFile<-read.table(\"../out/Output/Genotypes/$FileName\", header=TRUE, sep = \"\\t\")\n";
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
		
		
		if ($ambig == 1) {	#print with ambiguity codes
		
			my $SNPFileCounter = 0;
			my $LineCounter = 1;
			my @CurrentAllele1Array = ();
			my @CurrentAllele2Array = ();
		
			while(<SNPFILE>)  {
				chomp($_);
				
				if ($SNPFileCounter == 0)  {
					#print SNPOUT "$_";
					$SNPFileCounter++;
					next;
				}
			
				elsif ($LineCounter == 1) {
					my @TempCurrentAllele1Array = split(/\t/, $_);
					print SNPOUT ">$TempCurrentAllele1Array[0]\n";
					shift(@TempCurrentAllele1Array);
					
					foreach my $base (@TempCurrentAllele1Array)  {
						if ($base =~ /NA/) {
							push (@CurrentAllele1Array, "N");
						}
						
						else {
							push (@CurrentAllele1Array, $base);
						}
					}	
					
					
					$LineCounter++;
					next;
				}	
				
				elsif ($LineCounter == 2) {
					my @TempCurrentAllele2Array = split(/\t/, $_);
					shift(@TempCurrentAllele2Array);
					
					foreach my $base (@TempCurrentAllele2Array)  {
						if ($base =~ /NA/) {
							push (@CurrentAllele2Array, "N");
						}
						
						else {
							push (@CurrentAllele2Array, $base);
						}
					}
					
					my @TempConsensusRead = ();
					my $Element = 0;
					
					foreach my $snp (@CurrentAllele1Array) {
						
						if ($CurrentAllele1Array[$Element] eq $CurrentAllele2Array[$Element])  {
							push (@TempConsensusRead,$CurrentAllele1Array[$Element])
						}
				
						elsif ( ($CurrentAllele1Array[$Element] eq "A" && ($CurrentAllele2Array[$Element] eq "T")) || ($CurrentAllele1Array[$Element] eq "T" && ($CurrentAllele2Array[$Element] eq "A")))  {
							push (@TempConsensusRead, "W")
						}
				
						elsif ( ($CurrentAllele1Array[$Element] eq "A" && ($CurrentAllele2Array[$Element] eq "C")) || ($CurrentAllele1Array[$Element] eq "C" && ($CurrentAllele2Array[$Element] eq "A")))  {
							push (@TempConsensusRead, "M")
						}
				
						elsif ( ($CurrentAllele1Array[$Element] eq "A" && ($CurrentAllele2Array[$Element] eq "G")) || ($CurrentAllele1Array[$Element] eq "G" && ($CurrentAllele2Array[$Element] eq "A")))  {
							push (@TempConsensusRead, "R")
						}
				
						elsif ( ($CurrentAllele1Array[$Element] eq "C" && ($CurrentAllele2Array[$Element] eq "G")) || ($CurrentAllele1Array[$Element] eq "G" && ($CurrentAllele2Array[$Element] eq "C")))  {
							push (@TempConsensusRead, "S")
						}
				
						elsif ( ($CurrentAllele1Array[$Element] eq "C" && ($CurrentAllele2Array[$Element] eq "T")) || ($CurrentAllele1Array[$Element] eq "T" && ($CurrentAllele2Array[$Element] eq "C")))  {
							push (@TempConsensusRead, "Y")
						}
				
						elsif ( ($CurrentAllele1Array[$Element] eq "G" && ($CurrentAllele2Array[$Element] eq "T")) || ($CurrentAllele1Array[$Element] eq "T" && ($CurrentAllele2Array[$Element] eq "G")))  {
							push (@TempConsensusRead, "K")
						}
						$Element++;
					}
					print SNPOUT "@TempConsensusRead";
					print SNPOUT "\n";
					$LineCounter = 1;
					@CurrentAllele1Array = ();
					@CurrentAllele2Array = ();
					
				}	
					
					
			}
		}
		
		
		
	
		else {		#don't print with ambiguity codes - print two lines per sample.
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
		}	
					
		
		
	}
	
	
	else {		#print all SNPs - not just unlinked SNPs
		
		
		open SNPFILE, "../out/Output/Genotypes/$FileName" or die$!;
		
		open SNPOUT, ">TempFiles/FASTARaw.txt" or die$!;
		
		
		if ($ambig == 1)  {
			
			my $SNPFileCounter = 0;
			my $LineCounter = 1;
			my @CurrentAllele1Array = ();
			my @CurrentAllele2Array = ();
		
			while(<SNPFILE>)  {
				chomp($_);
				
				if ($SNPFileCounter == 0)  {
					#print SNPOUT "$_";
					$SNPFileCounter++;
					next;
				}
			
				elsif ($LineCounter == 1) {
					my @TempCurrentAllele1Array = split(/\t/, $_);
					print SNPOUT ">$TempCurrentAllele1Array[0]\n";
					shift(@TempCurrentAllele1Array);
					
					foreach my $base (@TempCurrentAllele1Array)  {
						if ($base =~ /NA/) {
							push (@CurrentAllele1Array, "N");
						}
						
						else {
							push (@CurrentAllele1Array, $base);
						}
					}	
					
					$LineCounter++;
					next;
				}	
				
				elsif ($LineCounter == 2) {
					my @TempCurrentAllele2Array = split(/\t/, $_);
					shift(@TempCurrentAllele2Array);
					
					foreach my $base (@TempCurrentAllele2Array)  {
						if ($base =~ /NA/) {
							push (@CurrentAllele2Array, "N");
						}
						
						else {
							push (@CurrentAllele2Array, $base);
						}
					}
					
					
					my @TempConsensusRead = ();
					my $Element = 0;
					
					foreach my $snp (@CurrentAllele1Array) {
						
						if ($CurrentAllele1Array[$Element] eq $CurrentAllele2Array[$Element])  {
							push (@TempConsensusRead,$CurrentAllele1Array[$Element]);
						}
				
						elsif ( ($CurrentAllele1Array[$Element] eq "A" && ($CurrentAllele2Array[$Element] eq "T")) || ($CurrentAllele1Array[$Element] eq "T" && ($CurrentAllele2Array[$Element] eq "A")))  {
							push (@TempConsensusRead, "W");
						}
				
						elsif ( ($CurrentAllele1Array[$Element] eq "A" && ($CurrentAllele2Array[$Element] eq "C")) || ($CurrentAllele1Array[$Element] eq "C" && ($CurrentAllele2Array[$Element] eq "A")))  {
							push (@TempConsensusRead, "M");
						}
				
						elsif ( ($CurrentAllele1Array[$Element] eq "A" && ($CurrentAllele2Array[$Element] eq "G")) || ($CurrentAllele1Array[$Element] eq "G" && ($CurrentAllele2Array[$Element] eq "A")))  {
							push (@TempConsensusRead, "R");
						}
				
						elsif ( ($CurrentAllele1Array[$Element] eq "C" && ($CurrentAllele2Array[$Element] eq "G")) || ($CurrentAllele1Array[$Element] eq "G" && ($CurrentAllele2Array[$Element] eq "C")))  {
							push (@TempConsensusRead, "S");
						}
				
						elsif ( ($CurrentAllele1Array[$Element] eq "C" && ($CurrentAllele2Array[$Element] eq "T")) || ($CurrentAllele1Array[$Element] eq "T" && ($CurrentAllele2Array[$Element] eq "C")))  {
							push (@TempConsensusRead, "Y");
						}
				
						elsif ( ($CurrentAllele1Array[$Element] eq "G" && ($CurrentAllele2Array[$Element] eq "T")) || ($CurrentAllele1Array[$Element] eq "T" && ($CurrentAllele2Array[$Element] eq "G")))  {
							push (@TempConsensusRead, "K");
						}
						
						$Element++;
					}
					
					print SNPOUT "@TempConsensusRead";
					print SNPOUT "\n";
					$LineCounter = 1;
					@CurrentAllele1Array = ();
					@CurrentAllele2Array = ();
				}
				
				
				
					
					
			}
		}


		
		else {		#don't print with ambiguity codes - print two lines per sample.
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
		}



		
		
	}
	
	close SNPOUT;

	open SNPIN, "TempFiles/FASTARaw.txt" or die$!;
	$FileName =~ s/\.txt//;
	
	open OUTFINAL, ">$FileName.fasta" or die$!;
	
	while(<SNPIN>)  {
		
		if ($_ =~ /[A-Za-z0-9]/)  {
		
			$_ =~ s/"//g;
			$_ =~ s/\s//g;
			print OUTFINAL "$_\n";
		}	
	}
	
	
	close SNPIN;
	close OUTFINAL;
}


else {	#printing all sites, including monomorphic
	
	
	
	
	
	
	if ($ambig == 1)  {
	
		
		open INFILE, "../out/Output/Genotypes/$FileName" or die$!;
		open OUTFILE, ">TempFiles/FASTA_In_Temp.txt" or die$!;
		
		#Put locus names wanted in a hash (LocusNamesWanted).  These are the locus names in the haplotypes file.
		
		my $LineCounter = 0;
		my @LocusNamesArray;
		my @SampleNames = ();
		
		while (<INFILE>)  {
			chomp($_);
			
			if ($LineCounter == 0)  {
				@LocusNamesArray = split(/\t/, $_);
				shift(@LocusNamesArray);
				my $Length = @LocusNamesArray;
				$LineCounter++;
				
			}	
			
		
			else {
				
				if ($_ =~ /[A-Za-z0-9]/)  {
					my @TempArray = split(/\t/, $_);
					push(@SampleNames, $TempArray[0]);
				}	
			}
			
			
		}
		
		close INFILE;
		
		#my $NumIndividuals = @SampleNames/2;
		#print "Recognized $NumIndividuals individuals\n";
		
		my %LocusNamesWanted;
		
		foreach my $name (@LocusNamesArray)  {
			$LocusNamesWanted{$name} = 1;
		}
		
		
		
		
		
		
		
		
		#Put all locus names from GenotypesUpdate file in array Loci. 
		#We'll reference these against the hash of locus names from the SNPMatrix file.  
		
		open GENOTYPES, "../out/TempFiles/GenotypesUpdate.txt" or die$!;
		
		my @Loci;
		
		$LineCounter = 0;
		
		while(<GENOTYPES>)  {
			if ($LineCounter == 1)  {
				last;
			}
			
			@Loci = split(/\t/, $_);
			
			$LineCounter++;
		}
		
		pop(@Loci);
		
		close GENOTYPES;
		
		
		
		
		#Get the lengths of each locus, including indels, to account for missing data.
		
		my @LocusLengths = ("samplename");
		my @LocusLengthCheck = ("samplename");
		my $NumLociInArray = @Loci;
		
		foreach my $number  (1..$NumLociInArray) {
			
			open GENOTYPES, "../out/TempFiles/GenotypesUpdate.txt" or die$!;
		
			print "Locus number $number\n";
			
			my $printed = 0;
			
			while(<GENOTYPES>)  {
				if ($_ =~ /Locus/) {
					next;
				}	
				
				my @TempArray = split(/\t/, $_);
				if ($TempArray[$number] ne 'NA')  {
					my $TempLength = length($TempArray[$number]);
					push(@LocusLengths, $TempLength);
					print "Locus length is $TempLength\n\n";
					push(@LocusLengthCheck, $TempArray[$number]);
					$printed = 1;
					last;
				}
			}
			
			close GENOTYPES;
			
			if ($printed == 0)  {
				push(@LocusLengths, 'NA');
				push(@LocusLengthCheck, 'NA');
			}	
		
		}	
		
		my $LocusLengthsArrayLength = @LocusLengths;
		
		print "\nLength of locus lengths array is $LocusLengthsArrayLength\n\n";
		
		
		
		open GENOTYPES, "../out/TempFiles/GenotypesUpdate.txt" or die$!;
		
		
		my $CurrentFirstAlleleLine;
		my $CurrentSecondAlleleLine;
		
		my @CurrentFirstLine;
		my @CurrentSecondLine;
		
		my @ConsensusRead;
		my $CurrentSample;
		$LineCounter = 0;
		
		my $CharacterLength;
		
		my @Read1 = ();
		my @Read2 = ();
		
		
		while (<GENOTYPES>)  {
			chomp($_);
			
			if ($LineCounter == 0)  {	#Skip the locus names row
				$LineCounter++;
				next;
			}
		
		
			else {
				
				if ($LineCounter == 1)  {	#Assign first row of alleles to CurrendFirstAlleleLine
					chomp ($_);
					$LineCounter++;
					$CurrentFirstAlleleLine = $_; 
					next;
				}
				
				elsif ($LineCounter == 2)  {	#Assign second row of alleles to CurrendSecondAlleleLine
					chomp ($_);
					$LineCounter = 1;
					$CurrentSecondAlleleLine = $_;
					
				
			
					my @CurrentFirstLine = split (/\t/, $CurrentFirstAlleleLine);
					my @CurrentSecondLine = split (/\t/, $CurrentSecondAlleleLine);
					$CurrentSample = $CurrentFirstLine[0];
			
					my $ElementCounter = 0;
				
				
					my $NumLociElements = @Loci;
					my $NumCurrentFirstLineElements = @CurrentFirstLine;
				
					foreach my $CurrentLocus (@Loci)  {
					
						#my $IndelSeq = 0;
						my @TempConsensusRead = ();
						
						if (exists($LocusNamesWanted{$CurrentLocus}))  {
							
							my $CurrentAllele1 = $CurrentFirstLine[$ElementCounter];
							my $CurrentAllele2 = $CurrentSecondLine[$ElementCounter];
							
							if ($CurrentAllele1 eq "NA")  {
								my $CurrentLocusLength = $LocusLengths[$ElementCounter];
								print "Element counter $ElementCounter\t$LocusLengthCheck[$ElementCounter]\n";
								$CurrentAllele1 = ("N" x $CurrentLocusLength);
								$CurrentAllele2 = ("N" x $CurrentLocusLength);
							}
							
							my @CurrentAllele1Array = split(//,$CurrentAllele1);
							my @CurrentAllele2Array = split(//,$CurrentAllele2);
							my $ReadLength = @CurrentAllele1Array;
							
							
							foreach my $Element (0..$ReadLength-1)  {
								
								
								if ( ($CurrentAllele1Array[$Element] =~ /-/) &&  ($CurrentAllele2Array[$Element] =~ /-/))  {
									#$IndelSeq = 1;
									push (@TempConsensusRead,"-");
									
								}
								
								elsif ($CurrentAllele1Array[$Element] eq "-") {
									push (@TempConsensusRead, $CurrentAllele2Array[$Element]);
									
								}
								
								elsif ($CurrentAllele2Array[$Element] eq "-") {
									push (@TempConsensusRead, $CurrentAllele1Array[$Element]);
									
								}
								
								
								elsif ($CurrentAllele1Array[$Element] eq $CurrentAllele2Array[$Element])  {
									push (@TempConsensusRead,$CurrentAllele1Array[$Element])
								}
						
								elsif ( ($CurrentAllele1Array[$Element] eq "A" && ($CurrentAllele2Array[$Element] eq "T")) || ($CurrentAllele1Array[$Element] eq "T" && ($CurrentAllele2Array[$Element] eq "A")))  {
									push (@TempConsensusRead, "W")
								}
						
								elsif ( ($CurrentAllele1Array[$Element] eq "A" && ($CurrentAllele2Array[$Element] eq "C")) || ($CurrentAllele1Array[$Element] eq "C" && ($CurrentAllele2Array[$Element] eq "A")))  {
									push (@TempConsensusRead, "M")
								}
						
								elsif ( ($CurrentAllele1Array[$Element] eq "A" && ($CurrentAllele2Array[$Element] eq "G")) || ($CurrentAllele1Array[$Element] eq "G" && ($CurrentAllele2Array[$Element] eq "A")))  {
									push (@TempConsensusRead, "R")
								}
						
								elsif ( ($CurrentAllele1Array[$Element] eq "C" && ($CurrentAllele2Array[$Element] eq "G")) || ($CurrentAllele1Array[$Element] eq "G" && ($CurrentAllele2Array[$Element] eq "C")))  {
									push (@TempConsensusRead, "S")
								}
						
								elsif ( ($CurrentAllele1Array[$Element] eq "C" && ($CurrentAllele2Array[$Element] eq "T")) || ($CurrentAllele1Array[$Element] eq "T" && ($CurrentAllele2Array[$Element] eq "C")))  {
									push (@TempConsensusRead, "Y")
								}
						
								elsif ( ($CurrentAllele1Array[$Element] eq "G" && ($CurrentAllele2Array[$Element] eq "T")) || ($CurrentAllele1Array[$Element] eq "T" && ($CurrentAllele2Array[$Element] eq "G")))  {
									push (@TempConsensusRead, "K")
								}
									
								
								
							}
							@CurrentAllele1Array = ();
							@CurrentAllele2Array = ();
					
						}	
					
					
						$ElementCounter++;
						my $ReadWOIndels = join("", @TempConsensusRead);
						push (@ConsensusRead, $ReadWOIndels);
						
					}
					
				}
			
				
				
				print OUTFILE ">$CurrentSample\n";
				my $JoinedRead = join("", @ConsensusRead);
				print OUTFILE "$JoinedRead\n";
				$CharacterLength = length($JoinedRead);
				@CurrentFirstLine = ();
				@CurrentSecondLine = ();
				@ConsensusRead = ();
		
			}
		
				
		}
		
	
		close GENOTYPES;
		close OUTFILE;
		
		
		open INFILE, "TempFiles/FASTA_In_Temp.txt" or die$!;
		$FileName =~ s/.txt$//;
		open OUTFILE, ">$FileName.fasta" or die$!;
		
		while (<INFILE>)  {
			print OUTFILE "$_";
		}
		
		close INFILE;
		close OUTFILE;
		
		system "rm -r TempFiles";
					
					
					
					
					
		
	}	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	else  {		#printing all sites, including monomorphic, but not as ambiguous.
	
		
		open INFILE, "../out/Output/Genotypes/$FileName" or die$!;
		open OUTFILE, ">TempFiles/FASTA_In_Temp.txt" or die$!;
		
		#Put locus names wanted in a hash (LocusNamesWanted).  These are the locus names in the SNPMatrix file.
		
		my $LineCounter = 0;
		my @LocusNamesArray;
		my @SampleNames = ();
		
		while (<INFILE>)  {	# Go into the SNP Matrix file to grab the names of the loci
			chomp($_);
			
			if ($LineCounter == 0)  {
				@LocusNamesArray = split(/\t/, $_);
				shift(@LocusNamesArray);
				my $Length = @LocusNamesArray;
				$LineCounter++;
				
			}	
			
		
			else {
				
				if ($_ =~ /[A-Za-z0-9]/)  {
					my @TempArray = split(/\t/, $_);
					push(@SampleNames, $TempArray[0]);
				}	
			}
			
			
		}
		
		close INFILE;
	
		my %LocusNamesWanted;
		
		foreach my $name (@LocusNamesArray)  {
			$LocusNamesWanted{$name} = 1;
		}
		
		
		
		
		
		
		
		
		#Put all locus names from GenotypesUpdate file in array Loci. 
		#We'll reference these against the hash of locus names from the Haplotypes file. 

		## GenotypesUpdate contains the actual sequence for every locus. -- HOWEVER, we're just getting locus names here.
		
		open GENOTYPES, "../out/TempFiles/GenotypesUpdate.txt" or die$!;	
		
		my @Loci;
		
		$LineCounter = 0;
		
		while(<GENOTYPES>)  {
			if ($LineCounter == 1)  {
				last;
			}
			
			@Loci = split(/\t/, $_);
			
			$LineCounter++;
		}
		
		pop(@Loci);
		
		close GENOTYPES;
		
		
		
		#Get the lengths of each locus, including indels, to account for missing data.
		
		my @LocusLengths = ("samplename");
		my @LocusLengthCheck = ("samplename");
		my $NumLociInArray = @Loci;	### WHY IS THIS NUMBER DIFFERENT FROM THE NUMBER IN THE SNP MATRIX? 367 vs. 333/334 (because of pop)
		
		foreach my $number  (1..$NumLociInArray) {
			
			open GENOTYPES, "../out/TempFiles/GenotypesUpdate.txt" or die$!;
		
			my $printed = 0;
			
			while(<GENOTYPES>)  {
				if ($_ =~ /Locus/) {
					next;
				}	
				
				my @TempArray = split(/\t/, $_);
				if ($TempArray[$number] ne 'NA')  {
					my $TempLength = length($TempArray[$number]);
					push(@LocusLengths, $TempLength);
					push(@LocusLengthCheck, $TempArray[$number]);
					$printed = 1;
					last;
				}
			}
			
			close GENOTYPES;
			
			if ($printed == 0)  {
				push(@LocusLengths, 'NA');
				push(@LocusLengthCheck, 'NA');
			}	
		
		}	
		
		my $LocusLengthsArrayLength = @LocusLengths;
		
		
		open GENOTYPES, "../out/TempFiles/GenotypesUpdate.txt" or die$!;
		
		
		my $CurrentFirstAlleleLine;
		my $CurrentSecondAlleleLine;
		my $LocusNameForHash;
		
		my @CurrentFirstLine;
		my @CurrentSecondLine;
		
		my @ConsensusRead;
		my $CurrentSample;
		$LineCounter = 0;
		my $MigrateLineCounter = 0;	# This counter will keep a cumulative count of each line we go through in the GenotypesUpdate file. By the end it should be 2N+1... 1 being the first row of locus names from that file. We will use this to refer back to a file that lists the populations in order of how the samples appear in the GenotypesUpdate file.
		
		my $CharacterLength;
		
		my @Read1 = ();
		my @Read2 = ();
		
		my %LociHash1;
		my %LociHash2;
		my %LociHash3;
		my %LociHash4;
		my @PopsForMigrate;
		my @NamesForMigrate;
		
		open POPSFORMIGRATE, "PopulationsForMigrate.txt" or die$!;
		while (<POPSFORMIGRATE>) {
			chomp ($_);
			push (@PopsForMigrate, $_);
		}
				
		while (<GENOTYPES>)  {
			
			if ($LineCounter == 0)  {	#Skip the locus names row
				$LineCounter++;
				$MigrateLineCounter++;
				next;
			}
		
		
			else { # Now you're entering into the lines which contain 2N individuals with the alleles at each locus in every column
				
				if ($LineCounter == 1)  {	#Assign first row of alleles to CurrentFirstAlleleLine
					chomp ($_);
					$LineCounter++;
					$MigrateLineCounter++;
					# At this point we will know what row we are on at all times in order to reference what population we are dealing with for that individual 
					$CurrentFirstAlleleLine = $_; 

					next;
				}
				
				elsif ($LineCounter == 2)  {	#Assign second row of alleles to CurrendSecondAlleleLine
					chomp ($_);
					$LineCounter = 1;
					$CurrentSecondAlleleLine = $_;
					
					my @CurrentFirstLine = split (/\t/, $CurrentFirstAlleleLine);
					my @CurrentSecondLine = split (/\t/, $CurrentSecondAlleleLine);
					$CurrentSample = $CurrentFirstLine[0];	# This is where we'll know what the sample name is for the migrate file.
					push(@NamesForMigrate, $CurrentSample);
					
					my $ElementCounter = 0;	# As we go through the haplotypes file in the foreach loop below looking for whether loci were genotyped at a certain percent we will be referencing with this counter which locus we are dealing with.
					
					
					my $NumLociElements = @Loci;	# Get number of loci in Genotypes file
					my $NumCurrentFirstLineElements = @CurrentFirstLine;
					
					
					foreach my $CurrentLocus (@Loci)  {	# Now going through and grabbing every locus name for the ifexists question below.
						
									
						my $CurrentAllele1;	# Will actually be holding the allele that we'll want to put as the last value in our multidimensional hash.
						my $CurrentAllele2;	# Ditto for the other allele at said locus.
					
						if (exists($LocusNamesWanted{$CurrentLocus}))  {	# Check to see if the locus is in the Haplotypes file
							
							$CurrentAllele1 = $CurrentFirstLine[$ElementCounter];	# Grab the allele for the locus from the first line array
							$CurrentAllele2 = $CurrentSecondLine[$ElementCounter];	# Grab the other allele for the second line array... same locus
							
							if ($CurrentAllele1 eq "NA")  {
								my $CurrentLocusLength = $LocusLengths[$ElementCounter];
								$CurrentAllele1 = ("N" x $CurrentLocusLength);
								$CurrentAllele2 = ("N" x $CurrentLocusLength);
							}
							
							$LociHash1{$CurrentLocus}=$CurrentSample;
							$LociHash2{$CurrentLocus}{$CurrentSample}=$CurrentAllele1;
							$LociHash3{$CurrentLocus}{$PopsForMigrate[$MigrateLineCounter]}{$CurrentSample}=$CurrentAllele1;
							$LociHash4{$CurrentLocus}{$PopsForMigrate[$MigrateLineCounter]}{$CurrentSample}{"Allele1"}=$CurrentAllele1;
							$LociHash4{$CurrentLocus}{$PopsForMigrate[$MigrateLineCounter]}{$CurrentSample}{"Allele2"}=$CurrentAllele2;
							
							
							# print $LociHash1{$CurrentLocus};
							# print $LociHash2{$CurrentLocus}{$CurrentSample};
							# print $LociHash3{$CurrentLocus}{$PopsForMigrate[$MigrateLineCounter]}{$CurrentSample};
							# print $LociHash4{$CurrentLocus}{$PopsForMigrate[$MigrateLineCounter]}{$CurrentSample}{"Allele1"};
							# print $LociHash4{$CurrentLocus}{$PopsForMigrate[$MigrateLineCounter]}{$CurrentSample}{"Allele2"};
														
							push (@Read1, $CurrentAllele1); # This is the part you are concatenating... Probably where we'd split for Migrate
							push (@Read2, $CurrentAllele2);
						}		
						
						$ElementCounter++;
						
					}
					
					$MigrateLineCounter++;
					
					print OUTFILE ">$CurrentSample\n";
					my $JoinedRead1 = join("", @Read1);
					print OUTFILE "$JoinedRead1\n";
					
					print OUTFILE ">$CurrentSample";
					print OUTFILE "_2\n";
					my $JoinedRead2 = join("", @Read2);
					print OUTFILE "$JoinedRead2\n";
					
					@CurrentFirstLine = ();
					@CurrentSecondLine = ();
					@ConsensusRead = ();
				
				}
				
				
				@Read1 = ();
				@Read2 = ();
								
			}
		}	
			
		close GENOTYPES;
		close OUTFILE;
		
		
		###########################################################################################
		### Practice prints to make sure you're hashes are working correctly. I used these to debug and fix the large loop below. This is the syntax in which we make calls to the hashes in the various loops below... so might be wise to comment out the larger block below and see if you get the expected return values here.
		# print  $LociHash4{"LocusNumber9"}{"Killdeer"}{"IndividualOH_KLDR1006"}{"Allele2"};	# An allele
		# print length($LociHash4{"LocusNumber9"}{"Killdeer"}{"IndividualOH_KLDR1006"}{"Allele2"});	# Size of that allele.
		############################################################################################
		####### Let's make the Migrate file here...
		open MIGRATEFILE, ">../out/formatted_files/Migrate_Infile.txt";	# Create the file that we will be writing to.
		my @unique_pops = uniq(@PopsForMigrate);	# Lets grab all of the unique populations for the population file you have in the directory.
		shift(@unique_pops); # We "shift" off the first element of the array... which should just be "Population"
		my $num_of_pops = @unique_pops;	# Get the total number of populations that are used in the creation of this file.
		print MIGRATEFILE "$num_of_pops\t";	# Print that number in the first line of the migrate file.
		
		my $num_of_loci = keys %LociHash4;	# This will give you the total number of loci from the Hapltoypes file that you will be working with.
		print MIGRATEFILE "$num_of_loci\t";	# Print this on the same line.
		
		print MIGRATEFILE "Title:\n";	# Per migrate requirement we need to list some text here for the Title of the analysis.
		
		my $dummypop = $unique_pops[0];
		
		my $dummyind = $NamesForMigrate[0]; 
			
		foreach my $locus (sort keys %LociHash4) {	# Now we are going to loop over each locus.
			print MIGRATEFILE length($LociHash4{$locus}{$dummypop}{$dummyind}{"Allele1"})."\t";	# print the length of a representative allele. if you ran in Aftrrad then they should all be the same length for a locus.
		}
		print MIGRATEFILE "\n";	# End that line.
		
		
		foreach my $pop (@unique_pops) {	#	Loop over every population
			
			foreach my $locus (sort keys %LociHash4) {	# First, we need to list the number of individuals genotyped at each locus within the population. As we loop through every locus we will take the scalar value of the keys of that hash.
				my $ind_num_for_loop = scalar keys (%{$LociHash4{$locus}{$pop}});
				print MIGRATEFILE ($ind_num_for_loop*2)."\t";	# Get the number of values from the internal population hash.
			}
			print MIGRATEFILE "$pop\n";	# Print the population at the end of all the n values and next line it.
			foreach my $locus (sort keys %LociHash4) {	# Again, loop over the loci.
				foreach my $individual (@NamesForMigrate) {	# Now loop over each individual.
					if (exists($LociHash4{$locus}{$pop}{$individual})) {	# If that individuals exists as a value within the population hash then start printing the name of individual and allele for the locus you are on.  
						print MIGRATEFILE substr($individual,-9)."\t";
						print MIGRATEFILE $LociHash4{$locus}{$pop}{$individual}{"Allele1"}."\n";
						print MIGRATEFILE substr($individual,-9)."\t";
						print MIGRATEFILE $LociHash4{$locus}{$pop}{$individual}{"Allele2"}."\n";
					}
				# Then you'll go to the next individual and so forth.
				}
			# Then you'll go to the next locus.
			}
		}
		#########################################################################################
		#########################################################################################

		
		
		system "rm -r TempFiles";
			
	}
	
}

print "Completed OutputMigrate.pl.  A Migrate datafile named 'Migrate_Infile.txt' has been printed to the out/formatted_files directory.\n";
print "This file is in Sequence data format for Migrate\n";
