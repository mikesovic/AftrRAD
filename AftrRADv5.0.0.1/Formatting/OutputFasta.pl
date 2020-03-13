#! /usr/bin/perl
use warnings;
use strict;



#Set parameters for the run

# -SNPsOnly  Flag to print only SNPs (variable sites) to the fasta file.
# -unlinked  Flag to print only unlinked SNPs.  Only applies if SNPsOnly flag is set to 1.
# -ambig     Flag to print one line per sample, with ambiguity codes used for heterozygous sites.

my %RunArguments = ();

#Defaults
$RunArguments{SNPsOnly} = 0;
$RunArguments{unlinked} = 0;
$RunArguments{ambig} = 0;



for my $argument (@ARGV)  {
	my @TempArgs = split(/-/,$argument);
	
	#print help information
	if (($TempArgs[1] =~ /^h$/)|| ($TempArgs[1] =~ /^help$/)) {
		
		print "Information for OutputFasta.pl...\n\n";
		print "This script creates a Fasta file, with the option of including all sites, or just variable sites.\n";  
		
		print "\n\nCommand line arguments available for this script...\n\n";
		print "SNPsOnly\nFlag to print only SNPs (variable sites) to the fasta file.\n";
		print "Set this to 1 to print SNPs only.\n";
		print "Default is 0 (all sites are included)\n";
		print "unlinked\nFlag to print only unlinked SNPs.\n";
		print "Only applies if SNPsOnly flag is set to 1.\n";
		print "Default is 0 (all variable sites are included)\n";
		print "ambig\nFlag to print only line per sample, with ambiguity codes used for heterozygous sites.\n";
		print "In the case of sites heterozygous for an indel, the base is printed, and the indel is ignored.\n";
		print "Default is 0 (two lines are printed for each sample in the fasta file).\n";
		
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
	print "\nEnter the name of the SNPMatrix file you want to use to create the Fasta infile\n";

	$FileName = <STDIN>;
	chomp($FileName);
}



mkdir "TempFiles" unless (-d "TempFiles");
	

open INFILE, "../Output/Genotypes/$FileName" or die$!;


if ($SNPsOnly == 1)  {
	
	if ($unlinked == 1)  {
		
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
		
		
		open SNPFILE, "../Output/Genotypes/$FileName" or die$!;
		
		open SNPOUT, ">TempFiles/FASTARaw.txt" or die$!;
		
		
		if ($ambig == 1)  {
			
			my $SNPFileCounter = 0;
			my $LineCounter = 1;
			my @CurrentAllele1Array = ();
			my @CurrentAllele2Array = ();
		
			while(<SNPFILE>)  {
				chomp($_);
				
				if ($SNPFileCounter == 0)  {
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
	
		
		open INFILE, "../Output/Genotypes/$FileName" or die$!;
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
		
		
		my %LocusNamesWanted;
		
		foreach my $name (@LocusNamesArray)  {
			$LocusNamesWanted{$name} = 1;
		}
		
		
		
		
		
		
		
		
		#Put all locus names from GenotypesUpdate file in array Loci. 
		#We'll reference these against the hash of locus names from the SNPMatrix file.  
		
		open GENOTYPES, "../TempFiles/GenotypesUpdate.txt" or die$!;
		
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
			
			open GENOTYPES, "../TempFiles/GenotypesUpdate.txt" or die$!;
		
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

		
		open GENOTYPES, "../TempFiles/GenotypesUpdate.txt" or die$!;
		
		
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
					
						my @TempConsensusRead = ();
						
						if (exists($LocusNamesWanted{$CurrentLocus}))  {
							
							my $CurrentAllele1 = $CurrentFirstLine[$ElementCounter];
							my $CurrentAllele2 = $CurrentSecondLine[$ElementCounter];
							
							if ($CurrentAllele1 eq "NA")  {
								my $CurrentLocusLength = $LocusLengths[$ElementCounter];
								$CurrentAllele1 = ("N" x $CurrentLocusLength);
								$CurrentAllele2 = ("N" x $CurrentLocusLength);
							}
							
							my @CurrentAllele1Array = split(//,$CurrentAllele1);
							my @CurrentAllele2Array = split(//,$CurrentAllele2);
							my $ReadLength = @CurrentAllele1Array;
							
							
							foreach my $Element (0..$ReadLength-1)  {
								
								
								if ( ($CurrentAllele1Array[$Element] =~ /-/) &&  ($CurrentAllele2Array[$Element] =~ /-/))  {
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
	
		
		open INFILE, "../Output/Genotypes/$FileName" or die$!;
		open OUTFILE, ">TempFiles/FASTA_In_Temp.txt" or die$!;
		
		#Put locus names wanted in a hash (LocusNamesWanted).  These are the locus names in the SNPMatrix file.
		
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
		
		my %LocusNamesWanted;
		
		foreach my $name (@LocusNamesArray)  {
			$LocusNamesWanted{$name} = 1;
		}
		
		
		
		
		
		
		
		
		#Put all locus names from GenotypesUpdate file in array Loci. 
		#We'll reference these against the hash of locus names from the Haplotypes file.  
		
		open GENOTYPES, "../TempFiles/GenotypesUpdate.txt" or die$!;
		
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
			
			open GENOTYPES, "../TempFiles/GenotypesUpdate.txt" or die$!;
		
			
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
		
		
		
		open GENOTYPES, "../TempFiles/GenotypesUpdate.txt" or die$!;
		
		
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
					
					
					my $NumLociElements = @Loci;	#Get number of loci in Genotypes file
					my $NumCurrentFirstLineElements = @CurrentFirstLine;
					
					
					foreach my $CurrentLocus (@Loci)  {
									
						my $CurrentAllele1;
						my $CurrentAllele2;
					
						if (exists($LocusNamesWanted{$CurrentLocus}))  {	#Check to see if the locus is in the Haplotypes file
							
							$CurrentAllele1 = $CurrentFirstLine[$ElementCounter];
							$CurrentAllele2 = $CurrentSecondLine[$ElementCounter];
							
							if ($CurrentAllele1 eq "NA")  {
								my $CurrentLocusLength = $LocusLengths[$ElementCounter];
								$CurrentAllele1 = ("N" x $CurrentLocusLength);
								$CurrentAllele2 = ("N" x $CurrentLocusLength);
							}	
							
							push (@Read1, $CurrentAllele1);
							push (@Read2, $CurrentAllele2);
						}		
						
						$ElementCounter++;
						
						
					}
					
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
	
}	
			
			
print "File $FileName.fasta is available in Formatting directory\n";			
			
		
