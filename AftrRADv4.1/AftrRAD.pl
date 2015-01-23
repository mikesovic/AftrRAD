#! /usr/bin/perl
use warnings;
use strict;

##################################################################################################################################################################################################
#Set parameters for the run


# -re Restriction enzyme recognition site
# -minQual Minimum Phred score
# -minDepth Minimum mean nonzero read depth to retain read
# -minIden  Minimum percent identity to consider two reads allelic
# -numIndels  Maximum number of indels allowed between each pair of sequences after alignment for the reads to be assigned to the same locus
# -P2  Beginning of the P2 adaptor sequence
# -stringLength  The maximum length of a homopolymer allowed in a read (i.e. AAAAAAAAAAAAAAAAA).  Reads containing homopolymers longer than this will be removed.
# -minParalog  Minimum number of reads at a third allele in an individual to flag locus as paralogous
# -Phred  33 or 64
# -DataPath Path to directory containing fastq data files for the run
# -BarcodePath  Path to directory containing barcode files for the run
# -dplexedData  Indicator (0/1) of whether the data are already demultiplexed (separate fastq file exists for each sample to be analyzed).


my %RunArguments = ();

#Defaults
$RunArguments{re} = 'TGCAGG';
$RunArguments{minQual} = 20;
$RunArguments{minDepth} = 5;
$RunArguments{minIden} = 90;
$RunArguments{P2} = 'ATTAGATC';
$RunArguments{minParalog} = 5;
$RunArguments{Phred} = 33;
$RunArguments{DataPath} = 'Data/';
$RunArguments{BarcodePath} = 'Barcodes/';
$RunArguments{Help} = 0;
$RunArguments{numIndels} = 3;
$RunArguments{dplexedData} = 0;
$RunArguments{stringLength} = 15;
$RunArguments{MaxH} = 90;

for my $argument (@ARGV)  {
	my @TempArgs = split(/-/,$argument);
	
	#print help information
	if (($TempArgs[1] =~ /^h$/)|| ($TempArgs[1] =~ /^help$/)) {
		
		print "Information for AftrRADv3.4...\n\n";
		print "AftrRAD is a bioinformatic pipeline written in Perl and R for processing and analyzing RADseq data.\n";  
		print "It is currently available for Mac and Linux systems, but runs most efficiently on Macs. \n";
		print "More details on running AftrRAD can be found in the Getting Started folder.";
		
		
		print "\n\nCommand line arguments available in AftrRAD...\n\n";
		print "re\nThe restriction enzyme recognition sequence that occurs in your reads.\n";
		print "If there is no restriction enzyme recognition sequence in you reads, enter 0.\nDefault is TGCAGG, which corresponds to digestion with SbfI. \n\n";
		print "minQual\nMinimum quality (Phred) score for retaining reads.\n";
		print "Reads that contain bases below this value are removed from the analysis.\nDefault is 20.\n\n";
		
		print "minDepth\nThe minimum mean nonzero read depth to retain a read.\n";
		print "Reads with values less than this threshold are eliminated as error reads.\nDefault is 5.\n\n";
		print "minIden\nThe minimum percent identity (after alignment) to consider two reads to be alternative alleles from the same locus.\n";
		print "Default is 90\%.\n\n";
		print "numIndels\nThe maximum number of indels allowed between any two reads to consider the reads to be alternative alleles from the same locus.\n";
		print "Default is 3.\n\n";
		print "P2\nThe beginning of the P2 adaptor sequence.\n";
		print "Reads containing this string are removed from the analysis.\nDefault is ATTAGATC.\n\n";
		print "minParalog\nThe minimum number of reads that must occur at a third allele in an individual to flag a locus as paralogous.\n";
		print "Default is 5.\n\n";
		print "stringLength\nMaximum length of a homopolymer allowed in a read (i.e. AAAAAAAAAAAAAAA).\n";  
		print "Reads with homopolymers exceeding this length will be removed.\nDefault is 15.\n\n";
		print "Phred\nThe quality score methodology used in the sequencing.\n";
		print "Most sequencers, including Illumina HiSeq, use Phred33. An alternative is Phred64.\nDefault is Phred33\n\n";
		print "MaxH\nMaximum proportion of samples allowed to be heterozygous at a locus.\n";
		print "Loci for which heterozygosity exceeds this value will be flagged as paralogous.\nDefault is 90\%.\n\n";
		print "dplexedData\nFlag to indicate data are already demultiplexed.\n";
		print "If the data are already demultiplexed, set this to \"1\", and put all demultiplexed fastq files in a directory named \"DemultiplexedFiles\".\n";
		print "Default is \"0\", indicating data need demultiplexed.\n\n";
		print "DataPath\nPath to the directory containing fastq data files for the run.\n";
		print "The names of these fastq files must be identical to their associated barcode files.  Note that any path specified here must end with a /.\n";
		print "Default is the Data folder in the working AftrRAD directory (Data/).\n\n";
		print "BarcodePath\nPath to the directory containing barcode files for the run.\n";
		print "This directory can contain only barcode files for the current run, and each must have an associated fastq file in the data directory with an identical name.\n";
		print "Default is the Barcodes folder in the working AftrRAD directory (Barcodes/).\n\n\n";
		print "To change defaults, enter the appropriate argument and value, separated by a dash, when executing AftrRAD.pl.\n";
		print "For example, the command 'perl AftrRAD.pl re-ATCACG numIndels-4' would change the re and numIndels defaults.\n";
		exit;
		
	}
	
	#update default values with entered parameters as appropriate and continue with run
	elsif (exists($RunArguments{$TempArgs[0]}))  {
		$RunArguments{$TempArgs[0]} = $TempArgs[1];
	}
}



#Check to see if data are already demultiplexed.  If so, format the data to fit the expected format for AftrRAD.

if ($RunArguments{dplexedData} == 1) {
	
	if (-d "DemultiplexedFiles")  {
	
	
		my $QualityScore;
	
	
		if ($RunArguments{Phred} == 33)  {
			$QualityScore = "I";
		}
	
		if ($RunArguments{Phred} == 64)  {
			$QualityScore = "g";
		}
	
		#assign high Phred scores to the inline barcodes that will be assigned to each individual.
		$QualityScore = $QualityScore.$QualityScore.$QualityScore.$QualityScore.$QualityScore.$QualityScore;
	
		my @RandomBarcodes = ();
	
		my @Bases = ("A", "T", "G", "C");
		my %Barcodes = ();
	
		#generate 1500 random, unique 6 bp barcodes - this assumes the dataset contains no more than 1500 individuals, and actually probably fewer than this, due to the nature of the hash these barcodes will be stored in.
		
		foreach my $number (1..1500)  {
		
			my $random1 = int(rand(4));
			my $random2 = int(rand(4));
			my $random3 = int(rand(4));
			my $random4 = int(rand(4));
			my $random5 = int(rand(4));
			my $random6 = int(rand(4));
		
			my $CurrentBarcode = $Bases[$random1].$Bases[$random2].$Bases[$random3].$Bases[$random4].$Bases[$random5].$Bases[$random6];
		
			$Barcodes{$CurrentBarcode} = 1;
		
		}
		
		my $NumberBarcodes = 0;
	
		for(keys(%Barcodes)) {
			push(@RandomBarcodes, $_);
		}
	
	
		opendir DATA, "DemultiplexedFiles";
		my @FileNames = grep { $_ ne '.' && $_ ne '..' && $_ ne '.DS_Store' } readdir(DATA);
		close DATA;
	
		mkdir "Data" unless(-d "Data");
		mkdir "Barcodes" unless(-d "Barcodes");
		
		open OUT, ">Data/AllSamples.txt" or die$!;
		open OUTBARCODES, ">Barcodes/AllSamplesRaw.txt" or die$!;
	
		my $CurrentBarcodeNumber = 0;
	
		for my $file (@FileNames)  {
		
			my $CurrentDplexedFile = $file;
			$CurrentDplexedFile =~ s/.txt//;
			$CurrentDplexedFile =~ s/.fastq//;
		
			print "Formatting data file $file\n";
			$CurrentBarcodeNumber++;
		
			print OUTBARCODES "$RandomBarcodes[$CurrentBarcodeNumber]\t$CurrentDplexedFile\n";
		
			open FILE, "DemultiplexedFiles/$file" or die$!;
		
			my $CurrentLine = 0;
			
			while(<FILE>)  {
			
				if ($CurrentLine == 0)  {
					print OUT "$_";
					$CurrentLine++;
					next;
				}
			
				if ($CurrentLine == 1)  {
					my $Line = $RandomBarcodes[$CurrentBarcodeNumber].$_;
					print OUT "$Line";
					$CurrentLine++;
					next;
				}
			
				if ($CurrentLine == 2)  {
					print OUT "$_";
					$CurrentLine++;
					next;
				}
			
				if ($CurrentLine == 3)  {
					my $Line = $QualityScore.$_;
					print OUT "$Line";
					$CurrentLine = 0;
					next;
				}
			}
	
			close FILE;
		}
	
		
		close OUT;
		close OUTBARCODES;
	
		open BARCODES, "Barcodes/AllSamplesRaw.txt" or die$!;
		open OUT, ">Barcodes/AllSamples.txt" or die$!;
	
		my $PrintedFirst = 0;
	
		while(<BARCODES>)  {
			if ($_ =~ /[A-Za-z1-9]/)  {
				chomp($_);
			
				if ($PrintedFirst == 0)  {
					print OUT "$_";
					$PrintedFirst = 1;
				}
			
				else {
					print OUT "\n$_";
				}	
				
			}
		}
	
		close BARCODES;
		close OUT;
			
		system "rm Barcodes/AllSamplesRaw.txt";	
	
	
	
	
	}
	
	else {
		print "No directory named DemultiplexedFiles found.  If your data are already demultiplexed, put the fastq files in a folder named DemultiplexedFiles and try again.\n";
		exit;
	}	
	
}	





#Check to see if the data contain a restriction enzyme recognition sequence.  If so, get the expected sequence, and its length.

my $RELength;
my $RESeq;
my $RE = $RunArguments{re};

if ($RE !~ /^0/)  {
	$RESeq = $RE;
	$RELength = length($RESeq);
}

else {
	$RELength = 0;
}	


my $MinQualScore = $RunArguments{minQual};
my $MinNonZeroMean = $RunArguments{minDepth};
my $AllelicCriticalValue = $RunArguments{minIden};
my $P2AdaptorSeq = $RunArguments{P2};
my $MinForParalog = $RunArguments{minParalog};
my $QualityScoreType = $RunArguments{Phred};
my $numIndelsAllowed = $RunArguments{numIndels};
my $stringLength = $RunArguments{stringLength};
my $MaxProportionHeterozygotes = $RunArguments{MaxH};
$MaxProportionHeterozygotes = $MaxProportionHeterozygotes/100;
my $MinAlignmentLength;


print "\n\n";


print "Arguments entered are...\n";
for (keys %RunArguments) {
	print "$_\t$RunArguments{$_}\n";
}	


print "\n\nRunning AftrRAD...\n\n";

##################################################################################################################################################################################################
##################################################################################################################################################################################################

#Create directories for some of the AftrRAD output

mkdir "Output" unless (-d "Output");
mkdir "Output/RunInfo" unless (-d "Output/RunInfo");
mkdir "CandidateLoci" unless (-d "CandidateLoci");
mkdir "TempFiles" unless (-d "TempFiles");




#Go through the Barcodes directory to get the number and identity of data files to be included in the run.  

opendir BARCODEDIR, "$RunArguments{BarcodePath}";
my @FileNamesWExtension = grep { $_ ne '.' && $_ ne '..' && $_ ne '.DS_Store' } readdir(BARCODEDIR);
close BARCODEDIR;

my $NumOfRuns = @FileNamesWExtension;					
print "Recognized $NumOfRuns data files.\n\n";
print "File names are...\n @FileNamesWExtension \n\n";

my @BarcodeFileNames = ();

foreach my $filename (@FileNamesWExtension)  {
	$filename =~ s/\s//g;
	push (@BarcodeFileNames, $filename);
}	








#Create a master barcode file and hash that have all barcodes associated with their sample names.  Also, get the lengths of the longest and shortest barcodes.

my $MinBarcodeLength;
my $MaxBarcodeLength;
my %MasterBarcodes;
my @UnsortedBarcodeLengths = ();


if ($NumOfRuns == 1)  {
			
	system "cp $RunArguments{BarcodePath}$BarcodeFileNames[0] TempFiles/MasterBarcodeFile.txt";
	
	open MASTERBARCODEFILE, "TempFiles/MasterBarcodeFile.txt" or die$!;
	
	while(<MASTERBARCODEFILE>)  {
		
		if ($_ =~ /[A-Za-z]/) {
				chomp($_);
				my @TempArray = split(/\t/, $_);
				
				$MasterBarcodes{$TempArray[1]} = $TempArray[0];						my $CurrentBarcodeLength = length($TempArray[0]);
				push(@UnsortedBarcodeLengths, $CurrentBarcodeLength);
		}
	}
	close MASTERBARCODEFILE;


	my @SortedBarcodeLengths = sort {$a <=> $b} @UnsortedBarcodeLengths;			#Get min and max barcode lengths

	$MinBarcodeLength = $SortedBarcodeLengths[0];
	$MaxBarcodeLength = $SortedBarcodeLengths[-1];
	
	
}	


else  {
	
	
	foreach my $filename (@BarcodeFileNames)  {
		open TEMPFILE, "$RunArguments{BarcodePath}$filename" or die$!;
		while (<TEMPFILE>)  {
			if ($_ =~ /[A-Za-z]/) {
				chomp($_);
				my @TempArray = split(/\t/, $_);
				
				$MasterBarcodes{$TempArray[1]} = $TempArray[0];		#Note that there may be a rare scenario that this causes a problem in which the same individual is run with two different barcodes of different lengths, as the hash will only keep track of one of the barcodes.

			}
		}	
	}			
				
	open MASTERBARCODEFILE, ">TempFiles/MasterBarcodeFile.txt" or die$!;
	
	
	#Get min and max lengths of barcodes and shorten barcodes to shortest one.
		


	for (keys %MasterBarcodes)  {								#Concatenate all barcodes across all runs into MasterBarcode file.
		print MASTERBARCODEFILE "$MasterBarcodes{$_}\t$_\n";
		my $CurrentBarcodeLength = length($MasterBarcodes{$_});
		
		if ($CurrentBarcodeLength > 1)  {
				push (@UnsortedBarcodeLengths, $CurrentBarcodeLength);

		}
		
	}	
	
	my @SortedBarcodeLengths = sort {$a <=> $b} @UnsortedBarcodeLengths;			#Get min and max barcode lengths

		$MinBarcodeLength = $SortedBarcodeLengths[0];
		$MaxBarcodeLength = $SortedBarcodeLengths[-1];


	

	close MASTERBARCODEFILE;	
}		










#Create array that will store sample names for all individuals in the analysis (across all runs).
my @AllSampleNames = ();

#Store length of reads (all reads analyzed should have the same length)
my $ReadLength = ();





##################################################################################################################################################################################################
##################################################################################################################################################################################################
#Filter reads based on quality and demultiplex samples for each of the runs, one at a time.
##################################################################################################################################################################################################
##################################################################################################################################################################################################


foreach my $FileNumber (1..$NumOfRuns)  {
		
	
	
		#Determine what characters are "bad" for checking quality scores later on.

		my $Phred33String = '#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ';
		my $Phred64String = 'BCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh';
		my $BadCharacters;
					
		if (($QualityScoreType == 33) || ($QualityScoreType eq "Phred33") || ($QualityScoreType eq "phred33") || ($QualityScoreType eq "PHRED33"))  {
			$BadCharacters = substr($Phred33String,0,$MinQualScore);
		}
					
		elsif (($QualityScoreType == 64) || ($QualityScoreType eq "Phred64") || ($QualityScoreType eq "phred64") || ($QualityScoreType eq "PHRED64")) {
			$BadCharacters = substr($Phred64String,0,$MinQualScore);
		}
		
		else {
			print "Warning...did not recognize quality score type (i.e. Phred33 or Phred64).\n\n";
			next;
		}	

	
	
	
	
		#Open file that contains the barcodes and sample names for the current run and push the barcodes to array "CurrentFileBarcodes" and the sample names to array "CurrentFileSampleNames".
		
		
		open CURRENTBARCODEFILE, "$RunArguments{BarcodePath}$BarcodeFileNames[($FileNumber-1)]" or die$!;
		
		my @CurrentFileBarcodes = ();
		my @CurrentFileSampleNames = ();
		my %CurrentFileBarcodesAndSampleNames = ();
		my $BarcodeFileLineCounter = 0;

		while (<CURRENTBARCODEFILE>)  {
			
			$BarcodeFileLineCounter++;
			chomp ($_);
			my @TabSeparatedArray = split (/\t/, "$_");
			
			push (@CurrentFileBarcodes, $TabSeparatedArray[0]);
			push (@CurrentFileSampleNames, $TabSeparatedArray[1]);
			push (@AllSampleNames, $TabSeparatedArray[1]);
			
			$CurrentFileBarcodesAndSampleNames{$TabSeparatedArray[0]} = $TabSeparatedArray[1];
		}
		
		
		if ($BarcodeFileLineCounter <= 1)  {
			print "\nWarning...only one line recognized in barcode file $BarcodeFileNames[($FileNumber-1)].  Check format of this file if it contains more than one sample.\n";
		}
		
		
		#Create a hash that has all of the user barcodes for the current file as the keys.  We'll search these later and keep track of number of hits for each barcode.
		my %UserBarcodes;
			
			foreach my $Barcode (@CurrentFileBarcodes)  {
				
				$UserBarcodes{$Barcode} = 1;
			}	
		
		
		close CURRENTBARCODEFILE;
		
		
		
		print "The barcodes you entered for data file $FileNumber (file $FileNamesWExtension[$FileNumber-1]) are...\n @CurrentFileBarcodes \n\n";	#Print to check to make sure user barcodes are input correctly.
		print "The names you entered for data file $FileNumber (file $FileNamesWExtension[$FileNumber-1]) are...\n @CurrentFileSampleNames \n\n";
		
##################################################################################################################################################################################################
		#Open the fastq file for the current run. 
		
	
		open CURRENTDATAFILE, "$RunArguments{DataPath}$FileNamesWExtension[$FileNumber-1]" or die$!;
		open GOODSEQSGOODBARCODES, ">TempFiles/AllGoodSeqsWithAssignableBarcodes$FileNumber.txt" or die$!;
		open GOODSEQSNOBARCODES, ">>TempFiles/GoodSeqsNoBarcodes.txt" or die$!;
		
		
		my @FastqLines = ();
		my $FastqCounter = 0;
		my $TotalSequenceNumberCounter = 0;
		my $ReadBarcode;
		my $ReadRecogSite;
		my $TotalBarcodeNonMatches;
		my $NumTrueBarcodeMatches = 0;
		my @TempNearMatchesForUnidentifiedBarcode = ();
		my @SplitUnidentifiedBarcode = ();
		my %ReplacedBarcodes = ();
		my $NumberOfReplacedBarcodes = 0;
		my $BarcodesMoreThanTwoBasesAway = 0;
		my $BarcodesCloseToTwoOrMore = 0;
		my @EditedSeq = ();
		my $RawReadBase = 0;
		my $RetainedSequenceNumberCounter = 0;
		my $LowQualitySeqs = 0;
		my %BarcodeReplacementCounts = ();
		my $SequencesWithHomopolymersRemoved = 0;
		my $NumReadsWithP2AdaptorSeq = 0;
		my $NumReadsDiscardedForBadRecognitionSite = 0;
		my %UniqueBarcodesWithCounts;
		my %BarcodesTwoAway;
		my %BarcodesOneAwayFromTwo;
		my $ReadsWithBadBarcodesAndRESites = 0;

		
		
		#Join barcodes to search for at the beginning of the reads
		my @JoinedBarcodes = join("|",@CurrentFileBarcodes);
		
		#Create array that stores barcodes that have all been trimmed to the same length (in case barcodes have different lengths), along with their associated sample names.
		my @ShortenedBarcodes = ();

		for (keys %CurrentFileBarcodesAndSampleNames)  {
			my $CurrentShortenedBarcode = substr ($_, 0, $MinBarcodeLength);
			push (@ShortenedBarcodes, $CurrentShortenedBarcode);
			push (@ShortenedBarcodes, $CurrentFileBarcodesAndSampleNames{$_});
		}

		my %CurrentFileShortenedBarcodesAndSampleNames = @ShortenedBarcodes;	#Store the shortened barcodes in a hash with their associated sample names.
	



		

		
		
		

		#Begin running through fastq file

		while (<CURRENTDATAFILE>)  { 		
				chomp ($_);
				push (@FastqLines, $_);				#Get the 4 lines for the current read in an array
				$FastqCounter++;
  
				if ($FastqCounter == 4)  {			#Have the 4 lines in the array - now evaluate this current read
     
					$TotalSequenceNumberCounter++;
				
					if ($TotalSequenceNumberCounter =~ /000000$/)  {
						print "Filtering sequence $TotalSequenceNumberCounter.\n";
					}	
				
					
					my $CurrentSeq = $FastqLines[1];
					my $CurrentQuality = $FastqLines[3];
				
##################################################################################################################################################
				
				
					#If the read has the P2 adaptor, throw it out and go to next read.
					if ($CurrentSeq =~ /($P2AdaptorSeq)/i)  {
						$NumReadsWithP2AdaptorSeq++;
						@FastqLines = ();
						$FastqCounter = 0;
						next;
					}
				
				
##################################################################################################################################################


					#Check to see if the read contains an exact match to a barcode
			
					my @MatchedBarcodes = ($CurrentSeq =~ /^(@JoinedBarcodes)/);
			
					if (@MatchedBarcodes)  {						#The read does contain an exact barcode match
						my $CurrentMatchedBarcode = $MatchedBarcodes[0];
					
						my $CurrentValue = $UserBarcodes{$CurrentMatchedBarcode};      #Gets the number of times we've seen this barcode prior to this observation.
						$CurrentValue++;
						$UserBarcodes{$CurrentMatchedBarcode} = $CurrentValue;		#Keep track of the number of matches for each barcode
						$NumTrueBarcodeMatches++;					#Keep track of the total number of matches for all barcodes
									
				
						my $CurrentBarcodeLength = length($CurrentMatchedBarcode);
				
						#Now check restriction enzyme recognition sequence (if one exists in the data).  If it is 2 or more bases from what it's supposed to be, throw it out.
				
						if ($RELength != 0)  {
						
							$ReadRecogSite = substr($CurrentSeq, ($CurrentBarcodeLength), $RELength);
				
							if ($ReadRecogSite !~ m/$RESeq/i)  {
								my @SplitReadRecogSite = split (//, $ReadRecogSite);
								my @RESeq = split (//, $RESeq);
					
								my $CurrentElementNumberInComparison = 0;
								my $NumMatches = 0;
          
								foreach my $CurrentElement (@SplitReadRecogSite)  {
									$CurrentElement = uc($CurrentElement);
									if ((uc($RESeq[$CurrentElementNumberInComparison])) eq "$CurrentElement")  {
										$NumMatches++;
									}  
									$CurrentElementNumberInComparison++;
								}  
       
								if ($NumMatches < $RELength-1)  {	#There is more than one mismatch - go to next read in fastq file
									$NumReadsDiscardedForBadRecognitionSite++;
									@FastqLines = ();
									$FastqCounter = 0;
									next;
								}    
							}
						
						}
				
				
						#At this point, we have an exact barcode match, and we have an intact RE site (or at least within 1)
				
						#Determine if the read needs trimmed to account for barcodes of different sizes, and trim if necessary
				
						if ($CurrentBarcodeLength < $MaxBarcodeLength)  {
						
							my $Difference = $CurrentBarcodeLength-$MaxBarcodeLength;
							$CurrentSeq = substr($CurrentSeq,0,$Difference);
							$CurrentQuality = substr($CurrentQuality,0,$Difference);
						
						}
					
						#Check Phred scores and for homopolymers
					
					
						my $BeginningReadPosition = $CurrentBarcodeLength+$RELength;
						my $SeqAfterBarcodeAndRESite = substr($CurrentSeq, $BeginningReadPosition);
						my $QualityAfterBarcodeAndRESite = substr($CurrentQuality, $BeginningReadPosition);
					
						if ($QualityAfterBarcodeAndRESite =~ /[$BadCharacters]/)  {  #The read contains at least one low quality base - we will eliminate the read.
							$LowQualitySeqs++;
							$FastqCounter = 0;
							@FastqLines = ();
							next;
						}
					
						elsif  (($SeqAfterBarcodeAndRESite =~ /A{$stringLength}/)|($SeqAfterBarcodeAndRESite =~ /T{$stringLength}/)|($SeqAfterBarcodeAndRESite =~ /G{$stringLength}/)|($SeqAfterBarcodeAndRESite =~ /C{$stringLength}/)) {
						
							$SequencesWithHomopolymersRemoved++;
							$FastqCounter = 0;
							@FastqLines = ();
							next;
						}
			
						else {	
							$RetainedSequenceNumberCounter++;		#This keeps track of the total number of sequences that will be included in the final dataset.
							my $BarcodeToPrint = substr($CurrentMatchedBarcode,0,$MinBarcodeLength);
							print GOODSEQSGOODBARCODES "$BarcodeToPrint$SeqAfterBarcodeAndRESite\n";
							print GOODSEQSNOBARCODES "$SeqAfterBarcodeAndRESite\n";
						}	
					}
			
			
			
			
			
			
			
			
			
			
			
			
					else {		#The read does not contain an exact barcode match
			
						

						
						$TotalBarcodeNonMatches++;
						my $ReadBarcodeLength;
						
						
						#If barcodes are different sizes, find the restriction site first, then take the bases before this as the barcode
				
						
						if ($MaxBarcodeLength != $MinBarcodeLength)  {
						
						
							my $MaxPositionForRESite = $MaxBarcodeLength+$RELength;
				
							my $SubstringForRESiteSearch = substr($CurrentSeq,0,$MaxPositionForRESite);
				
				
							if ($SubstringForRESiteSearch !~ m/$RESeq/i)  { 	#this would mean we have a mismatch at both the barcode and the restriction site.
								$ReadsWithBadBarcodesAndRESites++;
								$FastqCounter = 0;
								@FastqLines = ();
								next;
							}
				
							
						
							#Alteratively, we have a mismatched barcode, but an intact RE site.
					
							$SubstringForRESiteSearch =~ m/$RESeq/gi;
					
							my $EndLocationOfRESite = pos($SubstringForRESiteSearch);
					
							my $LengthOfBarcode = $EndLocationOfRESite-$RELength;
							my $BarcodeEndPosition = $LengthOfBarcode;
					
				
							#Get current barcode (contains at least one error)
					
							$ReadBarcode = substr($CurrentSeq,0,$BarcodeEndPosition);
							$ReadBarcodeLength = length($ReadBarcode);
				       			
				       		}	
				       			
								
				       		else {			#If all barcodes are the same size
				       			
				       			$ReadBarcode = substr($CurrentSeq,0,$MinBarcodeLength);
				       			$ReadBarcodeLength = $MinBarcodeLength;
				       		}	
				       			
				       			
				       		
				       		
				       		
						#Check the barcode to see if it can be assigned unambiguously to one of the expected barcode sequences.
					
					
						@TempNearMatchesForUnidentifiedBarcode = ();
						my $NumberOfMatchesBetweenUnidentifiedAndUserBarcode = 0;
						@SplitUnidentifiedBarcode = split (//, $ReadBarcode);
              
						foreach my $truebarcode (@CurrentFileBarcodes)  {
				       		my $CurrentFileUserBarcodeLength = length($truebarcode);
				        		
				     			if ($ReadBarcodeLength == $CurrentFileUserBarcodeLength)  {
				       				my @SplitUserBarcode = split (//, $truebarcode);
				       				my $SplitUnidentifiedBarcodeElementCounter = 0;
				       				$NumberOfMatchesBetweenUnidentifiedAndUserBarcode = 0;
                    
				       				foreach my $d (@SplitUserBarcode)  {
				       					if ($d eq $SplitUnidentifiedBarcode[$SplitUnidentifiedBarcodeElementCounter])  {
				        					$NumberOfMatchesBetweenUnidentifiedAndUserBarcode++;
				        				}
				        				$SplitUnidentifiedBarcodeElementCounter++;
				        			}
       			    
       			    
				        			if ($NumberOfMatchesBetweenUnidentifiedAndUserBarcode >= $ReadBarcodeLength-1)  {
				        				push (@TempNearMatchesForUnidentifiedBarcode, "$truebarcode");
				        			} 
				        			
				        		}	
				        		
				        		} #Closes foreach in CurrentFileBarcodes (after the unidentified barcode has been checked against each of the user barcodes, and the number of differences recorded.)
       		     
				        		my $NumberOfUserBarcodesWithNearMatch = @TempNearMatchesForUnidentifiedBarcode;
       		   
       		   
###################################################################################################################################################################################################################       		   
       		   
       		  		 	


							if ($NumberOfUserBarcodesWithNearMatch == 1)  { #If there's only one user barcode that is one base away from the unidentified barcode, the unidentified barcode gets edited.
								my $TempReadBarcode = $ReadBarcode;
								my $MatchToHashKey = 0;
								
								for (keys %BarcodeReplacementCounts)  {		
									
									if ($_ eq $TempReadBarcode)  {
										
										$MatchToHashKey = 1;
										
										my $CurrentValue = "$BarcodeReplacementCounts{$_}";
										$CurrentValue++;
										$BarcodeReplacementCounts{$_} = $CurrentValue;
										last;
									}
								}
							
								if ($MatchToHashKey == 0) {
									$BarcodeReplacementCounts{$TempReadBarcode} = 1;		#This has the barcode that actually was in the sequence, and the number of times that error barcode was read.
								}    
       		          
       		               
								my $RepairedBarcode = $TempNearMatchesForUnidentifiedBarcode[0];
								$ReplacedBarcodes{$TempReadBarcode} = "$RepairedBarcode";
								$NumberOfReplacedBarcodes++;
								
						
								
								my $CurrentValue = $UserBarcodes{$RepairedBarcode};
								$CurrentValue++;
								$UserBarcodes{$RepairedBarcode} = $CurrentValue;		#Keep track of the number of matches for each barcode.
									
								
								
								
								
								
								
       		  		 	   
								#At this point, we have a barcode that we were able to edit and assign.
				
								#Determine if the read needs trimmed to account for barcodes of different sizes, and trim if necessary
								
								my $CurrentBarcodeLength = length($RepairedBarcode);
							
								if ($CurrentBarcodeLength < $MaxBarcodeLength)  {
						
									my $Difference = $CurrentBarcodeLength-$MaxBarcodeLength;
									$CurrentSeq = substr($CurrentSeq,0,$Difference);
									$CurrentQuality = substr($CurrentQuality,0,$Difference);
						
								}
					
								#Check Phred scores and for homopolymers
					
					
								my $BeginningReadPosition = $CurrentBarcodeLength+$RELength;
								my $SeqAfterBarcodeAndRESite = substr($CurrentSeq, $BeginningReadPosition);
								my $QualityAfterBarcodeAndRESite = substr($CurrentQuality, $BeginningReadPosition);
					
								if ($QualityAfterBarcodeAndRESite =~ /[$BadCharacters]/)  {  #The read contains at least one low quality base
									$LowQualitySeqs++;
									$FastqCounter = 0;
									@FastqLines = ();
									next;
								}
					
								elsif  (($SeqAfterBarcodeAndRESite =~ /A{$stringLength}/)|($SeqAfterBarcodeAndRESite =~ /T{$stringLength}/)|($SeqAfterBarcodeAndRESite =~ /G{$stringLength}/)|($SeqAfterBarcodeAndRESite =~ /C{$stringLength}/)) {
						
									$SequencesWithHomopolymersRemoved++;
									$FastqCounter = 0;
									@FastqLines = ();
									next;
								}
			
								else {	
									$RetainedSequenceNumberCounter++;		#This keeps track of the total number of sequences that will be included in the final dataset.
									my $BarcodeToPrint = substr($RepairedBarcode,0,$MinBarcodeLength);
									print GOODSEQSGOODBARCODES "$BarcodeToPrint$SeqAfterBarcodeAndRESite\n";
									print GOODSEQSNOBARCODES "$SeqAfterBarcodeAndRESite\n";
								
								}
							
							
						
							}	
							
							
							
		
       		  		 
							elsif ($NumberOfUserBarcodesWithNearMatch == 0)  {  #If the barcode is at least two bases from each of the user barcodes, it is added to hash BarcodesTwoAway.
       		  		 	
								if (exists($BarcodesTwoAway{$ReadBarcode}))  {
									$BarcodesTwoAway{$ReadBarcode}++;
								}
       		  		 	 
								else {
									$BarcodesTwoAway{$ReadBarcode} = 1;
								}	 
       		  		 	 	 
							} #close elsif   
       		      
       		      
       		  		 	
							elsif ($NumberOfUserBarcodesWithNearMatch >= 2)  {  #If the unidentified barcode is one base away from more than one user barcode, it is added to an array if it's not already there.
							
								if (exists($BarcodesOneAwayFromTwo{$ReadBarcode}))  {
									$BarcodesOneAwayFromTwo{$ReadBarcode}++;
								}
       		  		 	 
								else {
									$BarcodesOneAwayFromTwo{$ReadBarcode} = 1;
								}	 
       		  		 	 
       	
							}
				
       		  		 
					}
       		  		 
         
					$FastqCounter = 0;
					@FastqLines = ();
	
				}  #Finished with the current sequence (set of 4 fastq lines)


		} #Finished with the current data file.
  
		   
		close CURRENTDATAFILE;
		close GOODSEQSGOODBARCODES;
		close GOODSEQSNOBARCODES;

	
		
#####################################################################################################################################################################################
		#Print a Report file for each Illumina run. 
#############################################################################################################################################################################
		
		open REPORT, ">Output/RunInfo/Report_$FileNamesWExtension[$FileNumber-1].txt"; 
		
		print REPORT "Run information for datafile $FileNamesWExtension[$FileNumber-1]\n\n";
		
		print REPORT "Parameters used for this run...\n\n";
		
		for (keys %RunArguments) {
			print REPORT "$_\t$RunArguments{$_}\n";
		}
		

		print REPORT "The total number of reads in the fastq file is $TotalSequenceNumberCounter.\n\n"; 

		my $RemainingNumberOfSeqs1a = "$TotalSequenceNumberCounter"-"$NumReadsWithP2AdaptorSeq";
		print REPORT "You entered $P2AdaptorSeq as the beginning of the P2 adaptor.  $NumReadsWithP2AdaptorSeq reads contained this string and were discarded, leaving $RemainingNumberOfSeqs1a reads for further analysis.\n\n";

		
		print REPORT "Of these $RemainingNumberOfSeqs1a sequences, $NumTrueBarcodeMatches reads were an exact match to one of your barcodes and $TotalBarcodeNonMatches were not.\n\n";
		print REPORT "Of the $TotalBarcodeNonMatches not matching a barcode exactly, $NumberOfReplacedBarcodes reads had barcodes that could be confidently assigned to one of your barcodes.\n\n";
		
		my $TotalSeqsIDd = $NumTrueBarcodeMatches+$NumberOfReplacedBarcodes;
		print REPORT "This left $TotalSeqsIDd sequences.\n\n";
		
		my $RemainingNumberOfSeqs1b = $TotalSeqsIDd - $NumReadsDiscardedForBadRecognitionSite;
		print REPORT "Of these $TotalSeqsIDd reads, $NumReadsDiscardedForBadRecognitionSite contained restriction enzymes recognition sites that were 2 or more bases away from your recognition site and were discarded, leaving $RemainingNumberOfSeqs1b reads for further analysis.\n\n";

		
		my $RemainingNumberOfSeqs3 = "$RemainingNumberOfSeqs1b"-"$LowQualitySeqs";
		print REPORT "Of these $RemainingNumberOfSeqs1b sequences, $LowQualitySeqs were removed because they contained at least one base that had a quality score below the minimum value of $MinQualScore, leaving $RemainingNumberOfSeqs3 sequences.\n\n";

		
		print REPORT "Of these $RemainingNumberOfSeqs3, $SequencesWithHomopolymersRemoved additional reads were removed from the dataset because they contained a homopolymer at least $stringLength bases long.\n\n";
		print REPORT "This left $RetainedSequenceNumberCounter sequences for further analysis.\n";

		
		
		
		
#####################################################################################################################################################################################
		#Print a Barcode file for each Illumina run. 
#############################################################################################################################################################################
		
		open BARCODES, ">Output/RunInfo/BarcodeInfo_$FileNamesWExtension[$FileNumber-1].txt";
		open BARCODECOUNTSTOPLOT, ">TempFiles/BarcodeCountsToPlot.txt" or die$!;
		open BARCODENAMESTOPLOT, ">TempFiles/BarcodeNamesToPlot.txt" or die$!;
		
		print BARCODES "Total number of matches to each barcode...\n";

		foreach my $CurrentKey (keys %UserBarcodes)  {
       		  	print BARCODES "$CurrentKey\t$CurrentFileBarcodesAndSampleNames{$CurrentKey}\t$UserBarcodes{$CurrentKey}\n";
       		  	print BARCODENAMESTOPLOT "$CurrentKey\t";
       		  	print BARCODECOUNTSTOPLOT "$UserBarcodes{$CurrentKey}\t";
			
       		}
		
		print BARCODES "\n";

		
		
		
		#Check to see if any barcodes that are more than 2 bases away from any of the expected barcodes are represented in large numbers in the dataset.
		
		my @UnexpectedBarcodes = ();
		
		for (keys %BarcodesTwoAway)  {		
									
			if ($BarcodesTwoAway{$_} > 5000)   {
				
				push (@UnexpectedBarcodes, $_);
			}
		}
		
		if (@UnexpectedBarcodes)  {
			
			print BARCODES "The following barcodes were read at least 5000 times, and were not included in the barcodes file.  Check to make sure they were not accidentally left out.\n";
			
			foreach my $unexpectedbarcode (@UnexpectedBarcodes)  {
				print BARCODES "$unexpectedbarcode\t$BarcodesTwoAway{$unexpectedbarcode}\n";
			}	
			
			#print BARCODES "@UnexpectedBarcodes\n\n";
		}
		
		else {
			print BARCODES "No unexpected barcodes identified in >5000 reads.\n\n";
		}		
		
	
	close BARCODES;
	close BARCODECOUNTSTOPLOT;
	close BARCODENAMESTOPLOT;
#############################################################################################################################################################################		
	#Plot barcode read counts.  #First, update the R script for the current file number.

	# open RSCRIPT, "RScripts/PlotBarcodeReadCounts.R" or die$!;
	# open RSCRIPTOUT, ">RScripts/PlotBarcodeReadCountsUpdate.R" or die$!;
# 
	# while(<RSCRIPT>)  {
		# chomp;
		# if ($_ =~ /^pdf/)  {
			# print RSCRIPTOUT "pdf(\"Output/RunInfo/BarcodeReadCounts_$FileNamesWExtension[$FileNumber-1].pdf\")\n";
		# }
	# 
		# else {
			# print RSCRIPTOUT "$_\n";
		# }
# 
	# }
# 
	# close RSCRIPT;
	# close RSCRIPTOUT;
	# close BARCODENAMESTOPLOT;
	# close BARCODECOUNTSTOPLOT;
# 
# 
	# system "R --vanilla --slave < RScripts/PlotBarcodeReadCountsUpdate.R";

############################################################################################################################################################################
		
	#Create a unique file in TempFiles directory for each individual in the dataset.  Individual names are used to name these files instead of barcodes, so that barcodes can be re-used across different datasets.
	#Go through sorted AllGoodSeqsWithAssignableBarcodes, and parse these seqs without the barcode into their appropriate Individual file. 
		
		print "\nDemultiplexing samples for data file $FileNamesWExtension[$FileNumber-1].\n";
		
		open UNSORTED, "TempFiles/AllGoodSeqsWithAssignableBarcodes$FileNumber.txt" or die$!;

		my @ToSort = ();

		while (<UNSORTED>)  {
			chomp($_);
			if ($_ =~ /[a-zA-Z]/)  {
				push (@ToSort, $_);
			}
		}

		my @Sorted = sort { lc($a) cmp lc($b) } @ToSort;

		close UNSORTED;

		open SORTED, ">TempFiles/GoodSeqsWithBarcodesSorted.txt" or die$!;

		foreach my $read (@Sorted)  {
			print SORTED "$read\n";
		}	

		close SORTED;
		
	
		open GOODSEQSFILE, "TempFiles/GoodSeqsWithBarcodesSorted.txt" or die$!;
		
		my $CurrentGoodBarcode = "RandomStartingBarcodeSeq";
		
		open CURRENTWRITEFILE, ">TempFiles/IndividualRandomStartingBarcodeSeq.txt" or die$!;
		
		while (<GOODSEQSFILE>)  {
			
			my $TempBarcode = substr $_,0, $MinBarcodeLength;
			
			if ($TempBarcode !~ /$CurrentGoodBarcode/)  {
				
				close CURRENTWRITEFILE;
				
				chomp($_);
				my $NewIndividualName = $CurrentFileShortenedBarcodesAndSampleNames{$TempBarcode};
				
				open CURRENTWRITEFILE, ">>TempFiles/Individual$NewIndividualName.txt" or die$!;
				
				my $TempStringToPrint = substr $_, $MinBarcodeLength;
				
				print CURRENTWRITEFILE "$TempStringToPrint\n";
				$CurrentGoodBarcode = $TempBarcode;
				$MinAlignmentLength = length($TempStringToPrint)-$numIndelsAllowed;
			}
			
			else {
				chomp($_);
				my $TempStringToPrint = substr $_, $MinBarcodeLength;
				
				print CURRENTWRITEFILE "$TempStringToPrint\n";
				$MinAlignmentLength = length($TempStringToPrint)-$numIndelsAllowed;
			}
		}
		
		system "rm TempFiles/IndividualRandomStartingBarcodeSeq.txt";
		close CURRENTWRITEFILE;	
		close REPORT; 
			
############################################################################################################################################################################
	#Get unique seqs for each individual and counts for each.

	open UNIQUESFORPLOT, ">TempFiles/UniqueCountsForPlot.txt" or die$!;
	close UNIQUESFORPLOT;

	print "\nIdentifying unique sequences for each individual.\n";
	
	open SAMPLENAMESFORPLOT, ">TempFiles/SampleNamesForPlot.txt" or die$!;

		foreach my $g (@CurrentFileSampleNames)  {
		
			
			print SAMPLENAMESFORPLOT "$g\t";
		
			my @TempArrayOfUniqueSeqsPerIndividual = ();
		
			
			if (-e "TempFiles/Individual$g.txt")  {			#Check to make sure the file exists (the individual had at least one unique read that passed all of the filtering).  Without this check, the program will crash if there is an individual with no reads (Individual$g file doesn't exist).
			
			
				open UNSORTED, "TempFiles/Individual$g.txt" or die$!;

				my @ToSort = ();

				while (<UNSORTED>)  {
					chomp($_);
					if ($_ =~ /[a-zA-Z]/)  {
						push (@ToSort, $_);
					}
				}

				my @Sorted = sort { lc($a) cmp lc($b) } @ToSort;

				close UNSORTED;

				open SORTED, ">TempFiles/SortedIndividual$g.txt" or die$!;

				foreach my $read (@Sorted)  {
					print SORTED "$read\n";
				}	

				close SORTED;
			
		
			
				system "uniq -c TempFiles/SortedIndividual$g.txt TempFiles/UniqueWithCountsIndividual$g.txt";
		
				open REPORT, ">>Output/RunInfo/Report_$FileNamesWExtension[$FileNumber-1].txt";  
		
				print REPORT "The number of unique sequences for individual $g is";
		
				close REPORT;
		
				system "wc -l TempFiles/UniqueWithCountsIndividual$g.txt | tee -a Output/RunInfo/Report_$FileNamesWExtension[$FileNumber-1].txt TempFiles/UniqueCountsForPlot.txt";
				open UNIQUESFORPLOT, ">>TempFiles/UniqueCountsForPlot.txt" or die$!;
				print UNIQUESFORPLOT "\t";
				close UNIQUESFORPLOT;
				
			}	
			
			
			else {
				open UNIQUESFORPLOT, ">>TempFiles/UniqueCountsForPlot.txt" or die$!;
				print UNIQUESFORPLOT "0\t";
				close UNIQUESFORPLOT;
			}	
			
		}

		close SAMPLENAMESFORPLOT;     
		
	system "rm TempFiles/AllGoodSeqsWithAssignableBarcodes$FileNumber.txt";

		
}


##################################################################################################################################################################################################


#Plot unique numbers of reads in each individual across all individuals analyzed (across all dataset, if more than one).

# open UNIQUESFORPLOT, "TempFiles/UniqueCountsForPlot.txt" or die$!;
# open UNIQUESFORPLOTEDIT, ">TempFiles/UniqueCountsForPlotEdit.txt" or die$!;
# 
# while (<UNIQUESFORPLOT>)  {
	# $_ =~ s/^\t+//;
	# $_ =~ s/^\s+//;
	# print UNIQUESFORPLOTEDIT "$_";
# }	
# 
# close UNIQUESFORPLOTEDIT;
# close UNIQUESFORPLOT;
# 
# open UNIQUESFORPLOTEDIT, "TempFiles/UniqueCountsForPlotEdit.txt" or die$!;
# open UNIQUESFORPLOTUPDATE, ">TempFiles/UniqueCountsForPlotUpdate.txt" or die$!;
# 
# while (<UNIQUESFORPLOTEDIT>)  {
	# if ($_ =~ /[0-9]/) {
	   # chomp($_);
	   # my @TempArray = split(/\s+/, $_);
	   # print UNIQUESFORPLOTUPDATE "$TempArray[0]\t";
	# }   
# }
# 
# close UNIQUESFORPLOT;
# close UNIQUESFORPLOTUPDATE;
	# 
# 
# system "R --vanilla --slave < RScripts/PlotUniqueCounts.R";








##################################################################################################################################################################################################
##################################################################################################################################################################################################
#Finished going through each individual data file (if more than one) to filter reads based on quality scores, demultiplex reads, etc.  
#Now, use all of the reads (across all of the datafiles, if more than one) to identify loci, and map alleles to these loci.
##################################################################################################################################################################################################
##################################################################################################################################################################################################



#Find all unique reads across the entire dataset as first step in locus identification.

print "\nIdentifying all unique sequences within the dataset.\n\n";

system "sort TempFiles/GoodSeqsNoBarcodes.txt -o TempFiles/AllSeqsSorted.txt";	
system "uniq TempFiles/AllSeqsSorted.txt TempFiles/AllUniquesSorted.txt";

system "rm TempFiles/GoodSeqsNoBarcodes.txt";
system "rm TempFiles/AllSeqsSorted.txt";

#Get length of each read
open SORTED, "TempFiles/AllUniquesSorted.txt";

while (<SORTED>)  {
	chomp($_);
	$ReadLength = length($_);
	last;
}






##################################################################################################################################################################################################

#Try to identify error reads that weren't identified using Phred scores.  These are likely reads that have had errors added in during PCR amplification, so Phred scores are expected to be high at these sites.
#We'll try to identify these based on the mean nonzero read count across all individuals in the dataset.


#Create a file called AllUniquesForErrorTest.txt that is just a copy of AllUniquesSorted.txt.  Then append the counts for each individual to this.
#Go through each individual file and create a hash that has the unique seqs and the counts.  Then append these to the AllUniquesForErrorTest file.  If the seq is missing, append "0".
#The completed file is called AllReadsAndDepths.txt.


mkdir "TempFiles/ErrorReadTest" unless(-d "TempFiles/ErrorReadTest");

system "cp TempFiles/AllUniquesSorted.txt TempFiles/ErrorReadTest/AllUniquesForErrorTest.txt";
print "Creating file to test mean read counts.\n";

my $IndividualCounter = 0;
my $Number = 0;

foreach my $Individual (@AllSampleNames) {
	
	if (-e "TempFiles/UniqueWithCountsIndividual$Individual.txt")  {
	
	  $Number++;
	  my $NumberPlus1 = $Number+1;
	
	  if ($IndividualCounter == 0) {			#First individual
		
		$IndividualCounter++;
		open CURRENTINDIVIDUAL, "TempFiles/UniqueWithCountsIndividual$Individual.txt" or die$!;
		
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
	
		open VARIANCEDATA, "TempFiles/ErrorReadTest/AllUniquesForErrorTest.txt" or die$!;
		open VARIANCEUPDATE, ">TempFiles/ErrorReadTest/ErrorUpdate$Number.txt" or die$!;
		
		while (<VARIANCEDATA>) {
		   chomp($_);
		   my @TempArray = split (/\t/, $_);
		   my $CurrentSeq = $TempArray[0];
		   my $Value = 0;
		
		   if (exists($TempHash{$CurrentSeq})) {
			my $CurrentCount = $TempHash{$CurrentSeq};
			print VARIANCEUPDATE "@TempArray\t$CurrentCount\n";
		   }
		
		   else {
			print VARIANCEUPDATE "@TempArray\t0\n";
		   }	
		}
	
		close VARIANCEDATA;
		close VARIANCEUPDATE;
		close CURRENTINDIVIDUAL;
	
	  }

	

	  else {				#Not the first individual
		
	
		open CURRENTINDIVIDUAL, "TempFiles/UniqueWithCountsIndividual$Individual.txt" or die$!;
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
		
		open VARIANCEUPDATE, "TempFiles/ErrorReadTest/ErrorUpdate$NumberMinus1.txt" or die$!;
		open VARIANCEWRITE, ">TempFiles/ErrorReadTest/ErrorUpdate$Number.txt" or die$!;
		
		while (<VARIANCEUPDATE>) {
		   chomp($_);
		   $_ =~ s/\s/\t/;
		   my @TempArray = split (/\t/, $_);
		   my $CurrentSeq = $TempArray[0];
		   my $Value = 0;
		
		   if (exists($TempHash{$CurrentSeq})) {
			my $CurrentCount = $TempHash{$CurrentSeq};
			print VARIANCEWRITE "@TempArray\t$CurrentCount\n";
		
		   }
		
		   else {
			print VARIANCEWRITE "@TempArray\t0\n";
			
		   }	
		}
	
		close VARIANCEUPDATE;
		close VARIANCEWRITE;
		close CURRENTINDIVIDUAL;
		system "rm TempFiles/ErrorReadTest/ErrorUpdate$NumberMinus1.txt"
	
	  }
	
	}	
}

	open VARIANCEWRITE, "TempFiles/ErrorReadTest/ErrorUpdate$Number.txt" or die$!;
	open VARIANCEFINAL, ">TempFiles/ErrorReadTest/AllReadsAndDepths.txt" or die$!;
	
	while (<VARIANCEWRITE>)  {
		chomp($_);
		$_ =~ s/\s/\t/g;
		print VARIANCEFINAL "$_\n";
	}

close VARIANCEWRITE;
close VARIANCEFINAL;

system "rm TempFiles/ErrorReadTest/AllUniquesForErrorTest.txt";
system "rm TempFiles/ErrorReadTest/ErrorU*.txt";



#Go through file AllReadsAndDepths and obtain the mean read count (not counting zeros) for each read.  This value is compared to the set threshold to determine whether the read is discarded as error.
#File ErrorTestOut.txt contains the reads that exceed the threshold, and presumably are real (non-error) reads.  They will be used in the next step to begin identifying loci.



open INFILE, "TempFiles/ErrorReadTest/AllReadsAndDepths.txt" or die$!;
open VARIANCEOUT, ">TempFiles/ErrorReadTest/ErrorTestOut.txt" or die$!;
open VARIANCEERRORS, ">TempFiles/ErrorReadTest/VarianceErrors.txt" or die$!;
open SEQSWITHZEROCOUNTS, ">TempFiles/ErrorReadTest/SeqsWithZeroCounts.txt" or die$!;

print "Testing nonzero means in read counts to identify error reads.\n\n";

while(<INFILE>) {
	
	my $Sum = 0;
	my $SampleSize = 0;
	my $Mean = 0;
		
	my @TempArray = split(/\t/, "$_");
	my $ArraySize = @TempArray;
	
		
	foreach my $i (@TempArray[1..$ArraySize-1])  {
			
		if ($i != 0) {
		  $Sum = $Sum+$i;
		  $SampleSize++;
		}
	}
	
		
	if ($SampleSize != 0) {
			
		   $Mean = $Sum/$SampleSize;
		  
	}	
		
	else {
		print SEQSWITHZEROCOUNTS "$TempArray[0]\n";
		next;
	}	
		
		
		if ($Mean >= $MinNonZeroMean)  {
			print VARIANCEOUT "$TempArray[0]\n";
		}
		
		else {
			print VARIANCEERRORS "$TempArray[0]\n";
		}	
		
}




open VARIANCEOUT, "TempFiles/ErrorReadTest/ErrorTestOut.txt" or die$!;
my @UniqueNonErrorSeqsInTotalDataset;

while (<VARIANCEOUT>)  {
	chomp($_);
	push (@UniqueNonErrorSeqsInTotalDataset, $_);
}	

my $NumOfUniqueNonErrorSeqsInTotalDataset = @UniqueNonErrorSeqsInTotalDataset;



my @UniqueNamesNonError = ();
my $NumOfUniqueNonErrorSeqsCounter = 0;


foreach (0..$NumOfUniqueNonErrorSeqsInTotalDataset-1)  {
  $NumOfUniqueNonErrorSeqsCounter++;
  push (@UniqueNamesNonError, ">Seq_$NumOfUniqueNonErrorSeqsCounter");
}  



######################################################################################################################################################################
######################################################################################################################################################################

#Perform heuristic search to find potentially allelic pairs that will subsequently be aligned.
#The goal here is to reduce the number of pairwise alignments that need to be performed by eliminating pairs of reads that are very unlikely to be allelic.
#The lengths of the subsequences, and the number that must match, depend on the total sequence length.  This likely could still be optimized from where it's currently at.

print "Performing heuristic search of all pairwise comparisons to identify potentially allelic pairs.\n";
my $ElementNumber = 0;
my @PotentialLoci = ();


if ($ReadLength < 60)  {		#The length of the matching subsets we're looking for depends on the seq length.


 foreach (0..$NumOfUniqueNonErrorSeqsInTotalDataset-1)  {
  
  if ($_ =~ /000$/)  {
  	  
  	  print "Rapidly searching unique sequence number $_ of $NumOfUniqueNonErrorSeqsCounter for potentially allelic read pairs.\n";
  }
  
  push (@PotentialLoci, "Potential_Locus");
  push (@PotentialLoci, ">$UniqueNamesNonError[$ElementNumber]");
  push (@PotentialLoci, "$UniqueNonErrorSeqsInTotalDataset[$ElementNumber]");
  
  
  
  	my @SplitSeq = split (//, "$UniqueNonErrorSeqsInTotalDataset[$ElementNumber]");
  	my @Beginning7 = @SplitSeq[0..5];
  	my $Beginning7 = join "", @Beginning7;
  	my @Middle7A = @SplitSeq[10..15];
  	my $Middle7A = join "", @Middle7A;
  	my @Middle7B = @SplitSeq[17..22];
  	my $Middle7B = join "", @Middle7B;
  	my @Middle7C = @SplitSeq[26..31];
  	my $Middle7C = join "", @Middle7C;
  	my @Ending7 = @SplitSeq[-6..-1];
  	my $Ending7 = join "", @Ending7;

     	my $CountCycle = 0;					
     
     	my $ElementNumberB;	
     
     	my $PopTester = 0;
     
     	foreach my $p (@UniqueNonErrorSeqsInTotalDataset[($ElementNumber+1)..($NumOfUniqueNonErrorSeqsInTotalDataset-1)])  {
         
         	$CountCycle++;					
         
         	$ElementNumberB = $ElementNumber+$CountCycle;
         
         	
         	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7A/) && ($p =~ /$Middle7B/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
           	
           	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7A/) && ($p =~ /$Middle7C/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
         	
         	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7A/) && ($p =~ /$Ending7/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
         	
         	if (($p =~ /$Middle7A/) && ($p =~ /$Middle7B/) && ($p =~ /$Middle7C/))  {
          	 push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	 push (@PotentialLoci, "$p");
           	 $PopTester = 1;
           	 next;
            	}
         	
           	if (($p =~ /$Middle7A/) && ($p =~ /$Middle7B/) && ($p =~ /$Ending7/))  {
          	 push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	 push (@PotentialLoci, "$p");
           	 $PopTester = 1;
           	 next;
            	}
           	
            	if (($p =~ /$Middle7B/) && ($p =~ /$Middle7C/) && ($p =~ /$Ending7/))  {
           	 push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	 push (@PotentialLoci, "$p");
           	 $PopTester = 1;
           	 next;
           	} 
           
           	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7B/) && ($p =~ /$Middle7C/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
           	
           	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7B/) && ($p =~ /$Ending7/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
           	
           	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7C/) && ($p =~ /$Ending7/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
           	
           	if (($p =~ /$Middle7A/) && ($p =~ /$Middle7C/) && ($p =~ /$Ending7/))  {
           	 push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	 push (@PotentialLoci, "$p");
           	 $PopTester = 1;
           	 next;
           	} 
           	
          	 
      	}
      
      if ($PopTester == 0)  {  
           pop (@PotentialLoci);
           pop (@PotentialLoci);
           pop (@PotentialLoci);
      }
      
      $ElementNumber++;
  }

}  




elsif ($ReadLength > 80) {

 
   foreach (0..$NumOfUniqueNonErrorSeqsInTotalDataset-1)  {
  
   if ($_ =~ /000$/)  {
  	  
  	  print "Rapidly searching unique sequence number $_ of $NumOfUniqueNonErrorSeqsCounter for potentially allelic read pairs.\n";
   }
  
   push (@PotentialLoci, "Potential_Locus");
   push (@PotentialLoci, ">$UniqueNamesNonError[$ElementNumber]");
   push (@PotentialLoci, "$UniqueNonErrorSeqsInTotalDataset[$ElementNumber]");

  	
   my @SplitSeq = split (//, "$UniqueNonErrorSeqsInTotalDataset[$ElementNumber]");
  	my @Beginning7 = @SplitSeq[0..10];
  	my $Beginning7 = join "", @Beginning7;
  	my @Middle7A = @SplitSeq[20.30];
  	my $Middle7A = join "", @Middle7A;
  	my @Middle7B = @SplitSeq[31..41];
  	my $Middle7B = join "", @Middle7B;
  	my @Middle7C = @SplitSeq[42..52];
  	my $Middle7C = join "", @Middle7C;
  	my @Ending7 = @SplitSeq[-10..-1];
  	my $Ending7 = join "", @Ending7;
   

     	my $CountCycle = 0;					
     
     	my $ElementNumberB;	
     
     	my $PopTester = 0;
     
     	foreach my $p (@UniqueNonErrorSeqsInTotalDataset[($ElementNumber+1)..($NumOfUniqueNonErrorSeqsInTotalDataset-1)])  {
         
         	$CountCycle++;					
         
         	$ElementNumberB = $ElementNumber+$CountCycle;
         
         	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7A/) && ($p =~ /$Middle7B/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
           	
           	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7A/) && ($p =~ /$Middle7C/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
         	
         	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7A/) && ($p =~ /$Ending7/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
         	
         	if (($p =~ /$Middle7A/) && ($p =~ /$Middle7B/) && ($p =~ /$Middle7C/))  {
          	 push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	 push (@PotentialLoci, "$p");
           	 $PopTester = 1;
           	 next;
            	}
         	
           	if (($p =~ /$Middle7A/) && ($p =~ /$Middle7B/) && ($p =~ /$Ending7/))  {
          	 push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	 push (@PotentialLoci, "$p");
           	 $PopTester = 1;
           	 next;
            	}
           	
            	if (($p =~ /$Middle7B/) && ($p =~ /$Middle7C/) && ($p =~ /$Ending7/))  {
           	 push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	 push (@PotentialLoci, "$p");
           	 $PopTester = 1;
           	 next;
           	} 
           
           	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7B/) && ($p =~ /$Middle7C/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
           	
           	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7B/) && ($p =~ /$Ending7/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
           	
           	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7C/) && ($p =~ /$Ending7/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
           	
           	if (($p =~ /$Middle7A/) && ($p =~ /$Middle7C/) && ($p =~ /$Ending7/))  {
           	 push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	 push (@PotentialLoci, "$p");
           	 $PopTester = 1;
           	 next;
           	}  
      	}
      
      if ($PopTester == 0)  {  
           pop (@PotentialLoci);
           pop (@PotentialLoci);
           pop (@PotentialLoci);
      }
      
      $ElementNumber++;
  }

}




else  {

	
   foreach (0..$NumOfUniqueNonErrorSeqsInTotalDataset-1)  {
  
   if ($_ =~ /000$/)  {
  	  
  	  print "Rapidly searching unique sequence number $_ of $NumOfUniqueNonErrorSeqsCounter for potentially allelic read pairs.\n";
   }
  
   push (@PotentialLoci, "Potential_Locus");
   push (@PotentialLoci, ">$UniqueNamesNonError[$ElementNumber]");
   push (@PotentialLoci, "$UniqueNonErrorSeqsInTotalDataset[$ElementNumber]");

  	
   
   	my @SplitSeq = split (//, "$UniqueNonErrorSeqsInTotalDataset[$ElementNumber]");
  	my @Beginning7 = @SplitSeq[0..8];
  	my $Beginning7 = join "", @Beginning7;
  	my @Middle7A = @SplitSeq[15.23];
  	my $Middle7A = join "", @Middle7A;
  	my @Middle7B = @SplitSeq[28..36];
  	my $Middle7B = join "", @Middle7B;
  	my @Middle7C = @SplitSeq[42..50];
  	my $Middle7C = join "", @Middle7C;
  	my @Ending7 = @SplitSeq[-9..-1];
  	my $Ending7 = join "", @Ending7;
   
   
   
   
   

     	my $CountCycle = 0;					
     
     	my $ElementNumberB;	
     
     	my $PopTester = 0;
     
     	foreach my $p (@UniqueNonErrorSeqsInTotalDataset[($ElementNumber+1)..($NumOfUniqueNonErrorSeqsInTotalDataset-1)])  {
         
         	$CountCycle++;					
         
         	$ElementNumberB = $ElementNumber+$CountCycle;
         
         	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7A/) && ($p =~ /$Middle7B/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
           	
           	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7A/) && ($p =~ /$Middle7C/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
         	
         	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7A/) && ($p =~ /$Ending7/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
         	
         	if (($p =~ /$Middle7A/) && ($p =~ /$Middle7B/) && ($p =~ /$Middle7C/))  {
          	 push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	 push (@PotentialLoci, "$p");
           	 $PopTester = 1;
           	 next;
            	}
         	
           	if (($p =~ /$Middle7A/) && ($p =~ /$Middle7B/) && ($p =~ /$Ending7/))  {
          	 push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	 push (@PotentialLoci, "$p");
           	 $PopTester = 1;
           	 next;
            	}
           	
            	if (($p =~ /$Middle7B/) && ($p =~ /$Middle7C/) && ($p =~ /$Ending7/))  {
           	 push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	 push (@PotentialLoci, "$p");
           	 $PopTester = 1;
           	 next;
           	} 
           
           	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7B/) && ($p =~ /$Middle7C/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
           	
           	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7B/) && ($p =~ /$Ending7/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
           	
           	if (($p =~ /$Beginning7/) && ($p =~ /$Middle7C/) && ($p =~ /$Ending7/)) {
          	push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	push (@PotentialLoci, "$p");
           	$PopTester = 1;
           	next;
           	}
           	
           	if (($p =~ /$Middle7A/) && ($p =~ /$Middle7C/) && ($p =~ /$Ending7/))  {
           	 push (@PotentialLoci, "$UniqueNamesNonError[$ElementNumberB]");
           	 push (@PotentialLoci, "$p");
           	 $PopTester = 1;
           	 next;
           	}  
      	}
      
      if ($PopTester == 0)  {  
           pop (@PotentialLoci);
           pop (@PotentialLoci);
           pop (@PotentialLoci);
      }
      
      $ElementNumber++;
  }

}




open POTENTIAL_LOCI, ">TempFiles/PotentialLoci.txt" or die$!;

foreach my $l (@PotentialLoci)  {
   print POTENTIAL_LOCI "$l\n";
   } 
   
print POTENTIAL_LOCI "Potential_Locus";        
        
close POTENTIAL_LOCI;




############################################################

#The file "PotentialLoci" contains each of the (non-error) reads that may align with at least one other read in the dataset.
#Underneath each read are all of the reads that may be allelic with that read.
#Each original read (reads in "UniqueSeqsInTotalDataset" array) is indicated with a >>, and is preceded by "New_Potential_Locus" on the preceeding line.  
#This script goes through the PotentialLoci file and prints the relevant pairwise comparisons for each potential locus to a new file "TempFasta.txt".
#This TempFasta.txt file stores all of the pairwise comparisons that will be aligned by ACANA (in Fasta format, which is required by ACANA).

open POTLOCI, "TempFiles/PotentialLoci.txt" or die$!;
open TEMPFASTA, ">TempFiles/TempFasta.txt" or die$!;

my @PotLocus = ();

my $FirstCycle = 1;

while (<POTLOCI>)  {
   
   if ($FirstCycle == 1)  {
   $FirstCycle = 0;
   next;
   }
   
    if ($_ =~ /Pot/) {
      my $NumElements = @PotLocus;							
      my $NumComparisonsToWrite = ($NumElements/2)-1;	
      
      my $ElementToPrint = 2;
      
      for (0..($NumComparisonsToWrite-1))  {
         print TEMPFASTA "$PotLocus[0]";
         print TEMPFASTA "$PotLocus[1]";
         print TEMPFASTA "$PotLocus[$ElementToPrint]";
         $ElementToPrint++;
         print TEMPFASTA "$PotLocus[$ElementToPrint]";
         $ElementToPrint++;
      }
         @PotLocus = ();
         next;
    }

       
   if ($_ =~ />/)  {
      push (@PotLocus, $_);
      next;
   }
      
   if ($_ =~ /[^ACTGN]/)  {
      push (@PotLocus, $_);
      next;
   }
}

###############################################################################################

#TempFasta.txt file is not truly Fasta, as some names start with a >>.  
#Find all >> and replace with >.
#Print everything to file "TempFastaB.txt".

open TEMPFASTA, "TempFiles/TempFasta.txt" or die$!;
open TEMPFASTAB, ">TempFiles/TempFastaB.txt" or die$!;

while (<TEMPFASTA>)  {
  $_ =~ s/>>/>/;
  print TEMPFASTAB "$_"
  }

close TEMPFASTA;
close TEMPFASTAB;


###############################################################################

#Create an array that contains each of the lines in the TempFastaB.txt file.
#Note that this contains 4 elements for each of the pairwise comparisons - two names and two sequence reads.      

open TEMPFASTAB, "TempFiles/TempFastaB.txt" or die$!;

my @ForAlignment = ();
my %AlignmentHash;

while (<TEMPFASTAB>)  {
  push (@ForAlignment, $_)
	
  }

#Get the total number of alignments that will be done in ACANA.
  
my $ComparisonElements = @ForAlignment;
my $Comparisons = $ComparisonElements/4;

print "\nNumber of pairwise comparisons to align and evaluate is $Comparisons\n";
close TEMPFASTAB;
#######################################################################

#Print each pairwise comparison for alignment to a temporary file "SeqsToAlign.txt" in Fasta format.  This temporary file is read and analyzed by ACANA.
#The output from each ACANA run is added to the end of the file "ACANAOutfileA.txt".

print "Aligning all potentially alleleic read pairs with ACANA.\n";

my $CurrentAlignmentNumber = 0;

open ACANAOUTFILEA, ">TempFiles/ACANAOutfileA.txt";

for (0..$Comparisons-1)  {

	$CurrentAlignmentNumber++;

	if ($CurrentAlignmentNumber =~ /000$/)  {
 	 print "Working on alignment $CurrentAlignmentNumber of $Comparisons.\n";
	}

    open SEQSTOALIGN, ">TempFiles/SeqsToAlign.txt" or die$!;
    
    print SEQSTOALIGN "$ForAlignment[0]";
    print SEQSTOALIGN "$ForAlignment[1]";
    print SEQSTOALIGN "$ForAlignment[2]";
    print SEQSTOALIGN "$ForAlignment[3]";
    
    system "./ACANA  -I TempFiles/SeqsToAlign.txt -O TempFiles/ACANAOut.txt" ;
    
    
    open TEMPALIGNFILE, "TempFiles/ACANAOut.txt";
    
    my $AlignFileCounter = 0;
    while (<TEMPALIGNFILE>)  {
      if ($AlignFileCounter < 4)  {
      print ACANAOUTFILEA "$_";
      $AlignFileCounter++;
      }
    }
    close TEMPALIGNFILE;

	shift(@ForAlignment);
	shift(@ForAlignment);
	shift(@ForAlignment);
	shift(@ForAlignment);
}

close ACANAOUTFILEA;  
###############################################################################  
#Clean up the alignments

open READFILE, "TempFiles/ACANAOutfileA.txt";  
open ACANAOUTFILE, ">TempFiles/FinalAlignments.txt";    

while (<READFILE>)  {

    if ($_ =~ /\t/) {
    	    my @TempTabArray = split (/\t/, $_);
    	    my $ElementToPrint = $TempTabArray[0];
    	    print ACANAOUTFILE "$ElementToPrint\n";
    }
    
    else {
    	    print ACANAOUTFILE "$_"
    }
}
close ACANAOUTFILE;
close READFILE;

###############################################################################

open FILE, "TempFiles/FinalAlignments.txt" or die$!;

my @Alignments = ();

while (<FILE>)  {
   push (@Alignments, $_)
}

#####################################################################################
#Score alignments - any sites where there is a dash are ignored in the scoring.

my $NumInAlignments = @Alignments;

my $Counter5 = 0;

my @AlignmentScores = ();  
my @AlignmentLengths = ();


for (0..$Comparisons-1)  {


    my $Score = 0;
    my $ElementToCompare = 0;

    my @SplitSeq1 = split (//, $Alignments[($Counter5+1)]);
    my $CurrentSeq1 = $Alignments[($Counter5)];
    my @SplitSeq2 = split (//, $Alignments[($Counter5+3)]);
    my $CurrentSeq2 = $Alignments[($Counter5+2)];
    
    my $AlignmentLength1 = @SplitSeq1;
    my $AlignmentLength2 = @SplitSeq2;
    
    $AlignmentLength1 = $AlignmentLength1-1;
    $AlignmentLength2 = $AlignmentLength2-1;
     
    foreach my $site (@SplitSeq1)  {
         
           if (($site eq "-") | ($SplitSeq2[$ElementToCompare] eq "-"))  {
             $AlignmentLength1--;
             $AlignmentLength2--;
             $ElementToCompare++;
             next;
           }   
         
           elsif (($SplitSeq2[$ElementToCompare]) eq $site)  {
             $Score++;
             $ElementToCompare++;
           }
           
           else {
             $ElementToCompare++;
             }
     }
           
     my $CalculatedScore = ($Score/$AlignmentLength1)*100;
     push (@AlignmentScores, $CalculatedScore);
     push (@AlignmentLengths, $AlignmentLength1);
      
     $Counter5++;
     $Counter5++;
     $Counter5++;
     $Counter5++;
}
 
#####################################################################################


















#Go through scored alignments and discard any that have a score below a given value that is input by the user.
#Pairwise alignments that exceed the threshold % identity are printed to GoodAlignments.txt.

 open GOODALIGNMENTS, ">TempFiles/GoodAlignments.txt" or die $!;
 
 print "\nScoring each pairwise alignment.  Alignments exceeding $AllelicCriticalValue % similarity will be retained.\n";
 
 my $Counter6 = 0;
 
 my @GoodAlignmentScores = ();
 my @GoodAlignments = ();
 my $AlignmentLengthElement = 0;
 
 foreach my $score (@AlignmentScores)  {
    
    if (($score >= $AllelicCriticalValue) && ($AlignmentLengths[$AlignmentLengthElement] >= $MinAlignmentLength))  {
     push (@GoodAlignmentScores, $score);
     push (@GoodAlignments, $Alignments[$Counter6]);
     push (@GoodAlignments, $Alignments[$Counter6+1]);
     push (@GoodAlignments, $Alignments[$Counter6+2]);
     push (@GoodAlignments, $Alignments[$Counter6+3]);
     print GOODALIGNMENTS "$Alignments[$Counter6]";
     print GOODALIGNMENTS "$Alignments[$Counter6+1]";
     print GOODALIGNMENTS "$Alignments[$Counter6+2]";
     print GOODALIGNMENTS "$Alignments[$Counter6+3]";
     $Counter6++;
     $Counter6++;
     $Counter6++;
     $Counter6++;
     $AlignmentLengthElement++;
    }
 
    else {
     $Counter6++;
     $Counter6++;
     $Counter6++;
     $Counter6++;
     $AlignmentLengthElement++;
    }
    
 }
 
 my $NumGoodAlignmentElements = @GoodAlignments;
 my $NumGoodAlignments = $NumGoodAlignmentElements/4;
 my $NumGoodAlignmentScores = @GoodAlignmentScores;

 close GOODALIGNMENTS;
######################################################################################
######################################################################################
#Call Loci
#Paired, aligned reads are put together at the same locus, and any other pairs that have at least one read in common with a previous pair are appended to the locus.

open GOODALIGNMENTS, "TempFiles/GoodAlignments.txt" or die$!;

print "Parsing aligned sequences into candidate loci.\n";

my @TempArray = ();
my $WhileLoopCounter = 0;
my %UsedNamesAndCandidateLocusNumbers = ();
my $NextCandidateLocusNumber = 0;
my $GoodAlignmentsLineCounter = 0;

while (<GOODALIGNMENTS>)  {
	
   $GoodAlignmentsLineCounter++;  
   chomp ($_);
  
     push (@TempArray, "$_");
     $WhileLoopCounter++;
  
  if ($WhileLoopCounter == 4)  { 
    
    if ($UsedNamesAndCandidateLocusNumbers{$TempArray[0]})  {	#First name in pair is defined (has been used already)
       if (exists $UsedNamesAndCandidateLocusNumbers{$TempArray[2]})  {
	 $WhileLoopCounter = 0;
         @TempArray = ();
         next;
       }
      

       else {

	 my $CurrentFileToPrintTo = $UsedNamesAndCandidateLocusNumbers{$TempArray[0]};
         open TEMPCANDIDATELOCUSFILE, ">>CandidateLoci/CandidateLocus$CurrentFileToPrintTo.txt" or die$!;
	 print TEMPCANDIDATELOCUSFILE "$TempArray[2]\n";
         print TEMPCANDIDATELOCUSFILE "$TempArray[3]\n";
         $UsedNamesAndCandidateLocusNumbers{$TempArray[2]} = "$CurrentFileToPrintTo";
         close TEMPCANDIDATELOCUSFILE;
         $WhileLoopCounter = 0;
         @TempArray = ();
         next; 	 
       }  
    }
    
    else  {   #First name in pair is undefined (hasn't been used yet)
       if ($UsedNamesAndCandidateLocusNumbers{$TempArray[2]})  {
	 my $CurrentFileToPrintTo = $UsedNamesAndCandidateLocusNumbers{$TempArray[2]};
         open TEMPCANDIDATELOCUSFILE, ">>CandidateLoci/CandidateLocus$CurrentFileToPrintTo.txt" or die$!;
         print TEMPCANDIDATELOCUSFILE "$TempArray[0]\n";
         print TEMPCANDIDATELOCUSFILE "$TempArray[1]\n";
         $UsedNamesAndCandidateLocusNumbers{$TempArray[0]} = "$CurrentFileToPrintTo";
         close TEMPCANDIDATELOCUSFILE;
         $WhileLoopCounter = 0;
         @TempArray = ();
         next;
       }
     
       else {
	 $NextCandidateLocusNumber++;
         open TEMPCANDIDATELOCUSFILE, ">CandidateLoci/CandidateLocus$NextCandidateLocusNumber.txt" or die$!;
         print TEMPCANDIDATELOCUSFILE "$TempArray[0]\n";
         print TEMPCANDIDATELOCUSFILE "$TempArray[1]\n";
         print TEMPCANDIDATELOCUSFILE "$TempArray[2]\n";
         print TEMPCANDIDATELOCUSFILE "$TempArray[3]\n";
         $UsedNamesAndCandidateLocusNumbers{$TempArray[0]} = "$NextCandidateLocusNumber";
         $UsedNamesAndCandidateLocusNumbers{$TempArray[2]} = "$NextCandidateLocusNumber";
         close TEMPCANDIDATELOCUSFILE;
         $WhileLoopCounter = 0;
         @TempArray = ();
         next;
       }
     } 
     
   }
}   




###################################################################################### 

system "rm TempFiles/SortedIndividual*.txt";
system "rm TempFiles/Individual*";

###################################################################################### 
#Realign all loci (global alignments - there may be more than 2 alleles at some loci).  This is why mafft is used instead of ACANA, as ACANA only does pairwise alignments.
#This is done with mafft.

 if (-d "CandidateLociRealigned") {
 	 system "rm -r CandidateLociRealigned";
 }
 
 mkdir "CandidateLociRealigned";
 
 
 my $GlobalAlignmentCounter = 0;
 
 foreach my $LocusNumber (1..$NextCandidateLocusNumber) {
 	 
 	 $GlobalAlignmentCounter++;
 	 
 	 if ($GlobalAlignmentCounter =~ /00$/)  {
 	 	 print "Globally aligning locus $GlobalAlignmentCounter of $NextCandidateLocusNumber.\n";
 	 }	 
	 
	 open CURRENTLOCUSFILE, "CandidateLoci/CandidateLocus$LocusNumber.txt" or die$!;
	 open TEMPFILETOREALIGN, ">TempFiles/SeqsToReAlign.txt" or die$!;
	 
	 while (<CURRENTLOCUSFILE>)  {
		 chomp($_);
		 
		 if ($_ =~ />/)  {
			 
			 print TEMPFILETOREALIGN "$_\n";
			 next;
		 }
       
		 else {
			 $_ =~ s/a/A/g;
			 $_ =~ s/c/C/g;
			 $_ =~ s/g/G/g;
			 $_ =~ s/t/T/g;
			 $_ =~ s/n/N/g;
			 $_ =~ s/-//g;
			 
			 print TEMPFILETOREALIGN "$_\n";
		 }	
	}
	 
	 system "mafft --quiet TempFiles/SeqsToReAlign.txt > CandidateLociRealigned/mafftReAlignOut.txt" ;
	 
	 open REALIGNEDTOPRINT, "CandidateLociRealigned/mafftReAlignOut.txt" or die$!;
	 open REALIGNED, ">CandidateLociRealigned/CandidateLocus$LocusNumber.txt" or die$!;
	 
	 while (<REALIGNEDTOPRINT>)  {
		 print REALIGNED "$_";
	 }
 }	
		
 close CURRENTLOCUSFILE;
 close TEMPFILETOREALIGN;
 close REALIGNEDTOPRINT;
 close REALIGNED;
	 
# ######################################################################################
if (-d "Output/PolymorphicLoci") {
	system "rm -r Output/PolymorphicLoci";
}

mkdir "Output/PolymorphicLoci";
 
 foreach my $LocusNumber (1..$NextCandidateLocusNumber) {
	 
	 open CURRENTLOCUSFILE, "CandidateLociRealigned/CandidateLocus$LocusNumber.txt" or die$!;
	 open WRITEFILE, ">Output/PolymorphicLoci/CandidateLocus$LocusNumber.txt" or die$!;
	 
	 my $Seq1 = 1;
	 
	 while (<CURRENTLOCUSFILE>)  {
		 chomp($_);
		 
		 if ($Seq1 == 1)  {
			 print WRITEFILE "$_\n";
			 $Seq1 = 0;
		 }	
			 
		 else {	
		 
			 if ($_ =~ />/)  {
				 print WRITEFILE "\n$_\n";
			 }	
	 
			 else {
				 print WRITEFILE "$_";
			 }
		 }
	 }	
 }	
 
 print "Finished writing the mafft output.\n";
 
 close CURRENTLOCUSFILE;
 close WRITEFILE;
 
 #system "rm CandidateLoci/CandidateLocus*";
 system "rm -r CandidateLoci";
 
 #system "rm CandidateLociRealigned/CandidateLocus*";
 #system "rm CandidateLociRealigned/ma*";
 system "rm -r CandidateLociRealigned";
 
 
#######################################################################################

#Generate raw genotype files for each individual
#These files contain all of the polymorphic loci in the dataset, and all of the alleles identified at each locus.  They also contain the read counts for each of these alleles for the respective individual.
#Note that the read counts are the only things that vary among these files.

my %TempSeqsAndCounts;

foreach my $a (@AllSampleNames)  {
	
	
	if (-e "TempFiles/UniqueWithCountsIndividual$a.txt")  {
	
		%TempSeqsAndCounts = ();

		my @TempSeqsAndCounts = ();

		open COUNTSANDALLELES, "TempFiles/UniqueWithCountsIndividual$a.txt" or die$!;
    
		open INDIVIDUALGENOTYPESFILE, ">TempFiles/TempRawCounts$a.txt" or die$!;
    
   
		while (<COUNTSANDALLELES>)  {
  	 
			if ($_ =~ /[a-zA-Z]/) {
    	    
				chomp ($_);
				$_ =~ s/^\s*//;
				my @SplitBetweenCountAndSeq = split (/\s/, $_);
				push (@TempSeqsAndCounts, "$SplitBetweenCountAndSeq[1]");
				push (@TempSeqsAndCounts, "$SplitBetweenCountAndSeq[0]");
  	    
    
			}
		}

		%TempSeqsAndCounts = @TempSeqsAndCounts;
	

		my $size = keys %TempSeqsAndCounts;
 
		foreach my $b (1..$NextCandidateLocusNumber)  {
 	
			if ($b == 1)  {
    	    
    	    
				print INDIVIDUALGENOTYPESFILE "LocusNumber$b";
				open CURRENTLOCUSFILE, "Output/PolymorphicLoci/CandidateLocus$b.txt" or print "Did not find CandidateLocus$b\n";
    
				while (<CURRENTLOCUSFILE>) {
    	    	    	  
					if ($_ =~ />/)  {
						next;
    	    	    	    	    	}
      
    	    	    	    	    	else {
    	    	    	    	    		my $TempSeqWithGaps;
    	    	    	    	    		my $TempSeq = "$_";
    	    	    	    	    		$TempSeq =~ s/a/A/g;
    	    	    	    	    		$TempSeq =~ s/c/C/g;
    	    	    	    	    		$TempSeq =~ s/g/G/g;
    	    	    	    	    		$TempSeq =~ s/t/T/g;
    	    	    	    	    		$TempSeq =~ s/n/N/g;
    	    	    	    	    		$TempSeqWithGaps = "$TempSeq";
    	    	    	    	    		chomp ($TempSeqWithGaps);
    	    	    	    	    		$TempSeq =~ s/-//g;
    	    	    	    	    		chomp ($TempSeq);
    	    	    	    	    		print INDIVIDUALGENOTYPESFILE "\n$TempSeqWithGaps\t";
      
    	    	    	    	    		if (exists $TempSeqsAndCounts{$TempSeq})  {
    	    	    	    	    			print INDIVIDUALGENOTYPESFILE "$TempSeqsAndCounts{$TempSeq}";
    	    	    	    	    		}
    	    	    	    	    	
    	    	    	    	    		else {
    	    	    	    	    			print INDIVIDUALGENOTYPESFILE "0";
    	    	    	    	    		}
      
    	    	    	    	    	}
    	    	    	    	}
    	    	    	}
    
    
    	    	    	else {
    	    	    		print INDIVIDUALGENOTYPESFILE "\nLocusNumber$b";
    	    	    		open CURRENTLOCUSFILE, "Output/PolymorphicLoci/CandidateLocus$b.txt" or print "Did not find CandidateLocus$b\n";
    
    	    	    		while (<CURRENTLOCUSFILE>) {
    	    	    			if ($_ =~ />/)  {
    	    	    	     	     next;
    	    	    	     	     	}
      
    	    	    	     	     	else {
    	    	    	     	     		my $TempSeqWithGaps;
    	    	    	     	     		my $TempSeq = "$_";
    	    	    	     	     		$TempSeq =~ s/a/A/g;
    	    	    	     	     		$TempSeq =~ s/c/C/g;
    	    	    	     	     		$TempSeq =~ s/g/G/g;
    	    	    	     	     		$TempSeq =~ s/t/T/g;
    	    	    	     	     		$TempSeq =~ s/n/N/g;
    	    	    	     	     		$TempSeqWithGaps = "$TempSeq";
    	    	    	     	     		chomp ($TempSeqWithGaps);
    	    	    	     	     		$TempSeq =~ s/-//g;
    	    	    	     	     		chomp ($TempSeq);
    	    	    	     	     		print INDIVIDUALGENOTYPESFILE "\n$TempSeqWithGaps\t";
      
    	    	    	     	     		if (exists $TempSeqsAndCounts{$TempSeq})  {
    	    	    	     	     			print INDIVIDUALGENOTYPESFILE "$TempSeqsAndCounts{$TempSeq}";
    	    	    	     	     		}
      
    	    	    	     	     		else {
    	    	    	     	     			print INDIVIDUALGENOTYPESFILE "0";
    	    	    	     	     		}
    
    	    	    	     	     	}
    	    	    	     	}
    	    	    	}    
 

    	    	    	close COUNTSANDALLELES;
 
    	    	    	close CURRENTLOCUSFILE;

    	    	}	

    	    	close INDIVIDUALGENOTYPESFILE;
	}  	
}




##################################################################################
#Edit genotype files

if (-d "TempFiles/RawReadCountFiles")  {
	system "rm -r TempFiles/RawReadCountFiles";
}		

mkdir "TempFiles/RawReadCountFiles";

foreach my $a (@AllSampleNames)  {

	
	if (-e "TempFiles/TempRawCounts$a.txt")  {
		
		open TEMPFILE, "TempFiles/TempRawCounts$a.txt" or die$!;
		open TEMPFILETWO, ">TempFiles/RawReadCountFiles/RawReadCounts_$a.txt" or die$!;
		while (<TEMPFILE>)  {
         
			if ($_ =~ /^\t/)  {
				next;
			}

			else {
				print TEMPFILETWO "$_";
			}
		}
		close TEMPFILE;
		close TEMPFILETWO;
	}	
      
}

system "rm TempFiles/TempRaw*";


#######################################################################################################
#Remove paralogs and print alleles two at a time to files for binomial test in R.

print "Identifying and removing paralogous loci.\n";


#Get the names of paralogous loci and store them in array ParalogousLociNames

my @ParalogousLociNames = ();
my @ParalogousLociNamesAltFixed = ();	#These are duplicated loci that are alternatively fixed within the genome - they don't get included in @ParalogousLocusNames, as they don't have a third allele.
my %PotentialAltFixedParalogousNames = ();
my %NumIndividualsThatCanBeScoredForAltFixedTest = ();


foreach my $name (@AllSampleNames)  {

  if (-e "TempFiles/RawReadCountFiles/RawReadCounts_$AllSampleNames[0].txt")  {
	
  	  open FILE, "TempFiles/RawReadCountFiles/RawReadCounts_$AllSampleNames[0].txt";
	
  	  while(<FILE>)  {
		
  	  	  if ($_ =~ /Locus/)  {
  	  	  	  
  	  	  	my @Locus = split (/r/, $_);  
			my $LocusNumber = $Locus[1];
			chomp($LocusNumber);
			$PotentialAltFixedParalogousNames{$LocusNumber} = 0;
			$NumIndividualsThatCanBeScoredForAltFixedTest{$LocusNumber} = 0;
		}
	}
	close FILE;
	
	last;
  }
}  





foreach my $name (@AllSampleNames)  {		#Go through each individual in the dataset
	
	if (-e 	"TempFiles/RawReadCountFiles/RawReadCounts_$name.txt")  {
	
	
		open READFILE, "TempFiles/RawReadCountFiles/RawReadCounts_$name.txt" or die$!;
	
		
		my $LocusNumber;
		my @TempArrayCounts;
		my @TempArraySeqs;
		my $Starter = 0;
		my $LastAnalyzed = 0;
   
		while (<READFILE>)  {	#open while 1
        
			if ($_ =~ /[a-zA-Z]/) {	
	
	
				if ($Starter == 1)  {  #open if Starter
        	
					if ($_ =~ /Locus/)  {  #open if Locus		#Have reached a new locus in the current individual.  Arrays TempArraySeqs and TempArrayCounts are currently populated with seqs and counts from the previous locus.  Need to analyze these.						
						
						my @Locus = split (/r/, $_);
						$LocusNumber = ($Locus[1]-1);
						
						$LastAnalyzed = 1;
							
            
						my @SortedCounts = sort {$b <=> $a} @TempArrayCounts; 	#If more than two alleles, sort the values of the counts to get the third highest count.
            
						#First, get sum of count of top two reads
						
						my $SumCount = $SortedCounts[0]+$SortedCounts[1];
						
						if ($SumCount > 25)  {
							$NumIndividualsThatCanBeScoredForAltFixedTest{$LocusNumber}++;
						}	
						
						#Then, check the 2nd highest count for potentially alternatively fixed paralogous locus
						
						my $SecondCount = $SortedCounts[1];
						if (($SecondCount >= $MinForParalog) && ($SumCount >20))  {
							$PotentialAltFixedParalogousNames{$LocusNumber}++;
						}	
						
						
						
						my $NumberInArray = @TempArrayCounts;	#Determine how many unique alleles are at the locus being analyzed.
            
						if ($NumberInArray == 2) {			#If just two alleles, can move on to next locus.
							@TempArraySeqs = ();
							@TempArrayCounts = ();
							next;
						}
						
						
						
						my $ThirdCount = $SortedCounts[2];
						if ($ThirdCount >= $MinForParalog)  {			#If the third highest count at this locus in the current individual exceeds the minimum threshold, locus is flagged as paralogous.
            	  
							my $LocusName = "LocusNumber$LocusNumber";
							push (@ParalogousLociNames, $LocusName);
							@TempArraySeqs = ();
							@TempArrayCounts = ();
							next;  
						}
          
						else {
							@TempArraySeqs = ();
							@TempArrayCounts = ();
							next;
						}    
					}
          
          
					else {
						my @SplitLine = split (/\t/, $_);			#Get the sequences and counts for the current locus in the current individual pushed to arrays.
						push (@TempArraySeqs, "$SplitLine[0]");
						push (@TempArrayCounts, "$SplitLine[1]");
						$LastAnalyzed = 0;
					}  
				} #close if Starter = 1
       	
     
        
				else {
					$Starter = 1;
					next;
				}
			}	    
		}


		if ($LastAnalyzed == 0)  {					#Analyze the last locus in the file.
			$LocusNumber++;
    	
			
				my @SortedCounts = sort {$b <=> $a} @TempArrayCounts; 
            
				#First, get sum of count of top two reads
						
				my $SumCount = $SortedCounts[0]+$SortedCounts[1];
						
				if ($SumCount > 25)  {
					$NumIndividualsThatCanBeScoredForAltFixedTest{$LocusNumber}++;
				}
				
				
				#Then, check the 2nd highest count for potentially alternatively fixed paralogous locus
						
				my $SecondCount = $SortedCounts[1];
					if ($SecondCount >= $MinForParalog)  {
					$PotentialAltFixedParalogousNames{$LocusNumber}++;
				}
				
				
				
				my $NumberInArray = @TempArrayCounts;	#Determine how many unique alleles are at the locus being analyzed.
            
				if ($NumberInArray == 2) {			#If just two alleles, can move on to next locus.
					@TempArraySeqs = ();
					@TempArrayCounts = ();
					next;
				}
				
				
				
				my $ThirdCount = $SortedCounts[2];
				if ($ThirdCount >= $MinForParalog)  {
					my $LocusName = "LocusNumber$LocusNumber";
					push (@ParalogousLociNames, $LocusName);
					@TempArraySeqs = ();
					@TempArrayCounts = ();
				}
			#}    
		}                

		close READFILE;
  	
  	}	
}


#Check the PotentialAltFixedParalogousNames hash for loci that have counts greater than some min proportion of individuals (maybe around 90%).




for (keys %PotentialAltFixedParalogousNames) {
	
	
	#print "$_\t$PotentialAltFixedParalogousNames{$_}\n";
	#Get the total number of individuals that can be genotyped at this locus
	
	my $NumIndividualsThatCanBeScored = $NumIndividualsThatCanBeScoredForAltFixedTest{$_};
	my $CurrentThresholdNumIndividuals = $NumIndividualsThatCanBeScored * $MaxProportionHeterozygotes;
	
	
	
	if (($PotentialAltFixedParalogousNames{$_} >= $CurrentThresholdNumIndividuals) && ($CurrentThresholdNumIndividuals > 5))  {
		push (@ParalogousLociNamesAltFixed, $_);
	}
}

open ALTFIXED, ">TempFiles/AltFixedParalogs.txt" or die$!;

#Merge the arrays ParalogousLociNames and ParalogousLociNamesAltFixed;

my $NumInParalogousLociNamesAltFixed = @ParalogousLociNamesAltFixed;


foreach my $name (@ParalogousLociNamesAltFixed)  {
	
	my $Match = 0;
	
	foreach my $name2 (@ParalogousLociNames)  {
		
		if ($name2 =~ /[A-Za-z]/) {
			my @TempArray = split(/r/, $name2);
			my $TempName = $TempArray[1];
		
			if ($name eq $TempName)  {
				$Match = 1;
				last;
			}
		}	
	}
	
	if ($Match == 0)  {
		my $NewParalogousLocus = 'LocusNumber'.$name;
		push (@ParalogousLociNames, $NewParalogousLocus);
		print ALTFIXED "LocusNumber$name\n";
	}
	
}	

close ALTFIXED;

#Print all nonparalogous loci to file RawReadCounts_NonParalogous

foreach my $name (@AllSampleNames)  {
	
	
	if (-e "TempFiles/RawReadCountFiles/RawReadCounts_$name.txt")  {


		open READFILE, "TempFiles/RawReadCountFiles/RawReadCounts_$name.txt" or die$!;
		open WRITEFILE, ">TempFiles/RawReadCountFiles/RawReadCounts_NonParalogous$name.txt" or die$!;
		open PARALOGOUSLOCI, ">TempFiles/ParalogousLoci.txt" or die$!;

		my @TempArray = ();
		my $Starter = 0;
		my $CurrentLocus;
		my $Match = 0;
		my $LastAnalyzed = 0;

		while (<READFILE>)  {
			chomp($_);
	
			if ($Starter == 1)  {
	
				if ($_ =~ /Locus/) {
					$Match = 0;	
					$LastAnalyzed = 1;
	        
					foreach my $LocusName (@ParalogousLociNames)  {
			
						if ($CurrentLocus eq $LocusName)  {
							$Match = 1;
	
				
							print PARALOGOUSLOCI "$CurrentLocus\n";
				
							foreach my $element (@TempArray)  {
								print PARALOGOUSLOCI "$element\n";
							}
				
							last;
						}	
					}
		
					if ($Match == 0)  {
						print WRITEFILE "$CurrentLocus\n";
			
						foreach my $element (@TempArray)  {
							print WRITEFILE "$element\n";
						}
					
					}
		

					$CurrentLocus = $_;
					@TempArray = ();
				}
	
	
	
				else {
					push (@TempArray, $_);
					$LastAnalyzed = 0;
					next;
				}	

			}
	
			elsif ($Starter == 0)  {
				$Starter = 1;
				$CurrentLocus = $_;
			}
		}


		if ($LastAnalyzed == 0)  {
	
			$Match = 0;	
	
			foreach my $LocusName (@ParalogousLociNames)  {
			
				if ($CurrentLocus eq $LocusName)  {
					$Match = 1;
				
					print PARALOGOUSLOCI "$CurrentLocus\n";
				
					foreach my $element (@TempArray)  {
						print PARALOGOUSLOCI "$element\n";
					}
				
					last;
				}	
			}
		
			if ($Match == 0)  {
				print WRITEFILE "$CurrentLocus\n";
			
				foreach my $element (@TempArray)  {
					print WRITEFILE "$element\n";
				}
			
			}
		
			$CurrentLocus = $_;
			@TempArray = ();
		}	


		close READFILE;	
		close WRITEFILE;
		close PARALOGOUSLOCI;

	}
}






#Create file "ForBinomialTest" for each individual that contains the alleles with the top two read counts for each locus for that individual.

foreach my $name (@AllSampleNames)  {

	if (-e 	"TempFiles/RawReadCountFiles/RawReadCounts_NonParalogous$name.txt") {
	
		open READFILE, "TempFiles/RawReadCountFiles/RawReadCounts_NonParalogous$name.txt" or die$!;
		open WRITEFILE, ">TempFiles/ForBinomialTest$name.txt" or die$!;


		my @TempSeqArray = ();
		my @TempCountArray = ();
		my $Starter = 0;
		my $CurrentLocus;
		my $PreviousLocusNumber;
		my $LastAnalyzed = 0;
		my $CurrentElementCounter;
	
		while (<READFILE>)  {
	
			if ($_ =~ /[a-zA-Z]/) {
	
				chomp($_);
	
				if ($Starter == 1)  {
	
					if ($_ =~ /Locus/) {
						$LastAnalyzed = 1;
						my $NumberOfElements = @TempCountArray;
	   	  
						if ($NumberOfElements == 2)  {	#Checks to see if the locus is biallelic.  If so, just print it as is.
							print WRITEFILE "LocusNumber$PreviousLocusNumber\tA\n";
							print WRITEFILE "$TempSeqArray[0]\t$TempCountArray[0]\n";
							print WRITEFILE "$TempSeqArray[1]\t$TempCountArray[1]\n";
							my @PreviousLocus = split (/r/, $_);
							$PreviousLocusNumber = $PreviousLocus[1];
							@TempSeqArray = ();
							@TempCountArray = ();
							next;
						}
	   	  
						else {	#The locus is not biallelic - find the top two counts and print these.
							my $MaxCount = 0;
							my $CurrentElementCounter = 0;
							my $MaxPosition = 0;
							my $AllZeros = 1;
	   	  	 
							foreach my $a (@TempCountArray)  {	
								if ($a > 0) {
									$AllZeros = 0;
								}
							}
	   	  	  
							if ($AllZeros == 1) {	
								print WRITEFILE "LocusNumber$PreviousLocusNumber\tA\n";
								print WRITEFILE "$TempSeqArray[0]\t$TempCountArray[0]\n";
								print WRITEFILE "$TempSeqArray[1]\t$TempCountArray[1]\n";
								my @PreviousLocus = split (/r/, $_);
								$PreviousLocusNumber = $PreviousLocus[1];
								@TempSeqArray = ();
								@TempCountArray = ();
								next;
							}	  
	   	  	  
	   	  	  
							#If the locus is not biallelic, and if it doesn't have all zeros, get the maximum value in the array; the value is in $MaxCount, and the position of this value in the array is $MaxPosition
							foreach my $element (@TempCountArray)  {	
								if ($element > $MaxCount)  {
									$MaxCount = $element;
									$MaxPosition = $CurrentElementCounter;
								}
								$CurrentElementCounter++;
	   	  	  	  
	   	  	 				}	  
	   	  	 
	   	  	  
	   	  	 				#Check to see if a second allele has the same count as the max allele.
	   	  	 				my $MaxCount2Equal = 0;	
	   	  	 				my $CurrentCounter = 0;
	   	  	 				my $TwoEqual = 0;
	   	  	  
	   	  	 				foreach my $b (@TempCountArray)  {
	   	  	 					if (($b == $MaxCount) && ($CurrentCounter != $MaxPosition))  {  #Have the same count for two alleles.
	   	  	 						my $IdenticalCountElement = $CurrentCounter;
	   	  	 						print WRITEFILE "LocusNumber$PreviousLocusNumber\tA\n";
	   	  	 						print WRITEFILE "$TempSeqArray[$MaxPosition]\t$TempCountArray[$MaxPosition]\n";
	   	  	 						print WRITEFILE "$TempSeqArray[$IdenticalCountElement]\t$TempCountArray[$IdenticalCountElement]\n";
	   	  	 						my @PreviousLocus = split (/r/, $_);
	   	  	 						$PreviousLocusNumber = $PreviousLocus[1];
	   	  	 						@TempSeqArray = ();
	   	  	 						@TempCountArray = ();
	   	  	 						$TwoEqual = 1;
	   	  	 						last;
	   	  	 					}	
	   	  	  	
	   	  	 					$CurrentCounter++;
	   	  	 				}
	   	  	  
	   	  	  
	   	  	 				if ($TwoEqual == 0)  {
	   	  	  	  
	   	  	 					#Get the second highest allele count and assign the count to $SecondMaxCount and the element position of this count to $SecondMaxPosition.
	   	  	 					my $SecondMaxCount = 0;	
	   	  	 					my $SecondMaxPosition = 0;
	   	  	 					my $CurrentElementCounter2 = 0;
	   	  	  
	   	  	 					foreach my $c (@TempCountArray)  {
	   	  	 						if (($c < $MaxCount) && ($c >= $SecondMaxCount)) {
	   	  	 							$SecondMaxCount = $c;
	   	  	 							$SecondMaxPosition = $CurrentElementCounter2;
	   	  	 						}
	   	  	 						$CurrentElementCounter2++;
	   	  	 					}
	   	  	  
	   	  	 					print WRITEFILE "LocusNumber$PreviousLocusNumber\tA\n";
	   	  	 					print WRITEFILE "$TempSeqArray[$MaxPosition]\t$TempCountArray[$MaxPosition]\n";
	   	  	 					print WRITEFILE "$TempSeqArray[$SecondMaxPosition]\t$TempCountArray[$SecondMaxPosition]\n";
	   	  	  	
	   	  	 				}
	   	  	 			}       
	   	  
	   	  	 			my @PreviousLocus = split (/r/, $_);
	   	  	 			$PreviousLocusNumber = $PreviousLocus[1];
	   	  	 			@TempSeqArray = ();
	   	  	 			@TempCountArray = ();
	   
	   	  	 		} 	  
	   	 
	   	  	 		
	   	  	 		
	   	  	 		else {
	   	  	 			my @SplitLine = split (/\t/, $_);
	   	  	 			push (@TempSeqArray, "$SplitLine[0]");
	   	  	 			push (@TempCountArray, "$SplitLine[1]");
	   	  	 			$LastAnalyzed = 0;
                 			} 
       	   
                 		}  
	   	   
	
                 		else {
                 			my @PreviousLocus = split (/r/, $_);
                 			$PreviousLocusNumber = $PreviousLocus[1];
           
                 			$Starter = 1;
                 			next;
                 		} 
                 	}	 
                }	 
	   	   

                
                
                if ($LastAnalyzed == 0)  {
	
                	my $NumberOfElements = @TempCountArray;
	   	  
                	if ($NumberOfElements == 2)  {	#Checks to see if the locus is biallelic.  If so, just print it as is.
	   	  	  print WRITEFILE "LocusNumber$PreviousLocusNumber\tA\n";
	   	  	  print WRITEFILE "$TempSeqArray[0]\t$TempCountArray[0]\n";
	   	  	  print WRITEFILE "$TempSeqArray[1]\t$TempCountArray[1]\n";
	   	  	}
	   	  
	   	  	else {	#The locus is not biallelic - find the top two counts and print these.
	   	  	  my $MaxCount = 0;
	   	  	  my $CurrentElementCounter = 0;
	   	  	  my $MaxPosition = 0;
	   	  	  my $AllZeros = 1;
	   	  	  
	   	  	  
	   	  	  #This prints the first two reads if all of the reads have count of zero.
	   	  	  foreach my $a (@TempCountArray)  {	
	   	  	  	  if ($a > 0) {
	   	  	  	  	  $AllZeros = 0;
	   	  	  	  }
	   	  	  }
	   	  	  
	   	  	  if ($AllZeros == 1) {	
	   	  	  	  print WRITEFILE "LocusNumber$PreviousLocusNumber\tA\n";
	   	  	  	  print WRITEFILE "$TempSeqArray[0]\t$TempCountArray[0]\n";
	   	  	  	  print WRITEFILE "$TempSeqArray[1]\t$TempCountArray[1]\n";
	   	  	  }	  
	   	  	  
	   	  	 
	   	  	  else {
	   	  	  	  
	   	  	  	  #If the locus is not biallelic, and if it doesn't have all zeros, get the maximum value in the array; the value is in $MaxCount, and the position of this value in the array is $MaxPosition
	   	  	  	  foreach my $element (@TempCountArray)  {	
	   	  	  	  	  if ($element > $MaxCount)  {
	   	  	  	  	  	  $MaxCount = $element;
	   	  	  	  	  	  $MaxPosition = $CurrentElementCounter;
	   	  	  	  	  }
	   	  	  	  	  $CurrentElementCounter++;
	   	  	  	  
	   	  	  	  }	  
	   	  	  
	   	  	  	  #Check to see if a second allele has the same count as the max allele.
	   	  	  	  my $MaxCount2Equal = 0;	
	   	  	  	  my $CurrentCounter = 0;
	   	  	  	  my $TwoEqual = 0;
	   	  	  
	   	  	  	  foreach my $b (@TempCountArray)  {
	   	  	  		if (($b == $MaxCount) && ($CurrentCounter =~ $MaxPosition))  {  #Have the same count for two alleles.
	   	  	 	 		my $IdenticalCountElement = $CurrentCounter;
	   	  	 	 		print WRITEFILE "LocusNumber$PreviousLocusNumber\tA\n";
	   	  	  			print WRITEFILE "$TempSeqArray[$MaxPosition]\t$TempCountArray[$MaxPosition]\n";
	   	  	  			print WRITEFILE "$TempSeqArray[$IdenticalCountElement]\t$TempCountArray[$IdenticalCountElement]\n";
	   	  	  			$TwoEqual = 1;
	   	  	  			
	   	  	  		}	
	   	  	  	
	   	  	  		$CurrentCounter++;
	   	  	  	  }
	   	  	  
	   	  	  	  if ($TwoEqual == 0)  {
	   	  	  	
	   	  	  
	   	  	  	  	  #Get the second highest allele count and assign the count to $SecondMaxCount and the element position of this count to $SecondMaxPosition.
	   	  	  	  	  my $SecondMaxCount = 0;	
	   	  	  	  	  my $SecondMaxPosition = 0;
	   	  	  	  	  my $CurrentElementCounter2 = 0;
	   	  	  
	   	  	  	  	  foreach my $c (@TempCountArray)  {
	   	  	  	  	  	  if (($c < $MaxCount) && ($c >= $SecondMaxCount)) {
	   	  	  	  	  	  	  $SecondMaxCount = $c;
	   	  	  	  	  	  	  $SecondMaxPosition = $CurrentElementCounter2;
	   	  	  	  	  	  }
	   	  	  	  	  	  $CurrentElementCounter2++;
	   	  	  	  	  }
	   	  	  
	   	  	  	  	  print WRITEFILE "LocusNumber$PreviousLocusNumber\tA\n";
	   	  	  	  	  print WRITEFILE "$TempSeqArray[$MaxPosition]\t$TempCountArray[$MaxPosition]\n";
	   	  	  	  	  print WRITEFILE "$TempSeqArray[$SecondMaxPosition]\t$TempCountArray[$SecondMaxPosition]\n";
	   	  	  	  }
	   	  	  }
	   	  
	   	  	} 
	   	  	  
	   	   }
	   	   close READFILE;
	   	   close WRITEFILE;
	}
}

print "\n\nLoci have been identified.  Aligned polymorphic loci are stored in the Output directory.  Read counts for the alleles at these loci can be found for each individual in the files stored in TempFiles/RawReadCountFiles.\n\n";
print "Demultiplexing and other general information can be found for each dataset analyzed in the Output/RunInfo directory.\n\n";
print "Next, run Genotype.pl, followed by FilterSNPs.pl to genotype all individuals and output SNPs for subsequent analyses.\n\n";

