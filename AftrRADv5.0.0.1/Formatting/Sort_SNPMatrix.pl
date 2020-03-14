#! usr/bin/perl
use strict;
use warnings;

mkdir "../out/formatted_files" unless(-d "../out/formatted_files");

my $FileName;

opendir GENOS, "../out/Output/Genotypes";
my @AllFiles = grep { $_ ne '.' && $_ ne '..' && $_ ne '.DS_Store' } readdir(GENOS);
close GENOS;


my @HaplotypeFileNames = ();

for my $name (@AllFiles)  {
	if ($name =~ /SNPMat/) {
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
	print "\nEnter the name of the SNPMatrix file you want to sort\n";

	$FileName = <STDIN>;
	chomp($FileName);
}



open HAPS, "../out/Output/Genotypes/$FileName" or die$!;

my $Counter = 0;

my %UnsortedHash = ();

my $LocusNames;

while(<HAPS>) {
	if ($Counter == 0) {
		$LocusNames = $_;
		$Counter++;
		next;
	}

	else {
		if ($_ =~ /[a-zA-Z0-9]/) {
			my @TempArray = split(/\t/, $_);
			my $TempName = $TempArray[0];
			shift(@TempArray);
			
			my $SNPs = join("\t", @TempArray);
			chomp($SNPs);
			
			if ($UnsortedHash{$TempName}) {
				$TempName = $TempName."_BCB23Zd";
			}	
			
			$UnsortedHash{$TempName} = $SNPs;
		}
	}
}

close HAPS;

my @SortedArray = sort {lc($a) cmp lc($b)} keys %UnsortedHash;

$FileName =~ s/.txt$//;
my $NewFileName = $FileName."_Sorted";

open OUT, ">../out/Output/Genotypes/$NewFileName.txt" or die$!;

print OUT "$LocusNames";

for my $name (@SortedArray) {
	my $SNPsToPrint = $UnsortedHash{$name};
	$name =~ s/_BCB23Zd$//;
	print OUT "$name\t";
	print OUT "$SNPsToPrint\n";
}

close OUT;
