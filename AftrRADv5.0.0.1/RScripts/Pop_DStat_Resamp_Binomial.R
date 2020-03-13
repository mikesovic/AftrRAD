
library(seqinr)


#############################################################################################################################
# D-statistic calculation based on allele frequencies - requires fasta file with samples grouped by population in order P1, P2, P3, outgroup (P4).
#Prior to running the script, do the following...
# 1.) update the SampleSizesVector	(line 17)
# 2.) update the file name to read in (line 30)
# 3.) set the number of simulated datasets to generate (default is 100).  (line 303)

#make sure the fasta file is in the working directory, and run the script in R.
#When script is finished running, type 'd' to get the observed d statistic, and 'sim.d' to get a set of 100 resampled d-statistics. No/little overlap with zero should be an indication of a significant population d-statistic.
#############################################################################################################################

#Must update samples sizes for 4 populations here, in same order they appear in fasta file (total number of chromosomes sampled)
SampleSizesVector<-c(16,24,26,6)

#get a vector with cumulative sample sizes - will use this to keep track of what population we're on as we work down through the matrix at various points in the script.
CumulativeSizesVector<-c(SampleSizesVector[1])

for (i in 2:length(SampleSizesVector)) {
	CumulativeSizesVector[i]<-CumulativeSizesVector[i-1]+SampleSizesVector[i]
}



###read in the alignment and load it into matrix 'alignment.matrix'.

alignment <- read.alignment(file="kzorresebe.fasta", format = "fasta")                          #  read in the alignment
alignment.matrix <- matrix(, length(alignment$nam), nchar(alignment$seq[[1]]))    #  make a matrix for the alignment
   
for(i in 1:length(alignment$nam)){
    alignment.matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))               #  fill in the matrix
}
  
  
#will update these values on a locus by locus basis and divide final values to get d statistic.   
dstatNumerator<-0
dstatDenominator<-0
	
#Define the number of samples and loci in the dataset
NumSamples<-length(alignment$nam)
NumLoci <- ncol(alignment.matrix)
	
#Keep identity of the ancestral and derived alleles at all biallelic loci
AncestralAlleles<-c()
DerivedAlleles<-c()
NumBiallelicLoci<-0

#Keep track of allele counts at each locus in each population - two vectors per population (ancestral and derived alleles).
pop1counts.a<-c()
pop1counts.d<-c()
pop2counts.a<-c()
pop2counts.d<-c()
pop3counts.a<-c()
pop3counts.d<-c()
pop4counts.a<-c()
pop4counts.d<-c()

TotalAllelesSampledPop1<-c()
TotalAllelesSampledPop2<-c()
TotalAllelesSampledPop3<-c()
TotalAllelesSampledPop4<-c()


#Keep track of observed derived allele freqs
Pop1DerivedFreqList<-c()
Pop2DerivedFreqList<-c()
Pop3DerivedFreqList<-c()
Pop4DerivedFreqList<-c()



	        
for (i in 1:NumLoci) { 	
	
	#first, check to see if the locus is biallelic
	num.matches<-0
	current.bases<-c()
	TempAlleles<-c(alignment.matrix[,i])
	
	for (b in (c("a","t","c","g"))) {
		if (b %in% TempAlleles) {
			current.bases<-c(current.bases,b)
			num.matches<-num.matches+1
		}
	}
		
	#if the locus is not biallelic, move to the next locus		
	if (num.matches != 2) {
	  next
	}
	
	NumBiallelicLoci<-NumBiallelicLoci+1
	
	#if the locus is biallelic, define ancestral and derived alleles. 
	AncestralAllele<-c()
	DerivedAllele<-c()
	Outgroup<-c(alignment.matrix[(CumulativeSizesVector[3]+1):CumulativeSizesVector[4],i])
	
	#Check to make sure the outgroup has at least one base scored (not all n's or -'s).  If not, define the ancestral allele arbitrarily.
	TF.vector<-c("a","t","c","g") %in% Outgroup
	if ("TRUE" %in% TF.vector) {
		OutGroupNumAllele1<-length(Outgroup[Outgroup==current.bases[1]])
		OutGroupNumAllele2<-length(Outgroup[Outgroup==current.bases[2]])	
	
		if (OutGroupNumAllele1 > OutGroupNumAllele2) {
			AncestralAllele<-current.bases[1]
			AncestralAlleles<-c(AncestralAlleles,current.bases[1])
			DerivedAllele<-current.bases[2]
			DerivedAlleles<-c(DerivedAlleles,current.bases[2])
		}
		
		else {
			AncestralAllele<-current.bases[2]
			AncestralAlleles<-c(AncestralAlleles,current.bases[2])
			DerivedAllele<-current.bases[1]
			DerivedAlleles<-c(DerivedAlleles,current.bases[1])	
		}
	}
	
	#if no information from outgroup, then define ancestral allele arbitrarily.
	
	else {
	 	AncestralAllele<-current.bases[1]
	 	AncestralAlleles<-c(AncestralAlleles,current.bases[1])
	 	DerivedAllele<-current.bases[2]
	 	DerivedAlleles<-c(DerivedAlleles,current.bases[2])	
	}
	
	
	#calculate the numerator and denominator for the d statistic calculation.
	
	AncestralAllelePop1Count<-0
	DerivedAllelePop1Count<-0
	AncestralAllelePop2Count<-0
	DerivedAllelePop2Count<-0
	AncestralAllelePop3Count<-0
	DerivedAllelePop3Count<-0
	AncestralAllelePop4Count<-0
	DerivedAllelePop4Count<-0
	NumAllelesSampledAtLocusPop1<-0
	NumAllelesSampledAtLocusPop2<-0
	NumAllelesSampledAtLocusPop3<-0
	NumAllelesSampledAtLocusPop4<-0
	
	for (j in 1:NumSamples) {			
		if (alignment.matrix[j,i] == AncestralAllele) {
					
			if (j <= CumulativeSizesVector[1]) {
				AncestralAllelePop1Count<-AncestralAllelePop1Count+1
				NumAllelesSampledAtLocusPop1<-NumAllelesSampledAtLocusPop1+1
			}
				
			else if (j <= CumulativeSizesVector[2]) {
				AncestralAllelePop2Count<-AncestralAllelePop2Count+1
				NumAllelesSampledAtLocusPop2<-NumAllelesSampledAtLocusPop2+1
			}
				
			else if (j <= CumulativeSizesVector[3]) {
				AncestralAllelePop3Count<-AncestralAllelePop3Count+1
				NumAllelesSampledAtLocusPop3<-NumAllelesSampledAtLocusPop3+1
			}
				
			else {
				AncestralAllelePop4Count<-AncestralAllelePop4Count+1
				NumAllelesSampledAtLocusPop4<-NumAllelesSampledAtLocusPop4+1
			}	
				
			
		}
						
		else if (alignment.matrix[j,i] == DerivedAllele) {
			
			if (j <= CumulativeSizesVector[1]) {
				DerivedAllelePop1Count<-DerivedAllelePop1Count+1
				NumAllelesSampledAtLocusPop1<-NumAllelesSampledAtLocusPop1+1
			}
				
			else if (j <= CumulativeSizesVector[2]) {
				DerivedAllelePop2Count<-DerivedAllelePop2Count+1
				NumAllelesSampledAtLocusPop2<-NumAllelesSampledAtLocusPop2+1
			}
				
			else if (j <= CumulativeSizesVector[3]) {
				DerivedAllelePop3Count<-DerivedAllelePop3Count+1
				NumAllelesSampledAtLocusPop3<-NumAllelesSampledAtLocusPop3+1
			}
				
			else {
				DerivedAllelePop4Count<-DerivedAllelePop4Count+1
				NumAllelesSampledAtLocusPop4<-NumAllelesSampledAtLocusPop4+1
			}

		}
	}
	
	
	
	
	Pop1AncestralFreq<-AncestralAllelePop1Count/NumAllelesSampledAtLocusPop1
	Pop1DerivedFreq<-DerivedAllelePop1Count/NumAllelesSampledAtLocusPop1
	Pop2AncestralFreq<-AncestralAllelePop2Count/NumAllelesSampledAtLocusPop2
	Pop2DerivedFreq<-DerivedAllelePop2Count/NumAllelesSampledAtLocusPop2
	Pop3AncestralFreq<-AncestralAllelePop3Count/NumAllelesSampledAtLocusPop3
	Pop3DerivedFreq<-DerivedAllelePop3Count/NumAllelesSampledAtLocusPop3
	Pop4AncestralFreq<-AncestralAllelePop4Count/NumAllelesSampledAtLocusPop4
	Pop4DerivedFreq<-DerivedAllelePop4Count/NumAllelesSampledAtLocusPop4
	
	#In cases of missing data, it's possible to not have samples from populations for some loci, causing division by zero above.  
	#In these cases, set the ancestral and derived frequencies to zero.
	
	if (NumAllelesSampledAtLocusPop1 == 0)  {
		Pop1AncestralFreq<-0
		Pop1DerivedFreq<-0
	}
	
	if (NumAllelesSampledAtLocusPop2 == 0)  {
		Pop2AncestralFreq<-0
		Pop2DerivedFreq<-0
	}
	
	if (NumAllelesSampledAtLocusPop3 == 0)  {
		Pop3AncestralFreq<-0
		Pop3DerivedFreq<-0
	}
	
	if (NumAllelesSampledAtLocusPop4 == 0)  {
		Pop4AncestralFreq<-0
		Pop4DerivedFreq<-0
	}
			
	
	#Update the running lists with info from the current locus
	
	Pop1DerivedFreqList<-c(Pop1DerivedFreqList,Pop1DerivedFreq)
	Pop2DerivedFreqList<-c(Pop2DerivedFreqList,Pop2DerivedFreq)
	Pop3DerivedFreqList<-c(Pop3DerivedFreqList,Pop3DerivedFreq)
	Pop4DerivedFreqList<-c(Pop4DerivedFreqList,Pop4DerivedFreq)
	
	TotalAllelesSampledPop1<-c(TotalAllelesSampledPop1,NumAllelesSampledAtLocusPop1)
	TotalAllelesSampledPop2<-c(TotalAllelesSampledPop2,NumAllelesSampledAtLocusPop2)
	TotalAllelesSampledPop3<-c(TotalAllelesSampledPop3,NumAllelesSampledAtLocusPop3)
	TotalAllelesSampledPop4<-c(TotalAllelesSampledPop4,NumAllelesSampledAtLocusPop4)
	
	#set 'p' for calculation of Equation 2 in Durand et al (2011) as the derived allele (or SNP) frequency, and calculate the numerator and denominator for the current locus.
		
	TempNumeratorValue<-((1-Pop1DerivedFreq)*Pop2DerivedFreq*Pop3DerivedFreq*(1-Pop4DerivedFreq))-(Pop1DerivedFreq*(1-Pop2DerivedFreq)*Pop3DerivedFreq*(1-Pop4DerivedFreq))
	TempDenominatorValue<-((1-Pop1DerivedFreq)*Pop2DerivedFreq*Pop3DerivedFreq*(1-Pop4DerivedFreq))+(Pop1DerivedFreq*(1-Pop2DerivedFreq)*Pop3DerivedFreq*(1-Pop4DerivedFreq))
			
	dstatNumerator<-dstatNumerator+TempNumeratorValue
	dstatDenominator<-dstatDenominator+TempDenominatorValue
	#print(c(dstatNumerator,dstatDenominator))
	
	pop1counts.a<-c(pop1counts.a,AncestralAllelePop1Count)
	pop1counts.d<-c(pop1counts.d,DerivedAllelePop1Count)
	pop2counts.a<-c(pop2counts.a,AncestralAllelePop2Count)
	pop2counts.d<-c(pop2counts.d,DerivedAllelePop2Count)
	pop3counts.a<-c(pop3counts.a,AncestralAllelePop3Count)
	pop3counts.d<-c(pop3counts.d,DerivedAllelePop3Count)
	pop4counts.a<-c(pop4counts.a,AncestralAllelePop4Count)
	pop4counts.d<-c(pop4counts.d,DerivedAllelePop4Count)
	
}
	
d<-dstatNumerator/dstatDenominator	
	
	
print(d)	
print(d)	
print(d)	
print(d)	
	
		


 
##########################################################################################
##########################################################################################
#Option for calculating confidence in D statistic by resampling from beta distribution at each locus at each population, as suggested by P. Blischak. 
 
#Generate k simulated datasets by sampling from a beta distribution based on the observed allele counts at each locus in each population. 		
#For each, calculate d and store this value in the vector 'sim.d'

sim.d<-vector()

for (k in 1:100)  {
	print(k)
	
	dstatNumerator<-0
	dstatDenominator<-0
	
	for (i in 1:NumBiallelicLoci) {
		#Sample from the beta distribution for each population for the current locus.
		#print(c("Locus",i))
		Pop1Count.d<-rbinom(1,TotalAllelesSampledPop1[i],Pop1DerivedFreqList[i])
		Pop2Count.d<-rbinom(1,TotalAllelesSampledPop2[i],Pop2DerivedFreqList[i])
		Pop3Count.d<-rbinom(1,TotalAllelesSampledPop3[i],Pop3DerivedFreqList[i])
		Pop4Count.d<-rbinom(1,TotalAllelesSampledPop4[i],Pop4DerivedFreqList[i])
		
		Pop1Proportion.d<-Pop1Count.d/TotalAllelesSampledPop1[i]
		Pop2Proportion.d<-Pop2Count.d/TotalAllelesSampledPop2[i]
		Pop3Proportion.d<-Pop3Count.d/TotalAllelesSampledPop3[i]
		Pop4Proportion.d<-Pop4Count.d/TotalAllelesSampledPop4[i]
		
		if (TotalAllelesSampledPop1[i] == 0)  {
			Pop1Proportion.d<-0
		}
		
		if (TotalAllelesSampledPop2[i] == 0)  {
			Pop2Proportion.d<-0
		}
		
		if (TotalAllelesSampledPop3[i] == 0)  {
			Pop3Proportion.d<-0
		}
		
		if (TotalAllelesSampledPop4[i] == 0)  {
			Pop4Proportion.d<-0
		}
	
		
		
		#print(c(Pop1Proportion.d,Pop2Proportion.d,Pop3Proportion.d,Pop4Proportion.d))
		#print(c(Pop1DerivedFreqList[i],Pop2DerivedFreqList[i],Pop3DerivedFreqList[i],Pop4DerivedFreqList[i]))
	
		#Calculate d for the current locus
	
		TempNumeratorValue<-((1-Pop1Proportion.d)*Pop2Proportion.d*Pop3Proportion.d*(1-Pop4Proportion.d))-(Pop1Proportion.d*(1-Pop2Proportion.d)*Pop3Proportion.d*(1-Pop4Proportion.d))
		TempDenominatorValue<-((1-Pop1Proportion.d)*Pop2Proportion.d*Pop3Proportion.d*(1-Pop4Proportion.d))+(Pop1Proportion.d*(1-Pop2Proportion.d)*Pop3Proportion.d*(1-Pop4Proportion.d))
			
		dstatNumerator<-dstatNumerator+TempNumeratorValue
		dstatDenominator<-dstatDenominator+TempDenominatorValue
		#print(c(dstatNumerator,dstatDenominator))
	
	}
	
	sim.d[k]<-dstatNumerator/dstatDenominator
	
}



sd.sim.d <- round(sqrt(var(sim.d)),5)
mn.sim.d <- round(mean(sim.d),5)
new.pval <- 2*(pnorm(-abs(d/sd.sim.d)))

NumBootstrapsFromPerl<-
string<-c(d,NumBootstrapsFromPerl,mn.sim.d,sd.sim.d,d/sd.sim.d,new.pval)

ResultsMatrix<-matrix(string,nrow=1,ncol=6,byrow=TRUE)
colnames(ResultsMatrix)<-c("d","num.resamples","mean.sim.d","sd.sim.d","d/sd.sim.d","p_value")

write.table(ResultsMatrix, file="DStat_pop_out.txt", sep="\t", row.names=FALSE, col.names=TRUE)























