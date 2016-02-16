
library(seqinr)



#############################################################################################################################
# D-statistic calculation based on allele frequencies - requires fasta file with samples grouped by population (two sequences per individual) in order P1, P2, P3, outgroup (P4).
#Reads the alignments into matrix 'alignment.matrix', which has dimensions [# chromosomes sampled, # of sites]
#Checks each site to see if 1.)site is biallelic, 2.) P1 and P2 have different alleles, 3.) P3 and P4 have different alleles
#If all of the above are true, then site qualifies for ABBA/BABA test and is scored as one of these.
#############################################################################################################################

###get a vector with cumulative sample sizes

SampleSizesVector<-c(26,26,8,2)

CumulativeSizesVector<-c(SampleSizesVector[1])

for (i in 2:length(SampleSizesVector)) {
	CumulativeSizesVector[i]<-CumulativeSizesVector[i-1]+SampleSizesVector[i]
}



###read in the alignment	

alignment <- read.alignment(file="SNPMatrix_100.32.fasta", format = "fasta")                          #  read in the alignment
alignment.matrix <- matrix(, length(alignment$nam), nchar(alignment$seq[[1]]))    #  make a matrix for the alignment
   
for(i in 1:length(alignment$nam)){
    alignment.matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))               #  fill in the matrix
}
   
dstatNumerator<-0
dstatDenominator<-0
	        
for (i in 1:nchar(alignment$seq[[1]])) {   		#for each locus
	allele1<-character()
	
	allele1totalcount<-0
	allele1pop1count<-0
	allele1pop2count<-0
	allele1pop3count<-0
	allele1pop4count<-0
	
	#keep track of sample sizes that don't include "n's"
	pop1samples<-0
	pop2samples<-0
	pop3samples<-0
	pop4samples<-0
	
	
	allele2totalcount<-0
	allele2pop1count<-0
	allele2pop2count<-0
	allele2pop3count<-0
	allele2pop4count<-0
	
	#get list of unique bases at the current site
	uniques<-unique(alignment.matrix[,i])
	
	#if the current site contains "n's", then looking for 3 total sites, otherwise, 2 (selecting biallelic sites)
	if ((('n' %in% uniques && (length(unique(alignment.matrix[,i]))==3)) || (!('n' %in% uniques) && (length(unique(alignment.matrix[,i]))==2))) && (!('-' %in% uniques)))  {
	
		#define allele 1 - can be anything but "n" or "-".
		
		for (j in 1:length(alignment$nam)) {
			if ((alignment.matrix[j,i] == 't') || (alignment.matrix[j,i] == 'c') || (alignment.matrix[j,i] == 'g') || (alignment.matrix[j,i] == 'a')) {
				allele1<-alignment.matrix[j,i]
				break
			}	
		}
		
		#count and store occurrences of allele1 at the current site for each population
		for (j in 1:length(alignment$nam)) {
			
			
			if (alignment.matrix[j,i] == allele1) {
				
				allele1totalcount<-allele1totalcount+1
				
				if (j <= CumulativeSizesVector[1]) {
					allele1pop1count<-allele1pop1count+1
					pop1samples<-pop1samples+1
				}
				
				else if (j <= CumulativeSizesVector[2]) {
					allele1pop2count<-allele1pop2count+1
					pop2samples<-pop2samples+1
				}
				
				else if (j <= CumulativeSizesVector[3]) {
					allele1pop3count<-allele1pop3count+1
					pop3samples<-pop3samples+1
				}
				
				else {
					allele1pop4count<-allele1pop4count+1
					pop4samples<-pop4samples+1
				}	
				
			
			}
			
			else {			#check to see if it is an 'n' and if not, count it as a sample
			
				if (j <= CumulativeSizesVector[1] && (alignment.matrix[j,i] != 'n')) {
					pop1samples<-pop1samples+1
				}
				
				else if (j <= CumulativeSizesVector[2] && (alignment.matrix[j,i] != 'n')) {
					pop2samples<-pop2samples+1
				}
				
				else if (j <= CumulativeSizesVector[3] && (alignment.matrix[j,i] != 'n')) {
					pop3samples<-pop3samples+1
				}
				
				else if (j <= CumulativeSizesVector[4] && (alignment.matrix[j,i] != 'n')) {
					pop4samples<-pop4samples+1
				}
			
			}
		
		}
	
		allele2totalcount<-CumulativeSizesVector[4]-allele1totalcount
		
		#set p as the least frequent allele in P4?
		
		if ((allele1pop4count/pop4samples) <= 0.5)  {			#going to use the pop1counts
			TempNumeratorValue<-((1-(allele1pop1count/pop1samples))*(allele1pop2count/pop2samples)*(allele1pop3count/pop3samples)*(1-(allele1pop4count/pop4samples)))-((allele1pop1count/pop1samples)*(1-(allele1pop2count/pop2samples))*(allele1pop3count/pop3samples)*(1-(allele1pop4count/pop4samples)))
			TempDenominatorValue<-((1-(allele1pop1count/pop1samples))*(allele1pop2count/pop2samples)*(allele1pop3count/pop3samples)*(1-(allele1pop4count/pop4samples)))+((allele1pop1count/pop1samples)*(1-(allele1pop2count/pop2samples))*(allele1pop3count/pop3samples)*(1-(allele1pop4count/pop4samples)))
			
			dstatNumerator<-dstatNumerator+TempNumeratorValue
			dstatDenominator<-dstatDenominator+TempDenominatorValue
	
		}
		
		else {
			
			allele2pop1count<-pop1samples-allele1pop1count
			allele2pop2count<-pop2samples-allele1pop2count
			allele2pop3count<-pop3samples-allele1pop3count
			allele2pop4count<-pop4samples-allele1pop4count
			
			TempNumeratorValue<-((1-(allele2pop1count/pop1samples))*(allele2pop2count/pop2samples)*(allele2pop3count/pop3samples)*(1-(allele2pop4count/pop4samples)))-((allele2pop1count/pop1samples)*(1-(allele2pop2count/pop2samples))*(allele2pop3count/pop3samples)*(1-(allele2pop4count/pop4samples)))
			TempDenominatorValue<-((1-(allele2pop1count/pop1samples))*(allele2pop2count/pop2samples)*(allele2pop3count/pop3samples)*(1-(allele2pop4count/pop4samples)))+((allele2pop1count/pop1samples)*(1-(allele2pop2count/pop2samples))*(allele2pop3count/pop3samples)*(1-(allele2pop4count/pop4samples)))
			
			dstatNumerator<-dstatNumerator+TempNumeratorValue
			dstatDenominator<-dstatDenominator+TempDenominatorValue
	
		
		}
			
	}
	
	
	else {
		next
	} 

}   
 
d<-dstatNumerator/dstatDenominator 
 
 
 
 ###Now do a bootstrap
 
 
 
 
sim.d<-vector()
foo <- ncol(alignment.matrix)
sim.matrix<-matrix(,length(alignment$nam),foo)
     
for(k in 1:1000){    															#doing 1000 bootstrap replicates, each time add the result to vector sim.d. 
       
    for(j in 1:length(alignment$nam)){	  															#there are 4 rows in the matrix (4 individuals in the test).
         #sim.matrix[j,1:foo] <-sample(alignment.matrix[j,1:foo],replace=T)
    	sim.matrix[j,1:foo] <-sample(c("a","t"),foo,replace=T)
    }


	dstatNumerator<-0
	dstatDenominator<-0
	        
	for (i in 1:nchar(alignment$seq[[1]])) {   		#for each locus
		allele1<-character()
	
		allele1totalcount<-0
		allele1pop1count<-0
		allele1pop2count<-0
		allele1pop3count<-0
		allele1pop4count<-0
	
		#keep track of sample sizes that don't include "n's"
		pop1samples<-0
		pop2samples<-0
		pop3samples<-0
		pop4samples<-0
	
	
		allele2totalcount<-0
		allele2pop1count<-0
		allele2pop2count<-0
		allele2pop3count<-0
		allele2pop4count<-0
	
		#get list of unique bases at the current site
		uniques<-unique(sim.matrix[,i])
	
		#if the current site contains "n's", then looking for 3 total sites, otherwise, 2 (selecting biallelic sites)
		if ((('n' %in% uniques && (length(unique(sim.matrix[,i]))==3)) || (!('n' %in% uniques) && (length(unique(sim.matrix[,i]))==2))) && (!('-' %in% uniques)))  {
			#define allele 1 - can be anything but "n"
		
			for (j in 1:length(alignment$nam)) {
				if ((sim.matrix[j,i] == 't') || (sim.matrix[j,i] == 'c') || (sim.matrix[j,i] == 'g') || (sim.matrix[j,i] == 'a')) {
					allele1<-sim.matrix[j,i]
					break
				}	
			}
		
			#count and store occurrences of allele1 at the current site for each population
			for (j in 1:length(alignment$nam)) {
			
			
				if (sim.matrix[j,i] == allele1) {
				
					allele1totalcount<-allele1totalcount+1
				
					if (j <= CumulativeSizesVector[1]) {
						allele1pop1count<-allele1pop1count+1
						pop1samples<-pop1samples+1
					}
				
					else if (j <= CumulativeSizesVector[2]) {
						allele1pop2count<-allele1pop2count+1
						pop2samples<-pop2samples+1
					}
				
					else if (j <= CumulativeSizesVector[3]) {
						allele1pop3count<-allele1pop3count+1
						pop3samples<-pop3samples+1
					}
				
					else {
						allele1pop4count<-allele1pop4count+1
						pop4samples<-pop4samples+1
					}	
				
			
				}
			
				else {			#check to see if it is an 'n' and if not, count it as a sample
			
					if (j <= CumulativeSizesVector[1] && (sim.matrix[j,i] != 'n')) {
						pop1samples<-pop1samples+1
					}
				
					else if (j <= CumulativeSizesVector[2] && (sim.matrix[j,i] != 'n')) {
						pop2samples<-pop2samples+1
					}
				
					else if (j <= CumulativeSizesVector[3] && (sim.matrix[j,i] != 'n')) {
						pop3samples<-pop3samples+1
					}
				
					else if (j <= CumulativeSizesVector[4] && (sim.matrix[j,i] != 'n')) {
						pop4samples<-pop4samples+1
					}
			
				}
		
			}
	
			allele2totalcount<-CumulativeSizesVector[4]-allele1totalcount
		
			#set p as the least frequent allele in P4?
		
			if ((allele1pop4count/pop4samples) <= 0.5)  {			#going to use the pop1counts
				TempNumeratorValue<-((1-(allele1pop1count/pop1samples))*(allele1pop2count/pop2samples)*(allele1pop3count/pop3samples)*(1-(allele1pop4count/pop4samples)))-((allele1pop1count/pop1samples)*(1-(allele1pop2count/pop2samples))*(allele1pop3count/pop3samples)*(1-(allele1pop4count/pop4samples)))
				TempDenominatorValue<-((1-(allele1pop1count/pop1samples))*(allele1pop2count/pop2samples)*(allele1pop3count/pop3samples)*(1-(allele1pop4count/pop4samples)))+((allele1pop1count/pop1samples)*(1-(allele1pop2count/pop2samples))*(allele1pop3count/pop3samples)*(1-(allele1pop4count/pop4samples)))
			
				dstatNumerator<-dstatNumerator+TempNumeratorValue
				dstatDenominator<-dstatDenominator+TempDenominatorValue
	
			}
		
			else {
			
				allele2pop1count<-pop1samples-allele1pop1count
				allele2pop2count<-pop2samples-allele1pop2count
				allele2pop3count<-pop3samples-allele1pop3count
				allele2pop4count<-pop4samples-allele1pop4count
			
				TempNumeratorValue<-((1-(allele2pop1count/pop1samples))*(allele2pop2count/pop2samples)*(allele2pop3count/pop3samples)*(1-(allele2pop4count/pop4samples)))-((allele2pop1count/pop1samples)*(1-(allele2pop2count/pop2samples))*(allele2pop3count/pop3samples)*(1-(allele2pop4count/pop4samples)))
				TempDenominatorValue<-((1-(allele2pop1count/pop1samples))*(allele2pop2count/pop2samples)*(allele2pop3count/pop3samples)*(1-(allele2pop4count/pop4samples)))+((allele2pop1count/pop1samples)*(1-(allele2pop2count/pop2samples))*(allele2pop3count/pop3samples)*(1-(allele2pop4count/pop4samples)))
			
				dstatNumerator<-dstatNumerator+TempNumeratorValue
				dstatDenominator<-dstatDenominator+TempDenominatorValue
	
		
			}
			
		}
		
	
		else {
			next
		}


	}

	sim.d[k] <- dstatNumerator/dstatDenominator

}

sd.sim.d <- round(sqrt(var(sim.d)),5)
mn.sim.d <- round(mean(sim.d),5)
new.pval <- 2*(pnorm(-abs(d/sd.sim.d)))



 

   ## NOW WE MAKE THE OUTPUTS  
cat("\nSites in alignment =", ncol(alignment.matrix))

 #  cat("\nNumber of sites with ABBA pattern =", abba)

  # cat("\nNumber of sites with BABA pattern =", baba)

cat("\n\nD raw statistic / Z-score = ", d, " / ", d/sd.sim.d)

string<-c(d,mn.sim.d,sd.sim.d,d/sd.sim.d,new.pval)

cat("\n\nResults from ", "NumBootstrapsFromPerl", "bootstraps")

ResultsMatrix<-matrix(string,nrow=1,ncol=5,byrow=TRUE)
colnames(ResultsMatrix)<-c("d","mean.sim.d","sd.sim.d","d/sd.sim.d","p_value")

write.table(ResultsMatrix, file="ABBA_BABA_pop_Out.txt", sep="\t", row.names=FALSE, col.names=TRUE)
