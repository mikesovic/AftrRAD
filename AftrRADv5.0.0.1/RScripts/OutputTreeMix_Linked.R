
#Read in data as data frame


DataFrame<-read.table("out/TempFiles/SingleSNPsAll.txt", header=TRUE, sep = "\t")

NumRowsInFrame<-nrow(DataFrame)		#this doesn't include the locus names
NumColsInFrame<-ncol(DataFrame)		#this does include the sample names



#Check each column to see if it is biallelic.  Create new data frame with these biallelics.

BiallelicFrame<-data.frame(DataFrame[,1])
colnames(BiallelicFrame)[1]<-"SampleNames"

CurrentBiallelicColumn<-1

for (i in 2:NumColsInFrame)  {

	if (length(levels(DataFrame[,i])) == 2)  {
		CurrentBiallelicColumn<-CurrentBiallelicColumn+1
		BiallelicFrame<-cbind(BiallelicFrame,DataFrame[,i])
		colnames(BiallelicFrame)[CurrentBiallelicColumn]<-colnames(DataFrame)[i]
	}
}	


#Create vector with the number of alleles sampled in each population, in order - this comes from perl script.
PopSizesVector<-c(12,6,6)

#Create the cumulative PopSizesVector

CumulativePopSizesVector<-c(0)

for (i in 1:length(PopSizesVector))   {
	LastCumulativeValue<-tail(CumulativePopSizesVector, n=1)
	CumulativePopSizesVector<-c(CumulativePopSizesVector,LastCumulativeValue+PopSizesVector[i])
}



#Now have data frame with only biallelic loci (BiallelicFrame).  Go through one locus/population at a time to count alleles for TreeMix file.

NumPops<-length(PopSizesVector)
#TreeMixMatrix<-matrix(0,ncol=NumPops,nrow=ncol(BiallelicFrame)-1,dimnames=list(c(1:(ncol(BiallelicFrame)-1))), c(1:NumPops))
TreeMixMatrix<-matrix(0,ncol=NumPops,nrow=ncol(BiallelicFrame)-1,dimnames=list(colnames(BiallelicFrame)[-1], c(1:NumPops)))

NumBiallelicColumns<-ncol(BiallelicFrame)
	
for (locus in 2:(NumBiallelicColumns)) {
	
	RefAllele=NULL
	
	#Go through samples at the current locus to find the first one that has a base call.
	for (sample in 1:CumulativePopSizesVector[(length(PopSizesVector)+1)]) {
		if ((is.na(BiallelicFrame[sample,locus])) || (BiallelicFrame[sample,locus] == "-")) {
		 	next
		}
		
		else {
			 RefAllele<-BiallelicFrame[sample,locus]	
		 	 break
		} 	 
	}
	
	for (popnumber in 1:length(PopSizesVector)) {
	
		
		RefCount<-0
			
		NonRefCount<-0
	
		for (row in (CumulativePopSizesVector[popnumber]+1):(CumulativePopSizesVector[popnumber+1]))  {
		
			if (is.na(BiallelicFrame[row,locus])) {
			   next
			}
			
			else {   
				if ((BiallelicFrame[row,locus] == RefAllele)) {
					RefCount = RefCount+1
				}
			
				else {
					
					NonRefCount = NonRefCount+1
			
					
				}
			}
		}
		
		ToPrint<-paste(RefCount,NonRefCount, sep=",")
		TreeMixMatrix[locus-1,popnumber] = ToPrint
	}
}



write.table(TreeMixMatrix, file = "TreeMix_Infile_AllSNPs.txt", sep = " ",row.names=FALSE, col.names=TRUE, quote=FALSE)





















