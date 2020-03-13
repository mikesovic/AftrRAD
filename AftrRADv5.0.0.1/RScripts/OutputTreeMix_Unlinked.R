#Read in data as data frame


DataFrame<-read.table("out/TempFiles/SingleSNPsAll.txt", header=TRUE, sep = "\t")


NumRowsInFrame<-nrow(DataFrame)		#this doesn't include the locus names
NumColsInFrame<-ncol(DataFrame)		#this does include the sample names


#Get unlinked BiallelicSNPs and store in data frame UnlinkedBiallelicFrame


UnlinkedFrame<-data.frame(DataFrame[,1])


LocusNames<-"starter"

CurrentLocus<-"aa"

CurrentColumn<-1

for (a in 2:(NumColsInFrame))  {	#Check each locus one at a time
	
	
	CurrentLocus<-colnames(DataFrame[a])
	CurrentLocus<-sub("\\.[0-9]+", "", CurrentLocus, perl=TRUE) 
	Match<-0
	
	for (c in LocusNames)  {
		if (CurrentLocus == c)  {
			Match<-1
			break
		}	
	}
	
	if (Match == 0) {
		CurrentColumn<-CurrentColumn+1
		LocusNames<-c(LocusNames,CurrentLocus)
			
		VectorToPrint<-c(as.character(DataFrame[,a]))
		UnlinkedFrame<-cbind(UnlinkedFrame,VectorToPrint)
		colnames(UnlinkedFrame)[CurrentColumn]<-CurrentLocus
		
	}
	
}






#Check each column to see if it is biallelic.  Create new data frame with these biallelics.

NumRowsInFrame<-nrow(UnlinkedFrame)		#this doesn't include the locus names
NumColsInFrame<-ncol(UnlinkedFrame)		#this does include the sample names

BiallelicFrame<-data.frame(UnlinkedFrame[,1])
colnames(BiallelicFrame)[1]<-"SampleNames"

CurrentBiallelicColumn<-1

for (i in 2:NumColsInFrame)  {

	if (length(levels(UnlinkedFrame[,i])) == 2)  {
		CurrentBiallelicColumn<-CurrentBiallelicColumn+1
		BiallelicFrame<-cbind(BiallelicFrame,UnlinkedFrame[,i])
		colnames(BiallelicFrame)[CurrentBiallelicColumn]<-colnames(UnlinkedFrame)[i]
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





#Now have data frame with only unlinked biallelic loci (BiallelicFrame).  Go through one locus/population at a time to count alleles for TreeMix file.

NumPops<-length(PopSizesVector)
UnlinkedTreeMixMatrix<-matrix(0,ncol=NumPops,nrow=ncol(BiallelicFrame)-1,dimnames=list(colnames(BiallelicFrame)[-1], c(1:NumPops)))

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
		UnlinkedTreeMixMatrix[locus-1,popnumber] = ToPrint
	}
}



write.table(UnlinkedTreeMixMatrix, file = "TreeMix_Infile_UnlinkedSNPs.txt", sep = " ",row.names=FALSE, col.names=TRUE, quote=FALSE)



