DataTableFromFile<-read.table("out/TempFiles/SingleSNPsAll.txt", header=TRUE, sep = "\t")

DataMatrix<-as.matrix(DataTableFromFile)
DataFrame<-data.frame(DataMatrix)
NumRowsInMatrix<-nrow(DataMatrix)
NumColsInMatrix<-ncol(DataMatrix)

NamesCol<-DataMatrix[,1]





UnlinkedMatrix<-matrix(c("",NamesCol), ncol=1, nrow=NumRowsInMatrix+1)

LocusNames<-"starter"

CurrentLocus<-"aa"

for (a in 2:(NumColsInMatrix))  {	#Check each locus one at a time
	
	
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
		LocusNames<-c(LocusNames,CurrentLocus)
			
		VectorToPrint<-c(CurrentLocus, as.character(DataFrame[,a]))
		UnlinkedMatrix<-cbind(UnlinkedMatrix,VectorToPrint)
	}
	
}

write.table(UnlinkedMatrix, file = "out/TempFiles/UnlinkedSNPsRaw.txt", sep = "\t",row.names=FALSE, col.names=FALSE)



