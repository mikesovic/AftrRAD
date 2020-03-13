DataTableFromFile<-read.table("out/TempFiles/FinalBiallelicSNPs.txt", header=TRUE, sep = "\t")

DataMatrix<-as.matrix(DataTableFromFile)
DataFrame<-data.frame(DataMatrix, stringsAsFactors=FALSE)
NumRowsInMatrix<-nrow(DataMatrix)
NumColsInMatrix<-ncol(DataMatrix)

NamesCol<-DataMatrix[,1]

#Create empty matrix and put sample names as first column.
SNAPPMatrix<-matrix(c("",NamesCol), ncol=1, nrow=NumRowsInMatrix+1)

#For each locus (column)

for (a in 2:(NumColsInMatrix))  {

	HaveFirstAllele<-0
	CurrentLocus<-colnames(DataFrame[a])
	
	for (b in 1:NumRowsInMatrix)  {
		
		if (as.character(DataFrame[b,a]) == "N") {
			DataFrame[b,a]<-"N"
			next
		}
		
		else if (HaveFirstAllele == 0)  {
			FirstAllele = DataFrame[b,a]
			HaveFirstAllele<-1
			DataFrame[b,a]<-"0"
		}
		
		else {
			
			if (DataFrame[b,a] == FirstAllele)  {
				DataFrame[b,a]<-"0"
			}
			
			else {
				DataFrame[b,a]<-"1"
			}
		}
	}
	
	CurrentLocus<-colnames(DataFrame[a])
	
	CurrentLocusColToPrint<-c(CurrentLocus, DataFrame[,a])
	
	SNAPPMatrix<-cbind (SNAPPMatrix, CurrentLocusColToPrint)

	
}	

write.table(SNAPPMatrix, file = "out/TempFiles/SNAPPMatrixRaw.txt", sep = "\t",row.names=FALSE, col.names=FALSE)
