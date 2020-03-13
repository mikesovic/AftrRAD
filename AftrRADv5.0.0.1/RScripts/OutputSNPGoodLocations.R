#This script goes through and removes SNPs that are in positions along the read that are higher than the max position allowed (input by user)
#This is to remove SNPs toward the ends of reads that are of lower quality.

LocationsToKeep<-scan(file="out/TempFiles/SitesToRetain.txt")

DataTableFromFile<-read.table("out/TempFiles/TempSNPMatrix.txt", header=TRUE, sep = "\t")  #####Have to edit this line.

DataMatrix<-as.matrix(DataTableFromFile)
DataFrame<-data.frame(DataMatrix)
NumRowsInMatrix<-nrow(DataMatrix)
NumColsInMatrix<-ncol(DataMatrix)




RetainedSNPsMatrix<-matrix(c("",DataMatrix[1:NumRowsInMatrix,1]),ncol=1,nrow=NumRowsInMatrix+1)

for (i in LocationsToKeep) {
	CurrentLocation<-i+1
	CurrentLocusName<-colnames(DataFrame[CurrentLocation])
	RetainedSNPsMatrix<-cbind(RetainedSNPsMatrix,c(CurrentLocusName,DataMatrix[,CurrentLocation]))
}

ncolInMatrix<-ncol(RetainedSNPsMatrix)

MatrixToPrint<-RetainedSNPsMatrix[,1:ncolInMatrix]

write.table(MatrixToPrint,file="out/TempFiles/SNPMatrix_GoodLocations.txt",sep="\t",row.names=FALSE, col.names=FALSE)