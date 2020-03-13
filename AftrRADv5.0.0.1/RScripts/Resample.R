DataTableFromFile<-read.table("../out/Output/Genotypes/SNPMatrix.txt", header=TRUE, sep = "\t")
DataMatrix<-as.matrix(DataTableFromFile)
Names<-c(colnames(DataTableFromFile))
DataMatrix<-rbind(Names,DataMatrix)
NumRowsInMatrix<-nrow(DataMatrix)
NumColsInMatrix<-ncol(DataMatrix)

NumLociToSample<-

columns<-c(sample(2:NumColsInMatrix, NumLociToSample, replace=F))

ResampledMatrix<-matrix(DataMatrix[,1])

for (i in columns)  {

     ResampledMatrix<-cbind(ResampledMatrix,DataMatrix[,i])
}


write.table(ResampledMatrix, file="ResampledDatasets/SNPMatrix_resamp.txt", sep="\t", row.names=FALSE, col.names=FALSE)
