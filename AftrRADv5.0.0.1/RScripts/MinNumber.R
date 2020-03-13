#We have genotyped all individuals at all loci in file GenotypesUpdate, however, this contains the entire read
#Here, we pull out the SNPs.  Some of the loci in GenotypesUpdate may actually be monomorphic, as they appeared polymorphic based on error, which was eliminated by the binomial test.  These will be eliminated here

DataTableFromFile<-read.table("out/TempFiles/SNPMatrix_GoodLocations.txt", header=TRUE, sep = "\t")  
SNPMatrix<-as.matrix(DataTableFromFile)
DataFrame<-data.frame(SNPMatrix)
NumRowsInMatrix<-nrow(SNPMatrix)
NumColsInMatrix<-ncol(SNPMatrix)

Subset<-subset(DataFrame[,1:(NumColsInMatrix-1)])
attach(Subset)

######################################################################################
#FirstColToSplit<-as.character(Subset[,1])
#FirstColSplit<-strsplit(FirstColToSplit,split="\t")
#SplitColDelisted<-unlist(FirstColSplit)
#SNPMatrix<-matrix(c("",SplitColDelisted),ncol=1,nrow=(NumRowsInMatrix+1))


#NumRowsInSNPMatrix<-nrow(SNPMatrix)
#NumColsInSNPMatrix<-ncol(SNPMatrix)




#Need to make this an option...
MinNumberAllelesGenotypedToRetainLocus<-NumRowsInMatrix*1
SNPMatrixWithLociWithMinGenotypes<-matrix(c("", SNPMatrix[,1]),ncol=1, nrow=(NumRowsInMatrix+1))

for (f in 2:NumColsInMatrix)  {
   
    Counter=0
    
    for (g in 1:NumRowsInMatrix)  {
       if (!is.na(SNPMatrix[g,f]))  {
         Counter<-Counter+1
       }
    }
   
    if (Counter>=MinNumberAllelesGenotypedToRetainLocus)  {
      SNPMatrixWithLociWithMinGenotypes<-cbind(SNPMatrixWithLociWithMinGenotypes, c(colnames(DataFrame[f]),SNPMatrix[,f]))
    }
}    
       
#This prints a file that contains each polymorphic site (from biallelic loci) not including indels, that were scored in the user-input proportion of the alleles in the total dataset.    
write.table(SNPMatrixWithLociWithMinGenotypes, file = "out/Output/Genotypes/SNPMatrixTemp_100B.txt", sep = "	",row.names=FALSE, col.names=FALSE)  
########################################################################################################################################################



