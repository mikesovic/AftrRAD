DataTableFromFile<-read.table("out/TempFiles/SingleSNPsAll.txt", header=TRUE, sep = "\t")

DataMatrix<-as.matrix(DataTableFromFile)
DataFrame<-data.frame(DataMatrix)
NumRowsInMatrix<-nrow(DataMatrix)
NumColsInMatrix<-ncol(DataMatrix)

NamesCol<-DataMatrix[,1]

BiallelicMatrix<-matrix(c("",NamesCol), ncol=1, nrow=NumRowsInMatrix+1)

for (a in 2:(NumColsInMatrix))  {	#Check each locus one at a time

     
      CurrentLocus<-colnames(DataFrame[a])
      
      UniqueHaplotypes<-"starter"
      
      
      
      for (b in 1:NumRowsInMatrix)  {	#Start going through all snps at this locus
      		Match<-0
      		
             	for (c in UniqueHaplotypes)  {
                  
                   if (DataFrame[b,a] == c) {
                       Match<-1
                       break
                   }
                }   
                   
                if (Match == 0)  {
                   
                   	UniqueHaplotypes<-c(UniqueHaplotypes,as.character(DataFrame[b,a]))
                }
               
                
      }		#Finished going through all snps at this locus
      
      
      NumberOfUniqueHaplotypes<-length(UniqueHaplotypes)
      VectorOfUniqueHaplotypes<-UniqueHaplotypes[2:NumberOfUniqueHaplotypes]

      
      if (NumberOfUniqueHaplotypes == 3)  {
      	   VectorToPrint<-c(CurrentLocus, as.character(DataFrame[,a]))
      	   BiallelicMatrix<-cbind(BiallelicMatrix,VectorToPrint)
      }	
      
}  #Finished with current locus

write.table(BiallelicMatrix, file = "TempFiles/AllBiallelicSNPsRaw.txt", sep = "\t",row.names=FALSE, col.names=FALSE)






UnlinkedBiallelicMatrix<-matrix(c("",NamesCol), ncol=1, nrow=NumRowsInMatrix+1)

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
		UnlinkedBiallelicMatrix<-cbind(UnlinkedBiallelicMatrix,VectorToPrint)
	}
	
}

write.table(UnlinkedBiallelicMatrix, file = "out/TempFiles/UnlinkedBiallelicSNPsRaw.txt", sep = "\t",row.names=FALSE, col.names=FALSE)


	
