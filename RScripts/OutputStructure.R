DataTableFromFile<-read.table("OutputStructure_Infile.txt", header=TRUE, sep = "\t")

DataMatrix<-as.matrix(DataTableFromFile)
DataFrame<-data.frame(DataMatrix)
NumRowsInMatrix<-nrow(DataMatrix)
NumColsInMatrix<-ncol(DataMatrix)

#Subset<-subset(DataFrame[,1:(NumColsInMatrix-1)])  #create data frame that has only the haplotypes (remove names col, and locus names are treated as header row)
#attach(Subset)

NamesCol<-DataMatrix[,1]

StructureMatrix<-matrix(c("",NamesCol), ncol=1, nrow=NumRowsInMatrix+1)

for (a in 2:(NumColsInMatrix))  {

     
      CurrentLocus<-colnames(DataFrame[a])
      
      UniqueHaplotypes<-"N"
      
      
	  #For the current locus (matrix column), determine how many unique haplotypes occur.
	
      for (b in 1:NumRowsInMatrix)  {
      		Match<-0
      		
             	for (c in UniqueHaplotypes)  {
                  
                   if (as.character(DataFrame[b,a]) == c) {
                       Match<-1
                       break
                   }
                }   
                   
                if (Match == 0)  {
                   
                   	UniqueHaplotypes<-c(UniqueHaplotypes,as.character(DataFrame[b,a]))
                }
                       
      }
      
      
      NumberOfUniqueHaplotypes<-length(UniqueHaplotypes)
      VectorOfUniqueHaplotypes<-UniqueHaplotypes[2:NumberOfUniqueHaplotypes]

      CurrentLocusHaps<-c(DataMatrix[,a])
      
      CurrentLocusColToPrint<-CurrentLocus
      
      for (hap in CurrentLocusHaps)  {
		  
		  #assign -9 to missing data.
		  
		  if (hap == "N") {
			  CurrentLocusColToPrint<-c(CurrentLocusColToPrint,"-9")
		  }	  
		  
		  
		  else {
			  
			  UniqueLocation<-0
      
			  
			  #Assign each haplotype a number, and print that to the StructureMatrix.
			  
			  for (unique in VectorOfUniqueHaplotypes)  {
              
				  
               if (hap == unique)  {
                   CurrentLocusColToPrint<-c(CurrentLocusColToPrint,UniqueLocation)
                   break
               }   
          
               else  {
               
                   UniqueLocation<-UniqueLocation+1
               }
              }    
		  } 
      }
      
      
      StructureMatrix<-cbind(StructureMatrix,CurrentLocusColToPrint)
      
} 
   

write.table(StructureMatrix, file = "StructureInput_Temp.txt", sep = "\t",row.names=FALSE, col.names=FALSE)
