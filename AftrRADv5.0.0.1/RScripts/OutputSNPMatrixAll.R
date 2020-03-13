DataTableFromFile<-read.table("out/TempFiles/Genotypes.txt", header=TRUE, sep = "\t")
DataMatrix<-as.matrix(DataTableFromFile)
DataFrame<-data.frame(DataMatrix)
NumRowsInMatrix<-nrow(DataMatrix)
NumColsInMatrix<-ncol(DataMatrix)

Subset<-subset(DataFrame[,1:(NumColsInMatrix-1)])
attach(Subset)

######################################################################################
FirstColToSplit<-as.character(Subset[,1])
FirstColSplit<-strsplit(FirstColToSplit,split="\t")
SplitColDelisted<-unlist(FirstColSplit)
SNPMatrix<-matrix(c("",SplitColDelisted),ncol=1,nrow=(NumRowsInMatrix+1))

for (a in 2:(NumColsInMatrix-1))  {
    Hit=0	

    
    
    for (b in 1:NumRowsInMatrix)  {
       
		if (!is.na(Subset[b,a]))  {
          TempToSplit<-as.character(Subset[b,a])
	      TempSplit<-strsplit(TempToSplit,split="")
	      Delisted<-unlist(TempSplit)
	      Length<-as.numeric(length(Delisted))
	      Hit=1
	      break
       }
    }
     
   if (Hit==1) {
    
      Matrix<-matrix(colnames(Subset[a]),nrow=1,ncol=Length)
      
      for (b in 1:NumRowsInMatrix)  {
	
	    if (is.na(Subset[b,a])) {
	        Matrix<-rbind(Matrix,Subset[b,a])
	    }
	  
	    else {
	        TempToSplit<-as.character(Subset[b,a])
	        TempSplit<-strsplit(TempToSplit,split="")  
	        Delisted<-unlist(TempSplit)
	        Matrix<-rbind(Matrix,Delisted)
	    }  
      }
      
      Matrix
      
      NumRowsInMatrix2<-nrow(Matrix)
      NumColsInMatrix2<-ncol(Matrix)
      
      for (c in 1:NumColsInMatrix2)  { 
      
      	  Cycle=0
      	  LetterToCheck=0
      	  
          for (d in 2:NumRowsInMatrix2)  {
             if (!is.na(Matrix[d,c]) && Matrix[d,c]!= "-") {
                CurrentBaseToCheckAgainst<-Matrix[d,c]
                LetterToCheck=1
                Cycle=d
                NumericCycle<-as.numeric(Cycle)
                NumericCyclePlusOne<-NumericCycle+1
                break
             }   
          }
        
          if (LetterToCheck==1)  {
             for (e in NumericCycle:NumRowsInMatrix2)  {
               if (!is.na(Matrix[e,c]) && Matrix[e,c]!=CurrentBaseToCheckAgainst && Matrix[e,c]!= "-")  {
                  SNPMatrix<-cbind(SNPMatrix,Matrix[,c])
                  break
               }   
             }
          }   
      }
   } 
    
} 


write.table(SNPMatrix, file = "out/TempFiles/SNPMatrixAll.txt", sep = "\t",row.names=FALSE, col.names=FALSE)       

########################################################################################################################################################
