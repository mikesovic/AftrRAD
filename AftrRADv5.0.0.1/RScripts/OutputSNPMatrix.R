#We have genotyped all individuals at all loci in file GenotypesUpdate, however, this contains the entire read
#Here, we pull out the SNPs.  Some of the loci in GenotypesUpdate may actually be monomorphic, as they appeared polymorphic based on error, which was eliminated by the binomial test.  These will be eliminated here

DataTableFromFile<-read.table("out/TempFiles/GenotypesUpdate.txt", header=TRUE, sep = "\t")  
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
SNPLocationsMatrix<-matrix(c("Locus","SNPLocation"),ncol=1,nrow=2)

for (a in 2:(NumColsInMatrix-1))  {		#for each locus
    Hit=0	

    
    for (b in 1:NumRowsInMatrix)  {		#get the length of this locus
       if (!is.na(Subset[b,a]))  {
          TempToSplit<-as.character(Subset[b,a])
	      TempSplit<-strsplit(TempToSplit,split="")
	      Delisted<-unlist(TempSplit)
	      Length<-as.numeric(length(Delisted))
	      Hit=1
	      break
       }
    }
     
	if (Hit==1) {		#means there is at least one individual at the current locus with a genotype
    
		Matrix<-matrix(colnames(Subset[a]),nrow=1,ncol=Length)	#create matrix to store genotypes at this locus (each site at the locus is a column in the matrix)
      
		for (b in 1:NumRowsInMatrix)  {					#for each individual in the dataset
	
			if (is.na(Subset[b,a])) {					#if the individual has NA, write "NA" to matrix at every site for this individual
			  Matrix<-rbind(Matrix,Subset[b,a])
		    }
			
			else {										#if the individual has a genotype, write it to the matrix
			  TempToSplit<-as.character(Subset[b,a])
			  TempSplit<-strsplit(TempToSplit,split="")  
			  Delisted<-unlist(TempSplit)
			  Matrix<-rbind(Matrix,Delisted)
		    }  
        }
      
      Matrix
      
      NumRowsInMatrix2<-nrow(Matrix)
      NumColsInMatrix2<-ncol(Matrix)
      
		
		
		for (c in 1:NumColsInMatrix2)  {						#for each site at the current locus
      
      	  Cycle=0
      	  LetterToCheck=0
      	  
			for (d in 2:NumRowsInMatrix2)  {					#for each individual, find the first base at the site
				
				if (!is.na(Matrix[d,c]) && Matrix[d,c]!= "-") {		#if there is a base called at the site (no NA or '-'), store it as CurrentBaseToCheckAgainst, and go to next if statement
                  CurrentBaseToCheckAgainst<-Matrix[d,c]
                  LetterToCheck=1
                  Cycle=d
                  NumericCycle<-as.numeric(Cycle)
                  NumericCyclePlusOne<-NumericCycle+1
                  break
                }   
            
			}
        
            if (LetterToCheck==1)  {							#start from the next individual in the matrix (indexed by NumericCycle) and check the remainder of the bases for variability
               
				for (e in NumericCycle:NumRowsInMatrix2)  {
                 
				    if (!is.na(Matrix[e,c]) && Matrix[e,c]!=CurrentBaseToCheckAgainst && Matrix[e,c]!= "-")  {		#if a new base (polymorphism), add this site as a column in SNPMatrix
                       SNPMatrix<-cbind(SNPMatrix,Matrix[,c])
					   CurrentPosition<-as.numeric(c)	
					   SNPLocationsMatrix<-cbind(SNPLocationsMatrix, c(colnames(Subset[a]),CurrentPosition))
                       break
                    }   
                }
            }   
        }
		
		
		
    } 
    
} 

write.table(SNPLocationsMatrix,file="out/TempFiles/SNPLocations.txt",row.names=FALSE, col.names=FALSE)
write.table(SNPMatrix,file="out/TempFiles/TempSNPMatrix.txt",row.names=FALSE, col.names=FALSE, sep="\t")