
Genotypes<-read.table(file="out/TempFiles/TempBiallelicLociForGenotyping.txt")
GenotypesMatrix<-as.matrix(Genotypes)

NumRowsInGenotypesMatrix=nrow(GenotypesMatrix)
NumColsInGenotypesMatrix=ncol(GenotypesMatrix)
NumLociInMatrix=NumRowsInGenotypesMatrix/3


for (j in 1:NumLociInMatrix) {
    CurrentThirdRow<-j*3
    CurrentSecondRow<-CurrentThirdRow-1
    CurrentFirstRow<-CurrentThirdRow-2
    
    NumTrials1 <- GenotypesMatrix[CurrentSecondRow,2]
    NumTrials1 <- as.numeric(NumTrials1)
    NumTrials2 <- GenotypesMatrix[CurrentThirdRow,2]
    NumTrials2 <- as.numeric(NumTrials2)
    NumTrials <- NumTrials1 + NumTrials2
    NumTrials <- as.numeric(NumTrials)
    
    if (NumTrials<5) {
       write(GenotypesMatrix[CurrentFirstRow,1],append=TRUE,file="TempGenotypes.txt") 
       write("NA",append=TRUE,file="TempGenotypes.txt") 
       write("NA",append=TRUE,file="TempGenotypes.txt")
    }
    else if (GenotypesMatrix[CurrentSecondRow,2]==0) {
       write(GenotypesMatrix[CurrentFirstRow,1],append=TRUE,file="TempGenotypes.txt") 
       write(GenotypesMatrix[CurrentThirdRow,1],append=TRUE,file="TempGenotypes.txt") 
       write(GenotypesMatrix[CurrentThirdRow,1],append=TRUE,file="TempGenotypes.txt")
    }
    else if (GenotypesMatrix[CurrentThirdRow,2]==0) {
       write(GenotypesMatrix[CurrentFirstRow,1],append=TRUE,file="TempGenotypes.txt") 
       write(GenotypesMatrix[CurrentSecondRow,1],append=TRUE,file="TempGenotypes.txt") 
       write(GenotypesMatrix[CurrentSecondRow,1],append=TRUE,file="TempGenotypes.txt")
    }
    
    else  {
       SecondRowCount<-as.numeric(GenotypesMatrix[CurrentSecondRow,2])
       ThirdRowCount<-as.numeric(GenotypesMatrix[CurrentThirdRow,2])
       Counts<-c(SecondRowCount,ThirdRowCount) 
       NumSuccesses<-min(Counts) 
       NumSuccesses<-as.numeric(NumSuccesses)
       CumProb<-pbinom(NumSuccesses,NumTrials,0.5) 
    
		if (NumTrials<100) {	
		
			if (CumProb<0.00001) { #Have a homozygoteLow
           
          if (SecondRowCount < ThirdRowCount)  {
		    write(GenotypesMatrix[CurrentFirstRow,1],append=TRUE,file="TempGenotypes.txt")        
            write(GenotypesMatrix[CurrentThirdRow,1],append=TRUE,file="TempGenotypes.txt")
            write(GenotypesMatrix[CurrentThirdRow,1],append=TRUE,file="TempGenotypes.txt")
		  }
          
		  else { 
           write(GenotypesMatrix[CurrentFirstRow,1],append=TRUE,file="TempGenotypes.txt")        
           write(GenotypesMatrix[CurrentSecondRow,1],append=TRUE,file="TempGenotypes.txt")
           write(GenotypesMatrix[CurrentSecondRow,1],append=TRUE,file="TempGenotypes.txt") 
		  }
		 }   
       
		 else  { #Have a heterozygote
           write(GenotypesMatrix[CurrentFirstRow,1],append=TRUE,file="TempGenotypes.txt") 
           write(GenotypesMatrix[CurrentSecondRow,1],append=TRUE,file="TempGenotypes.txt")
           write(GenotypesMatrix[CurrentThirdRow,1],append=TRUE,file="TempGenotypes.txt")
         }
		}
		
		else { 
			
			if (CumProb<0.000001) { #Have a homozygoteHigh
				
				if (SecondRowCount < ThirdRowCount)  {
					write(GenotypesMatrix[CurrentFirstRow,1],append=TRUE,file="TempGenotypes.txt")        
					write(GenotypesMatrix[CurrentThirdRow,1],append=TRUE,file="TempGenotypes.txt")
					write(GenotypesMatrix[CurrentThirdRow,1],append=TRUE,file="TempGenotypes.txt")
				}
				
				else { 
					write(GenotypesMatrix[CurrentFirstRow,1],append=TRUE,file="TempGenotypes.txt")        
					write(GenotypesMatrix[CurrentSecondRow,1],append=TRUE,file="TempGenotypes.txt")
					write(GenotypesMatrix[CurrentSecondRow,1],append=TRUE,file="TempGenotypes.txt") 
				}
			}   
			
			else  { #Have a heterozygote
				write(GenotypesMatrix[CurrentFirstRow,1],append=TRUE,file="TempGenotypes.txt") 
				write(GenotypesMatrix[CurrentSecondRow,1],append=TRUE,file="TempGenotypes.txt")
				write(GenotypesMatrix[CurrentThirdRow,1],append=TRUE,file="TempGenotypes.txt")
			}
		}
			
}
}
