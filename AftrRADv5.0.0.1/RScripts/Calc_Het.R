Monos<-read.table(file="R_Mono_Infile.txt", sep=" ", header=FALSE)

NumCols<-ncol(Monos)

MinReads<-edit

NumMonosScoredPerSample<-c()

for (i in 1:NumCols) {
  CurrentCol<-Monos[,i]
  Subset<-c(CurrentCol[CurrentCol>=MinReads])
  len<-length(Subset)
  NumMonosScoredPerSample<-c(NumMonosScoredPerSample,len)
}

write.table(NumMonosScoredPerSample, file="MonosOut.txt", row.names=FALSE, col.names=FALSE)
