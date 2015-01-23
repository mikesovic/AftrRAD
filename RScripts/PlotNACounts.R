Values<-scan(file="TempFiles/NACountsToPlot.txt")
Names<-scan(file="TempFiles/NASampleNamesToPlot.txt", what=character())
pdf("Output/RunInfo/MissingDataCounts.pdf")
barplot(Values,names.arg=Names,las=2)
dev.off()

