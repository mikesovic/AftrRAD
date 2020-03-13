Values<-scan(file="out/TempFiles/NACountsToPlot.txt")
Names<-scan(file="out/TempFiles/NASampleNamesToPlot.txt", what=character())
pdf("out/Output/RunInfo/MissingDataCounts.pdf")
barplot(Values,names.arg=Names,las=2)
dev.off()

