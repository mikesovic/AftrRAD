DepthsToPlot<-scan(file="out/TempFiles/DepthsToPlot.txt")
MaxLength<-100
Interval<-MaxLength/100
pdf("out/Output/RunInfo/ReadDepths.pdf")
hist(DepthsToPlot, breaks=seq(0,MaxLength,by=Interval))