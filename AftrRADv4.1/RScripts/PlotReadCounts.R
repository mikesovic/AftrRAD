DepthsToPlot<-scan(file="TempFiles/DepthsToPlot.txt")
MaxLength<-100
Interval<-MaxLength/100
pdf("Output/RunInfo/ReadDepths.pdf")
hist(DepthsToPlot, breaks=seq(0,MaxLength,by=Interval))