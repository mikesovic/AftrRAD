LocationsToPlot<-scan(file="out/TempFiles/SNPLocationsToPlot.txt")
MaxLength<-max(LocationsToPlot)
pdf("out/Output/RunInfo/SNPLocations.pdf")
hist(LocationsToPlot, breaks=seq(0,MaxLength,by=1))