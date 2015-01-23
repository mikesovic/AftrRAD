LocationsToPlot<-scan(file="TempFiles/SNPLocationsToPlot.txt")
MaxLength<-max(LocationsToPlot)
pdf("Output/RunInfo/SNPLocations.pdf")
hist(LocationsToPlot, breaks=seq(0,MaxLength,by=1))