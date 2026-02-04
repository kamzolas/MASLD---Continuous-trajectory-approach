# Load required scripts and libraries
source("~/Desktop/DESeq2_Analysis_Functions.R")
source("~/Desktop/Plotly_plots.R")
library(readxl)
library("ggpubr")
library("cowplot")

# Load datasets
patients_in_SWs = read.csv("data/Patients in SWs.csv")

template = read.csv(file = "data/metadata.csv")[, -c(1,3)]
cts_ucamsanyal <- read.csv(file = "data/mergedcounts.csv", sep = ",")
rownames(cts_ucamsanyal) = as.character(cts_ucamsanyal$X)
cts_ucamsanyal = cts_ucamsanyal[,-c(1)]
colnames(cts_ucamsanyal) = template$Sample.name
cts_ucamsanyal <- cts_ucamsanyal[(rowSums(cts_ucamsanyal)>dim(cts_ucamsanyal)[2]),] #Exclude low expressed counts

# Initialize vectors
meanNASinSWs = meanSteatosisinSWs = meanInflammationinSWs = meanBallooninginSWs = meanFibrosisinSWs = SW_number = rep(-2, dim(patients_in_SWs)[2])
semNAS = semSteatosis = semInflammation = semBallooning = semFibrosis = rep(NA, dim(patients_in_SWs)[2])

# Compute means and SEMs
for (i in 1:dim(patients_in_SWs)[2]) {
   patients_in_sw = patients_in_SWs[,i]
   patients_in_sw = patients_in_sw[patients_in_sw != ""]
   
   idx = is.element(template$Sample.name, patients_in_sw)
   
   meanNASinSWs[i] = mean(template$NAS[idx])
   semNAS[i] = sd(template$NAS[idx]) / sqrt(sum(idx))
   
   meanSteatosisinSWs[i] = mean(template$Fat[idx])
   semSteatosis[i] = sd(template$Fat[idx]) / sqrt(sum(idx))
   
   meanInflammationinSWs[i] = mean(template$Inflammation[idx])
   semInflammation[i] = sd(template$Inflammation[idx]) / sqrt(sum(idx))
   
   meanBallooninginSWs[i] = mean(template$Ballooning[idx])
   semBallooning[i] = sd(template$Ballooning[idx]) / sqrt(sum(idx))
   
   meanFibrosisinSWs[i] = mean(template$Fibrosis[idx])
   semFibrosis[i] = sd(template$Fibrosis[idx]) / sqrt(sum(idx))
}

# Plot NAS score with error bars
pdf("Average NAS score.pdf", width=8, height=3)
options(repr.plot.width=8, repr.plot.height=2)
range_meanNASinSWs = range(meanNASinSWs + semNAS, meanNASinSWs - semNAS)
range_meanNASinSWs[1] = range_meanNASinSWs[1] - .5
range_meanNASinSWs[2] = range_meanNASinSWs[2] + .5

plot(1, type="n", xlab="", ylab="", xlim=c(1, length(meanNASinSWs)), ylim=range_meanNASinSWs, main="", xaxt="n")

# Add lines and points
lines(meanNASinSWs, type="o", pch=16, col="darkblue", lwd=2)
points(meanNASinSWs, pch=16, col="darkblue", cex=1.5)

# Add error bars
arrows(x0=1:length(meanNASinSWs), y0=meanNASinSWs - semNAS,
       x1=1:length(meanNASinSWs), y1=meanNASinSWs + semNAS,
       angle=90, code=3, length=0.05, col="darkblue", lwd=1.5)

# Customize axis
axis(1, at=1:length(meanNASinSWs), labels=1:length(meanNASinSWs))
mtext("NAS score", side=2, line=2.5)

dev.off()




pdf("Average SBIF with error bars.pdf", width=8, height=3)
# Set up the plotting area
options(repr.plot.width=8, repr.plot.height=2)
plot(1, type="n", xlab="", ylab="", xlim=c(1, length(meanSteatosisinSWs)), ylim=c(0, 4), main="", xaxt="n")

# Add lines and points
lines(meanSteatosisinSWs, type="o", pch=16, col="#a6611a", lwd=2)
points(meanSteatosisinSWs, pch=16, col="#a6611a", cex=1.5)

lines(meanBallooninginSWs, type="o", pch=16, col="#80cdc1", lwd=2)
points(meanBallooninginSWs, pch=16, col="#80cdc1", cex=1.5)

lines(meanInflammationinSWs, type="o", pch=16, col="#dfc27d", lwd=2)
points(meanInflammationinSWs, pch=16, col="#dfc27d", cex=1.5)

lines(meanFibrosisinSWs, type="o", pch=16, col="#018571", lwd=2)
points(meanFibrosisinSWs, pch=16, col="#018571", cex=1.5)

# Add error bars
arrows(x0=1:13, y0=meanSteatosisinSWs - semSteatosis, y1=meanSteatosisinSWs + semSteatosis,
       code=3, angle=90, length=0.05, col="#a6611a")
arrows(x0=1:13, y0=meanBallooninginSWs - semBallooning, y1=meanBallooninginSWs + semBallooning,
       code=3, angle=90, length=0.05, col="#80cdc1")
arrows(x0=1:13, y0=meanInflammationinSWs - semInflammation, y1=meanInflammationinSWs + semInflammation,
       code=3, angle=90, length=0.05, col="#dfc27d")
arrows(x0=1:13, y0=meanFibrosisinSWs - semFibrosis, y1=meanFibrosisinSWs + semFibrosis,
       code=3, angle=90, length=0.05, col="#018571")

# Axes and labels
axis(1, at=1:length(meanNASinSWs), labels=1:length(meanNASinSWs))
mtext("Score", side=2, line=2.5)

# legend("topright", legend=c("Steatosis", "Ballooning", "Inflammation", "Fibrosis"),
#        col=c("#a6611a", "#80cdc1", "#dfc27d", "#018571"), pch=16, lwd=2, cex=0.8)

dev.off()








pdf("Average SBIF 4 panels with error bars.pdf", width=6, height=8)
par(mfrow=c(4,1), mar=c(3, 4, 2, 2))  # 4 rows, 1 column; adjust margins

# Steatosis
plot(meanSteatosisinSWs, type="o", pch=16, col="#a6611a", lwd=2,
     ylab="Steatosis", xlab="", ylim=c(0,4), xaxt="n")
arrows(x0=1:13, y0=meanSteatosisinSWs - semSteatosis, y1=meanSteatosisinSWs + semSteatosis,
       code=3, angle=90, length=0.05, col="#a6611a")
axis(1, at=1:13, labels=1:13)

# Ballooning
plot(meanBallooninginSWs, type="o", pch=16, col="#80cdc1", lwd=2,
     ylab="Ballooning", xlab="", ylim=c(0,4), xaxt="n")
arrows(x0=1:13, y0=meanBallooninginSWs - semBallooning, y1=meanBallooninginSWs + semBallooning,
       code=3, angle=90, length=0.05, col="#80cdc1")
axis(1, at=1:13, labels=1:13)

# Inflammation
plot(meanInflammationinSWs, type="o", pch=16, col="#dfc27d", lwd=2,
     ylab="Inflammation", xlab="", ylim=c(0,4), xaxt="n")
arrows(x0=1:13, y0=meanInflammationinSWs - semInflammation, y1=meanInflammationinSWs + semInflammation,
       code=3, angle=90, length=0.05, col="#dfc27d")
axis(1, at=1:13, labels=1:13)

# Fibrosis
plot(meanFibrosisinSWs, type="o", pch=16, col="#018571", lwd=2,
     ylab="Fibrosis", xlab="Sliding window", ylim=c(0,4))
arrows(x0=1:13, y0=meanFibrosisinSWs - semFibrosis, y1=meanFibrosisinSWs + semFibrosis,
       code=3, angle=90, length=0.05, col="#018571")
axis(1, at=1:13, labels=1:13)

dev.off()





