

library("preprocessCore")
#install.packages("plotly")
library(plotly)
library(tidyverse)
library("sva")
library(readxl)

#merged Ucam/VCU counts
merged_counts <- read.csv("../data/mergedcounts.csv")
rownames(merged_counts) = merged_counts$X
merged_counts = merged_counts[,-1]

#remove outlier - Sample 5
merged_counts = merged_counts[ , -5]

#merged template
template <- read.csv("../data/metadata.csv")[,-1]
template = template[-5,] #remove outlier
colnames(merged_counts) = template$Sample.name

#Convert counts to matrix and exlude low expressed counts
merged_counts = data.matrix(merged_counts)
merged_counts <- merged_counts[(rowSums(merged_counts)>dim(merged_counts)[2]),]

#Normalization
nrm=normalize.quantiles(log2(1+merged_counts))
dimnames(nrm) = dimnames(merged_counts)

batch_ = paste(template$Dataset, template$Sex)
mod_ = model.matrix(~as.factor(template$NAS))
corrected <- ComBat(dat=as.matrix(nrm), batch=batch_, mod=mod_, par.prior=TRUE, prior.plots=FALSE)
write.csv(corrected, file = "batch_corrected_counts_(dataset+gender).csv")

unique(template$Sample.name == colnames(corrected))

CTRL = which(template$SAF.score == "CTRL")
MASL = which(template$SAF.score == "MASL")#red
MASH_F012 = which(template$SAF.score == "MASH F0" | template$SAF.score == "MASH F1" | template$SAF.score == "MASH F2") #lightblue
MASH_F34 = which(template$SAF.score == "MASH F3" | template$SAF.score == "MASH F4") #black



#PCA plot --- After batch correction
pca = prcomp(t(corrected))
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) #How much variation in the original data each PC accounts for

condition_ = rep("#de2d26", dim(corrected)[2]) #MASH_F34
condition_[MASH_F012]="#fb6a4a" #MASH_F0123
condition_[MASL]="#fcae91" #MASL
condition_[CTRL]="#3182bd" #CTRL

pdf(file="Human_UCAM&VCU_PCA_afterbatchcorrection_excludingoutlier.pdf")
plot(pca$x[,1], pca$x[,2], col=condition_, xlab=paste("PC1"," - ", pca.var.per[1],"%"), ylab=paste("PC2"," - ", pca.var.per[2],"%"), main = "Batch Corrected Counts\n (Merged UCAM & VCU dataset - excluding outlier)", pch=16, cex = 1.2)
legend(x="topleft", legend = c("CTRL", "MASL","MASHF012","MASHF34"), col=c("#3182bd", "#fcae91", "#fb6a4a", "#de2d26"), pch=16, cex=1.2, box.lty = 0)
dev.off()

pdf(file="Human_ScreePlot.pdf")
barplot(pca.var.per, main="Scree Plot", xlab="Principal Components", ylab="Percent Variation", names.arg = pca.var.per)
dev.off()




