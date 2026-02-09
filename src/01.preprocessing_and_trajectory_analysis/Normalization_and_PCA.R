#Normalization, PCA and Screeplot
#Ucam & VCU (Sanyal) datasets (human MASLD datasets)

library("preprocessCore")
library(readxl)

#Convert counts to matrix and exlude low expressed counts
merged_counts = read.csv(file = "../data/mergedcounts.csv", sep = ",")
rownames(merged_counts) = merged_counts$X
merged_counts = as.matrix(merged_counts[,-1])

#Normalization
pdf(file="Human_Boxplot of counts.pdf")
par(mfrow=c(1,2), mar = c(4,3,2,4))
boxplot(log2(1+merged_counts), main="Un-Normalized Counts")
nrm=normalize.quantiles(log2(1+merged_counts))
dimnames(nrm) = dimnames(merged_counts)
boxplot((nrm), main="Normalized Counts")
write.csv(nrm, file ="Human_Normalised_counts.csv")
dev.off()

#Normalization, PCA and Screeplot
template = read.csv(file = "../data/metadata.csv", header = TRUE)[,-1]

unique(template$Sample.name == colnames(nrm))
CTRL = which(template$SAF.score == "CTRL")
MASL = which(template$SAF.score == "MASL")#red
MASH_F012 = which(template$SAF.score == "MASH F0" | template$SAF.score == "MASH F1" | template$SAF.score == "MASH F2") #lightblue
MASH_F34 = which(template$SAF.score == "MASH F3" | template$SAF.score == "MASH F4") #black

#PCA plot
#Points for the different stages of the PCA plot
pca = prcomp(t(nrm))
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) #How much variation in the original data each PC accounts for

condition_ = rep("#de2d26", dim(nrm)[2]) #MASH_F34
condition_[MASH_F012]="#fb6a4a" #MASH_F0123
condition_[MASL]="#fcae91" #MASL
condition_[CTRL]="#3182bd" #CTRL

pdf(file="Human_UCAM&VCU_PCA.pdf")
plot(pca$x[,1], pca$x[,2], col=condition_, xlab=paste("PC1"," - ", pca.var.per[1],"%"), ylab=paste("PC2"," - ", pca.var.per[2],"%"), main = "Normalised Counts\n (Merged UCAM & VCU dataset)", pch=16, cex = 1.2)
legend(x="top", legend = c("CTRL", "MASL","MASHF012","MASHF34"), col=c("#3182bd", "#fcae91", "#fb6a4a", "#de2d26"), pch=16, cex=1.2, box.lty = 0)
dev.off()

pdf(file="Human_ScreePlot.pdf")
barplot(pca.var.per, main="Scree Plot", xlab="Principal Components", ylab="Percent Variation", names.arg = pca.var.per)
dev.off()




