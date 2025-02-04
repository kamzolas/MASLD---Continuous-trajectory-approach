##Analysis was based on the example of this Vignette:
#https://www.bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("slingshot")

library("slingshot")
library("SingleCellExperiment")

counts=read.csv("../data/batch_corrected_counts_(dataset+gender).csv", sep = ",")
rownames(counts) = counts$X
counts = counts[ ,-c(1)]
counts = round(counts)
template=read.csv(file = "../data/metadata.csv",header = TRUE ,stringsAsFactors = FALSE,row.names = 1)
rownames(template) = template$Sample.name
colnames(counts) = template$Sample_names
sim <- SingleCellExperiment(assays = List(counts = as.matrix(counts)))

#Dimensionality Reduction
pca <- prcomp(t(assays(sim)$counts), scale. = FALSE)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) #How much variation in the original data each PC accounts for
rd1 <- pca$x[,1:2]

rescale <- function(rd1)
{
  max_ = max(abs(rd1[,2]))
  for (i in 1:length(rd1[,2]))
  {
    if(abs(rd1[i,2]) > max_/3) rd1[i,2] = rd1[i,2] * .2
    else if(abs(rd1[i,2]) > max_/5) rd1[i,2] = rd1[i,2] * .4
    else if(abs(rd1[i,2]) > max_/10) rd1[i,2] = rd1[i,2] * .6
  }
  return(rd1)
}

rd1 = rescale(rd1)
reducedDims(sim) <- SimpleList(PCA = rd1)

#Clustering Cells
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(sim)$GMM <- cl1

library(RColorBrewer)
#plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 7)$cluster
colData(sim)$kmeans <- cl2
#plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)


#Using Slingshot
#Uses either "GMM" or "kmeans" - kmeans will have a number of clusters = centers
sim <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA')
#sim <- slingshot(sim, clusterLabels = 'kmeans', reducedDim = 'PCA')

#change colors for plot
col_ = template$NAS
col_[template$NAS ==0] = "#081d58"
col_[template$NAS ==1] = "#253494"
col_[template$NAS ==2] = "#225ea8"
col_[template$NAS ==3] = "#1d91c0"
col_[template$NAS ==4] = "#41b6c4"
col_[template$NAS ==5] = "#7fcdbb"
col_[template$NAS ==6] = "#c7e9b4"
col_[template$NAS ==7] = "#edf8b1"
col_[template$NAS ==8] = "#ffffd9"

pdf(file="Pseudo_UCAM&VCU_aftercorrection.pdf")
plot(reducedDims(sim)$PCA[,1], reducedDims(sim)$PCA[,2] * (pca.var.per[2] / pca.var.per[1]), col = col_, pch=16, cex = .6,
     ylim=c( -max(abs(reducedDims(sim)$PCA[,1])) - 20, max(abs(pca$x[,1])) +10),
     main = "Pseudotemporal Ordering \n UCAM/VCU (After batch correction)", xlab="Pseudotime", ylab=" ")
lines(SlingshotDataSet(sim), lwd=2, col='black')
legend("top", legend=rev(c("0", "1", "2", "3", "4", "5", "6", "7")), horiz = T,
       col=rev(c("#081d58","#253494","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#edf8b1")), pch=20, cex=0.78, bty = "n")
dev.off()



