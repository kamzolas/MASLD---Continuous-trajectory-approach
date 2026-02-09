

library("preprocessCore")
#install.packages("plotly")
library(tidyverse)
library("sva")

library(readxl)


#merged Ucam/Sanyal
merged_counts <- read.csv("data/mergedcounts.csv", check.names = F)
rownames(merged_counts) = merged_counts[,1]
merged_counts = merged_counts[,-1]

#merged template
template <- read.csv("data/metadata.csv")[,-1]
table(is.element(colnames(merged_counts), template$Sample.name)) #when reading the csv file, it put a "." in the sample names that created confusion - I changed the names to agree!
template = template[ -which(is.element(template$Sample.name, c("Sample 41", "Sample 43"))), ]
merged_counts = merged_counts[ , template$Sample.name]

#Separate in males/females
counts_males = merged_counts[ , template$Sex == "M"]
template_males = template[ template$Sex == "M", ]
rm(merged_counts, template)

#Convert counts to matrix and exclude low expressed counts
counts_males = data.matrix(counts_males)
counts_males <- counts_males[(rowSums(counts_males)>dim(counts_males)[2]),]

#Normalization
nrm_males=normalize.quantiles(log2(1+counts_males))
dimnames(nrm_males) = dimnames(counts_males)



mod_males = model.matrix(~as.factor(template_males$NAS))
corrected_counts_males <- ComBat(dat=as.matrix(nrm_males), batch=template_males$Dataset, mod=mod_males, par.prior=TRUE, prior.plots=FALSE)
write.csv(corrected_counts_males, file = "Batch Corrected Counts (UCAM-Sanyal) - with mod_NAS_only males.csv")

rm(mod_males, counts_males, nrm_males)




# 
# #PCA Plots - Version1
# #####
pca = prcomp(t(corrected_counts_males))
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) #How much variation in the original data each PC accounts for


 
 
 
 
CTRL = which(template_males$SAF.score == "CTRL")
MASL = which(template_males$SAF.score == "MASL")#red
MASH_F012 = which(template_males$SAF.score == "MASH F0" | template_males$SAF.score == "MASH F1" | template_males$SAF.score == "MASH F2") #lightblue
MASH_F34 = which(template_males$SAF.score == "MASH F3" | template_males$SAF.score == "MASH F4") #black
 
 
condition_ = rep("#de2d26", dim(corrected_counts_males)[2]) #MASH_F34
condition_[MASH_F012]="#fb6a4a" #MASH_F0123
condition_[MASL]="#fcae91" #MASL
condition_[CTRL]="#3182bd" #CTRL
 
 

pdf(file="Human_UCAM&SANYAL_PCA_males.pdf")
plot(pca$x[,1], pca$x[,2], col=condition_, xlab=paste("PC1"," - ", pca.var.per[1],"%"), ylab=paste("PC2"," - ", pca.var.per[2],"%"), main = "Normalised Counts\n (Merged UCAM & Sanyal dataset)", pch=16, cex = 1.2)
legend(x="bottomright", legend = c("CTRL", "MASL","MASHF012","MASHF34"), col=c("#3182bd", "#fcae91", "#fb6a4a", "#de2d26"), pch=16, cex=1.2, box.lty = 0)
#text(pca$x[,1], pca$x[,2], colnames(nrm),pos = 1, offset=.1,cex=.2)
dev.off()
 

 
 
 
 



# Find confounding factors for all the batch_corrected counts together
###Part1
# First, I did it for the continuous variables, for "BW", "Liver", "LW.BW.", "Glucose", "ALT", "Insulin", 'TGs', "Cholesterol", "HDL", "LDL", "AST", "ALP", "Hepatic.Cholesterol", "Hepatic.TGs" , "Daily.Food.Intake..average."
# I correlated and created heatmaps using pvalues or cc (correlation value between 0 and 1)
# These are variables that (hopefully) could be correlated with the different PCs (Principal Components)
# I repeated the test for 4 different combinations (1.batch_corrected + phenotypes_without_log2_transformation, 2.batch_corrected + log2_transformed_(phenotypes), 3.normalised + phenotypes_without_log2_transformation), 4.normalised + log2_transformed_(phenotypes)
# Finally, I compare the non corrected vs the corrected  - for the log2 transformed - to verify that the biology is responsible for the variation in the principle components separation
###Part2
# I repeated the same (only the pvalues using anova) for the cotegorical values ("Partner","Diet.Group","Type.of.Cage","Sugar.Water.","Diet.started.at.age..week.","Duration.Experimental.Diet..weeks.")
# In that case I didn't use the log of the codes, cause they are categorial values
library("preprocessCore")


corrected <- read.csv(file = "Batch Corrected Counts (UCAM-Sanyal) - with mod_NAS_only males.csv", sep = ',', header = TRUE, check.names = F)
rownames(corrected) = corrected[,1]
corrected = corrected[,-1]

codes <- template_males
colnames(corrected) = codes$Sample.name
unique(colnames(corrected) == codes$Sample.name)
codes = codes[,-3] #remove the sex variable

for (j in 1:dim(codes)[2]) 
{
  print(colnames(codes)[j])
  print(summary(codes[,j]))
}


#For the continuous variables
codes = codes[ , c( "Age", "NAS", "Fibrosis", "Fat", "Ballooning", "Inflammation")]
#codes = log2(codes) #### ### ### ### ### ### To perform a parametric test, the distribution has to be Normal. A way to achive that sometimes is via log-transformation


#PCA Plots
pca = prcomp(t(corrected))

confactors= correlation = matrix(0, nrow = 10, ncol = dim(codes)[2])
for (i in 1:10)
  for (j in 1:dim(codes)[2])
  {
    a = cor.test(pca$x[,i], codes[,j], method = "pearson")
    correlation[i,j] = a$estimate
    confactors[i,j] = a$p.value
    #confactors[i,j] = - log10(a$p.value +.01)
  }
rownames(confactors) =  rownames(correlation) = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
colnames(confactors) = colnames(correlation)  =colnames(codes)




colored_pvals = confactors

library(pheatmap)
paletteLength <- 21
myColor <- colorRampPalette(c( "red","white"))(paletteLength)
myBreaks <- seq(0, 1, by = .01)
pheatmap(colored_pvals, color=myColor, cluster_rows = FALSE, breaks=myBreaks, display_numbers = round(confactors,3), cluster_cols=FALSE, fontsize_row = 10, fontsize_col = 10, main = "P-value for continuous variables\n(mod = NAS) - Males", file = "P-vals_males.png")

myColor <- colorRampPalette(c( "green","white", "red"))(paletteLength)
myBreaks <- seq(-1, 1, by = .1)
pheatmap(correlation, color=myColor, cluster_rows = FALSE, breaks=myBreaks, display_numbers = round(correlation,3), cluster_cols=FALSE, fontsize_row = 10, fontsize_col = 10, main = "Correlation Plot for continuous variables\n(mod = NAS) - Males", file = "Correlation_males.png")


#I performed the Shapiro-Wilk's method to test for Normality. If the p-value<.05 it is assumed that the distribution is not normal

for (i in 1:dim(codes)[2])
{
  print(colnames(codes)[i])
  test = shapiro.test(codes[,i])
  print(test$p.value)
  if(test$p.value >.05)
    print("Normal Distribution")
  else
    print("Non-Normality")
}










###For the categorical variable
codes = template_males[ , c("Age","NAS","T2DM","Dataset")]
allcodes =template_males

#
codes$Age[allcodes$Age<45] <- "young"
codes$Age[allcodes$Age>=45 & codes$Age<60] <- "middle"
codes$Age[allcodes$Age>=60] <- "old"


for (j in 1:dim(codes)[2]) 
{
  print(colnames(codes)[j])
  print(summary(codes[,j]))
  print(typeof(codes[,j]))
}


#PCA Plots
pca = prcomp(t(corrected))

confactors = matrix(2, nrow = 10, ncol = dim(codes)[2]) #The matrix contains the pvals
for (i in 1:10)
  for (j in 1:dim(codes)[2])
  {
    df = data.frame(codes[,j], pca$x[,i])
    age.aov = aov(pca$x[,i]~codes[,j], data=df)
    s = unlist(summary(age.aov))
    confactors[i,j] = s[9]
  }
rownames(confactors) = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
colnames(confactors) =colnames(codes)

library(pheatmap)
paletteLength <- 21
myColor <- colorRampPalette(c( "red","white"))(paletteLength)
myBreaks <- seq(0, 1, by = .01)
pheatmap(confactors, color=myColor, cluster_rows = FALSE, breaks=myBreaks, display_numbers = round(confactors,5), cluster_cols=FALSE, fontsize_row = 10, fontsize_col = 10, main = "P-value for categorical variables\n(mod = NAS) - Males", file = "P-vals_categorical variables_males.png")





