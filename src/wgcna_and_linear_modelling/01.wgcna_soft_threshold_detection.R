suppressMessages(library(readr)) # 2.1.4
suppressMessages(library(WGCNA)) # 1.72.1
suppressMessages(library(ggplot2)) # 3.4.2


output_dir <- '../../results/wgcna_and_linear_modelling/'
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

################################################################################ 
# Load ensembl mapping and get the rib, mito and unannotated genes
################################################################################
remove_genes <- FALSE
ensembl_mapping <- read.table('../../data/ensembl_mapping.tsv', sep='\t', header=TRUE, row.names=1)
if (remove_genes == TRUE) {
  ribo_genes <- ensembl_mapping[grepl('^RP', ensembl_mapping$external_gene_name),'ensembl_gene_id']
  mt_genes <- ensembl_mapping[grepl('^MT-', ensembl_mapping$external_gene_name),'ensembl_gene_id']
  mrp_genes <- ensembl_mapping[grepl('^MRP', ensembl_mapping$external_gene_name),'ensembl_gene_id']
  nan_genes <- ensembl_mapping[ensembl_mapping$external_gene_name == "",'ensembl_gene_id']
  to_remove_genes <- Reduce(union, list(ribo_genes, mt_genes, mrp_genes, nan_genes))
} else {
  to_remove_genes <- ensembl_mapping[ensembl_mapping$external_gene_name == "",'ensembl_gene_id']
}
print(length(to_remove_genes))


################################################################################
# Load the expression dataset and format it properly
################################################################################
GeneXData <- read_csv("../../data/batch_corrected_counts_(dataset+gender).csv")
colnames(GeneXData)[1] <- "GeneID"
gene.IDs <- GeneXData$GeneID
GeneXData <- as.data.frame(GeneXData)
row.names(GeneXData) <- gene.IDs
GeneXData$GeneName <- NULL
GeneXData$GeneID <- NULL
GeneXData <- as.data.frame(t(GeneXData))
colnames(GeneXData) <- gene.IDs


################################################################################
# Remove to_remove_genes
################################################################################
cols <- setdiff(colnames(GeneXData), to_remove_genes)
GeneXData <- GeneXData[, cols]


################################################################################
# Set up the parallel calculation back-end
################################################################################
enableWGCNAThreads()


################################################################################
# Apply the soft threshold function to find the optimal power for the adjacency 
# matrix
################################################################################
sft = pickSoftThreshold(GeneXData)


################################################################################
# Scatter plot for the model fitting of scale free topology in function with of 
# the adjacency power
################################################################################
p <- ggplot(data = sft$fitIndices, aes(x=Power, y=-sign(slope)*SFT.R.sq)) +
      geom_point(size=2, shape=4, col="#863434") +
      geom_line(color="#863434")+
      geom_hline(yintercept=0.9, linewidth = 1, col="#863434") + 
      scale_x_continuous(breaks = seq(2, max(sft$fitIndices$Power), 2)) +
      labs(x='Power of Adjacency Function', y=bquote('Scale Free Topology Model Fit '(R^2))) +
      theme(axis.text = element_text(size=10),
            axis.title = element_text(size=12))
ggsave(paste(output_dir, 'soft_threshold_fitting.png', sep=''), p, device="png", 
       dpi=600, height=5, width=8.3, units=c("in"))
ggsave(paste(output_dir, 'soft_threshold_fitting.tiff', sep=''), p, device="tiff", 
       dpi=300, height=5, width=8.3, units=c("in"))


################################################################################
# Scatter plot for the mean connectivity in function with of the adjacency power
################################################################################
p <- ggplot(data = sft$fitIndices, aes(x=Power, y=mean.k.)) +
  geom_point(size=2, shape=4, col="#863434") +
  geom_line(color="#863434") +
  scale_x_continuous(breaks = seq(2, max(sft$fitIndices$Power), 2)) +
  labs(x='Power of Adjacency Function', y=bquote('Mean Connectivity')) +
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=12))
ggsave(paste(output_dir, 'mean_connectivity.png', sep=''), p, device="png",
       dpi=600, height=5, width=8.3, units=c("in"))
ggsave(paste(output_dir, 'mean_connectivity.tiff', sep=''), p, device="tiff",
       dpi=300, height=5, width=8.3, units=c("in"))


################################################################################
# Save the results of soft threshold function
################################################################################
filename <- paste(output_dir, 'soft_threshold_results.tsv', sep='')
write.table(x = sft$fitIndices, file = filename, sep = '\t', quote = FALSE, 
            row.names = FALSE)


################################################################################
# Save the objects
################################################################################
filename <- paste(output_dir, "soft_threshold_results.RData", sep='')
save(GeneXData, sft, ensembl_mapping, file=filename)
#print(paste('power estimation:', sft$powerEstimate))
