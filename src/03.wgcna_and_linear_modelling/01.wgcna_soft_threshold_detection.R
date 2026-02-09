suppressMessages(library(readr)) # 2.1.6
suppressMessages(library(WGCNA)) # 1.73
suppressMessages(library(ggplot2)) # 4.0.1

################################################################################ 
# Description
################################################################################
# This is the first step of WGCNA analysis to detect the optimal 
# soft-threshold power for the correlation-based network construction. The
# optimal power is the minimum value which transforms the adjacency matrix
#into a network with scale free topology.
# Outputs:
# - Two plots to illustrate the main results of the analysis.
# - The main results of soft-threshold power detection in tsv and RData formats.
################################################################################

################################################################################ 
# Inputs
################################################################################
# Batch corrected expression dataset is the main input
counts_filename <- '../data/epos/batch_corrected_counts_matrix.csv'
# tsv file which contains the mapping of ensembl ids and gene symbols
ensembl_mapping <- read.table('../../data/ensembl_mapping.tsv', sep='\t', 
                              header=TRUE, row.names=1)
# Output directory to save the results of soft-threshold detection analysis
output_dir <- '../../results/epos/wgcna_and_linear_modelling/'
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

################################################################################ 
# 1. Load ensembl mapping and selecte the unannotated genes in order to remove them
# from the dataset
################################################################################
to_remove_genes <- ensembl_mapping[ensembl_mapping$external_gene_name == "",'ensembl_gene_id']

################################################################################
# 2. Load the expression dataset and format it properly
################################################################################
GeneXData <- read_csv(counts_filename)
colnames(GeneXData)[1] <- "GeneID"
gene.IDs <- GeneXData$GeneID
GeneXData <- as.data.frame(GeneXData)
row.names(GeneXData) <- gene.IDs
GeneXData$GeneName <- NULL
GeneXData$GeneID <- NULL
GeneXData <- as.data.frame(t(GeneXData))
colnames(GeneXData) <- gene.IDs

################################################################################
# 3. Remove the unannotated genes (to_remove_genes)
################################################################################
cols <- setdiff(colnames(GeneXData), to_remove_genes)
GeneXData <- GeneXData[, cols]

################################################################################
# 4. Set up the parallel calculation back-end
################################################################################
enableWGCNAThreads()

################################################################################
# 5. Apply the soft threshold function to find the optimal power for the 
# adjacency matrix
################################################################################
sft = pickSoftThreshold(GeneXData)

################################################################################
# 6. Scatter plot for the model fitting of scale free topology in function with 
# of the adjacency power
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
# 7. Scatter plot for the mean connectivity in function with of the adjacency 
# power
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
# 8. Save the results of soft threshold function in tsv and RData formats
################################################################################
filename <- paste(output_dir, 'soft_threshold_results.tsv', sep='')
write.table(x = sft$fitIndices, file = filename, sep = '\t', quote = FALSE, 
            row.names = FALSE)

filename <- paste(output_dir, "soft_threshold_results.RData", sep='')
save(GeneXData, sft, ensembl_mapping, file=filename)
#print(paste('power estimation:', sft$powerEstimate))
